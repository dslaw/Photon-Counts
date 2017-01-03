## Hierarchical model for photon count data
## Modeling the photon counts from stars in different galaxies
## where the photon count from each star is described by a Poisson
## distribution with parameter lambda_(j)
## y_ij = photon count for ith star in jth galaxy
## n_j = number of stars in jth galaxy
## J = number of galaxies
## lambda_j = mean photon count in jth galaxy (summing over i)

require(compiler); enableJIT(3)
require(MCMCpack)
require(ggplot2)

source("adaptive.R")

## Load and check data
data = read.table("stars.txt", header = TRUE)
J = length(unique(data[["galaxy"]]))

if (J != max(data[["galaxy"]]))
    is.missing = !(1:max(data[["galaxy"]]) %in% unique(data[["galaxy"]]))
    cat(which(is.missing), "\n")

data[data[["galaxy"]] < 0, "galaxy"] = NA

## Useful statistics
stars = split(data, data$galaxy)
j = as.numeric(names(stars))
n = unname(sapply(stars, nrow)) # sample size in jth galaxy
N = nrow(data) # total sample size
yj = lapply(stars, function(galaxy) galaxy[["photon.count"]])
lambda_j = sapply(stars, function(galaxy) mean(galaxy[["photon.count"]]))

## Plot data
gg_data = data.frame(Galaxy = factor(data[["galaxy"]]),
                     Photon.Count = data[["photon.count"]])

p <- ggplot(gg_data, aes(x = Galaxy, y = Photon.Count))
p + geom_boxplot() +
    xlab("Galaxy") +
    ylab("Photon count") +
    ggtitle("Distribution of photon counts by galaxy")

## Metropolis within Gibbs sampler
# Analyst defined prior parameters
# These parameters correspond to non-informative prior distributions
# alpha is the shape parameter and beta is the rate parameter
a.alpha = 0.001
a.beta = 0.001
b.alpha = 0.001
b.beta = 0.001

## Using Adaptive MWG
nsamples = 15000L

# Initialize parameter vector theta
# Note that all parameters must be greater than 0
theta = list(alpha = 1,
             beta = 1,
             lambda = rep(1, J))
nu = 2 # Std. deviation of proposal distribution (normal)

set.seed(13)
results = AdaptiveMWG(theta = theta, nu = nu,
                      n.samples = nsamples, periodicty = 100L,
                      rate = 1, bounds = c(0.4, 0.5),
                      J = J, a.alpha = a.alpha, a.beta = a.beta,
                      b.alpha = b.alpha, b.beta = b.beta,
                      lambda_j = lambda_j, yj = yj, n = n)
n_accepted = results$number.accepted
nu = results$nu
samples = results$samples; rm(results)
colnames(samples) = c("alpha", "beta", paste0("lambda", 1L:J))
cat(n_accepted / nsamples, "\n")

## Check that all samples are positive
if (any(samples <= 0))
    warning("Non-positive samples")

## Check trace plots for convergence and to determine burn in period
plot(mcmc(samples))
# No burn in period is necessary

## Diagnostic plots:
# Effective sample size
# Want ~1000 draws => effective sample size of 1000
# If the effective sample size is too low, nu should be tuned further or the
# number of samples taken should be increased
ess = effectiveSize(samples)

# Check for dependence
acf(samples[, c("alpha", "beta")])

## Compute summary statistics of posterior draws
summary_stats = apply(samples, 2, function(draws) {
                      summry = c(summary(draws), sd(draws));
                      names(summry)[length(summry)] = "Sd";
                      summry})
summary_stats = t(summary_stats)

## Compute Posterior Predictive Distributions (ppd)
# Use Monte Carlo simulations to sample from the new likelihood p(y_{~} | y)
# Draw m sets of samples from the posterior p(theta | y) (done with MWG)
# Use each set of sample parameters to simulate a new posterior predictive
# dataset that matches the dimensions of the original dataset

# Only need lambda samples to simulate the data
lambda_draws = samples[, grepl("lambda", colnames(samples))]
m = nrow(lambda_draws)

# Allocate storage: list of m matrices, each matrix with the same dimensions
# as the original dataset
ppdset = matrix(NA, nrow = N, ncol = ncol(data))
ppdset[, 1] = data[["galaxy"]]
colnames(ppdset) = colnames(data)
ppd_samples = list()

for (i in seq_len(m)) {

    ppd = lapply(seq_len(J), function(j)
                 rpois(n = n[j], lambda = lambda_draws[i, j]))
    ppdset[, 2] = unlist(ppd)
    ppd_samples[[i]] = ppdset
}

## Summary statistics of each PPD dataset
summary_ <- function(x){
    s = c(min(x), median(x), mean(x), max(x), sd(x))
    names(s) = c("Min", "Median", "Mean", "Max", "Sd")
    return(s)
}

ppd_summary_stats = sapply(ppd_samples, function(x)
                           summary_(x[, "photon.count"]))
ppd_summary_stats = t(ppd_summary_stats)

## Compare observed summary statistics to the distribution of PPD
## summary statistics
data_summary = summary_(data[["photon.count"]])

gg_summary_stats = data.frame(value = c(ppd_summary_stats, recursive = TRUE),
                              variable = c(rep("Min", m), rep("Median", m),
                                           rep("Mean", m), rep("Max", m),
                                           rep("Sd", m)))
gg_data_summary = data.frame(observed = data_summary,
                             variable = names(data_summary))

ggplot(gg_summary_stats, aes(x = value)) +
facet_wrap(~ variable, scales = "free_x") +
geom_histogram(aes(colour = "ppd")) +
ylab("Count") +
xlab("") +
ggtitle("Distributions of summary statistics\nfrom ppd datasets") +
geom_vline(data = gg_data_summary,
           aes(xintercept = observed, colour = "Observed")) +
scale_colour_manual(name = "", values = c("ppd" = NA, "Observed" = "blue"))

## Estimated model parameters
alpha_estimate = median(samples[, "alpha"])
beta_estimate = median(samples[, "beta"])
lambda_estimates = apply(lambda_draws, 2, median)

# 95% Credible Intervals
credible_intervals = apply(samples, 2, function(draws)
                           quantile(draws, probs = c(0.25, 0.975)))

ggdf = data.frame(Galaxy = rep(j, 2),
                  Lambda = c(lambda_estimates, lambda_j),
                  Group = c(rep("Estimated", J), rep("Observed", J)))

p <- ggplot(ggdf, aes(x = Galaxy, y = Lambda, group = Group, color = Group))
p + geom_point() + xlab("Galaxy") + ylab(substitute(lambda[j])) +
ggtitle("Model parameters")

