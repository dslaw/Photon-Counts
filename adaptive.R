# Adaptive Metropolis within Gibbs sampler.

adapt <- function(nu, j, acceptance.rate, rate,
                  bounds = c(.35, .5)){
    ## Adapt tuning paramter nu

    if (acceptance.rate < min(bounds)) {
        s = -1
    } else if (acceptance.rate > max(bounds)) {
        s = 1
    } else {
        s = 0
    }

    # Adjust nu based on acceptance rate
    # scale adjustment by jth adaptation to satisfy diminishing adaptation
    # property
    nu = exp(rate * log(nu) + s / j)
    return(nu)
}

alpha_conditional <- function(alpha, beta, J, lambda_j, a.alpha, b.alpha){
    ## Conditional distribution of alpha, log scale (no closed form)

    J * (alpha * log(beta) - lgamma(alpha)) +
    (alpha - 1) * sum(log(lambda_j)) +
    (a.alpha - 1) * log(alpha)
}

adaptive_mwg <- function(theta, nu, n.samples, periodicity = 50L,
                         rate = 1, bounds = c(0.4, 0.5), ...){
    ## Adaptive Metropolis within Gibbs

    if (nu < 0)
        stop("Tuning parameter must be positive")

    alpha = theta[["alpha"]]
    beta = theta[["beta"]]
    lambda = theta[["lambda"]]

    samples = matrix(NA, nrow = n.samples, ncol = length(unlist(theta)))
    accepted = logical(n.samples)

    for (i in seq_len(n.samples)) {
        # Sample from the conditional of beta
        beta = rgamma(n = 1,
                      shape = (J*alpha) + a.beta,
                      rate = sum(lambda_j) + b.beta)

        # Adaptive Metropolis sampler for alpha
        # using a truncated normal (alpha > 0) proposal distribution
        if (i %% periodicity == 0) {
            j = i / periodicity
            nu = adapt(nu = nu,
                       j = j,
                       acceptance.rate = sum(accepted) / i,
                       rate = rate,
                       bounds = bounds)
        }

        while (TRUE) {
            proposal = rnorm(n = 1, mean = alpha, sd = nu)
            if (proposal > 0) break
        }

        u = runif(n = 1, 0, 1)
        r = alpha_conditional(alpha = proposal, beta = beta,
                              J = J, lambda_j = lambda_j,
                              a.alpha = a.alpha, b.alpha = b.alpha) -
            alpha_conditional(alpha = alpha, beta = beta,
                              J = J, lambda_j = lambda_j,
                              a.alpha = a.alpha, b.alpha = b.alpha)
        if (log(u) <= r) {
            accepted[i] = TRUE
            alpha = proposal
        }

        # Sample from conditional of lambda_j for all j
        lambda = sapply(seq_along(lambda), function(j)
                        rgamma(n = 1,
                               shape = sum(yj[[j]] + alpha),
                               rate = beta + n[j]))

        samples[i, ] = c(alpha, beta, lambda)
    }

    list(number.accepted = sum(accepted),
         nu = nu,
         samples = samples)
}

