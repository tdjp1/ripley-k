#' Monte Carlo calculation for Ripley K function for 1D data
#'
#' Uses Ripley K function to analyse whether data are evenly dispersed,
#' random or clustered @param n number of values for analysis @param t vector of windows
#' @param xmin Minimum possible for data x.
#' @param xmax Maximum possible value for data x.
#' @param weighting A character string representing the type of weighting to apply
#' to mitigate against edge effects. Allowed values are "none" for no weighting, "reflect" to
#' increase the weight of points near xmin and xmax. "cycle" should added in thr future.
#' @param nsim Number of simulations to run
#'
#' @export
#'
ripley_k_mc <- function(
    n,
    t,
    xmin,
    xmax,
    weighting = "reflect",
    nsim = 100
) {
    # Some sanity checks
    if (n < 2) {
        stop("Too few values to analyse ", n)
    }

    if (anyNA(t)) {
        stop("Invalid vector of values for t")
    }

    # Start calculations
    sapply(1:nsim, function(i) {
        x <- as.integer(round(runif(n, min = xmin, max = xmax), 0))
        results <- ripley_k(
            x,
            t = t,
            xmin = xmin,
            xmax = xmax,
            weighting = weighting
        )
        diff <- sum(results$K - results$t)
    })
}
