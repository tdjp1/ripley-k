#' Statistics for Ripley K function in 1D
#'
#' This takes a dataset, calculates K and makes a comparison
#' with a comparable set of Monte Carlo simulations to
#' estimate significance
#'
#' @param x vector of values for analysis
#' @param t vector of windows
#' @param xmin Minimum possible for data x. If xmin is set to NA then the value is
#' calculated from x
#' @param xmax Maximum possible value for data x. If xmax i set to NA then the value is
#' calculcated from x
#' @param edge A character string representing the type of edge effect compensation
#' to apply. The allowed values are "none" for no compensation, "weight" to
#' increase the weight of points near xmin and xmax, and "periodic" to treat the
#' values as if they repeat periodically
#' @param nstep An integer value used to automate the creation of vector of t values
#'
#' @export
#'
ripley_k_stat <- function(
    x,
    t = NA,
    xmin = NA,
    xmax = NA,
    edge = "weight",
    nstep = NA,
    nsim = 1000
) {
    # Get length of input vector and apply sanity check
    n <- length(x)
    if (n < 2) {
        stop("Too few values to analyse ", n)
    }

    # Auto create xmin and xmax if not supplied
    if (is.na(xmin)) {
        xmin <- min(x)
    } else {
        if (min(x) < xmin) {
            stop("Data values lower than expected value ", xmin)
        }
    }
    if (is.na(xmax)) {
        xmax <- max(x)
    } else {
        if (max(x) > xmax) {
            stop("Data values higher than expected value ", xmax)
        }
    }

    # Create or check on t vector
    if (!is.na(nstep)) {
        t <- seq(
            from = (xmax - xmin) / (2 * nstep),
            to = (xmax - xmin) / 2,
            length.out = nstep
        )
    } else {
        if (anyNA(t)) stop("Invalid vector of values for t")
    }

    out <- ripley_k(x, t = t, xmin = xmin, xmax = xmax, edge = edge)
    score <- ripley_k_score(out)

    # Run simulation
    sim <- ripley_k_mc(
        n,
        t,
        xmin = xmin,
        xmax = xmax,
        edge = edge,
        nsim = nsim
    )
    q <- quantile(sim, c(0.05, 0.95))

    if (score < q[1]) {
        ret <- list(score, "diverse")
    } else if (score > q[2]) {
        ret <- list(score, "clustered")
    } else {
        ret <- list(score, "random")
    }
    ret
}
