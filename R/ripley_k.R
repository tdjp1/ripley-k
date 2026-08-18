#' Ripley K function for 1D data
#'
#' Uses Ripley K function to analyse whether data are evenly dispersed, random or clustered
#'
#' @param x vector of values for analysis
#' @param t vector of windows
#' @param xmin Minimum possible for data x. If xmin is set to NA then the value is
#' calculated from x
#' @param xmax Maximum possible value for data x. If xmax i set to NA then the value is
#' calculcated from x
#' @param weighting A character string representing the type of weighting to apply
#' to mitigate against edge effects. Allowed values are "none" for no weighting, "reflect" to
#' increase the weight of points near xmin and xmax. "cycle" should added in thr future.
#' @param nstep An integer value used to automate the creation of vector of t values
#'
#' @export
#'
ripley_k <- function(
    x,
    t = NA,
    xmin = NA,
    xmax = NA,
    weighting = "reflect",
    nstep = NA
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
    # Get linear density of points
    lindensity <- n / (xmax - xmin + 1)

    # Generate distance matrix
    dm <- dist(x)

    # Setup output data frame
    output <- data.frame(t = t, K = NA, L = NA)

    # Look at all values for t
    for (ix in 1:nrow(output)) {
        tt <- output$t[ix]
        if (weighting == "reflect") {
            w1 <- 1 - pmin(x - xmin - tt, 0) / tt
            w2 <- 1 - pmin(xmax - x - tt, 0) / tt
            w <- pmax(w1, w2)
        } else if (weighting == "none") {
            w <- rep(1, n)
        } else {
            stop("Invalid weighting mode ", weighting)
        }

        # Combine weights to form a matrix
        wm <- matrix(nrow = n, ncol = n)
        for (i in 1:(n - 1)) {
            for (j in i:n) {
                wm[j, i] <- max(w[i], w[j])
            }
        }
        wm <- as.dist(wm)
        dm_t <- as.dist(matrix(
            as.numeric(as.matrix(dm) < tt),
            ncol = n
        ))
        K <- sum(dm_t * wm) / n / lindensity
        output$K[ix] <- K
        output$L[ix] <- sqrt(K / pi)
    }
    output
}
