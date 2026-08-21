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
#' @param edge A character string representing the type of edge effect compensation
#' to apply. The allowed values are "none" for no compensation, "weight" to
#' increase the weight of points near xmin and xmax, and "periodic" to treat the
#' values as if they repeat periodically
#' @param nstep An integer value used to automate the creation of vector of t values
#'
#' @export
#'
ripley_k <- function(
    x,
    t = NA,
    xmin = NA,
    xmax = NA,
    edge,
    nstep = NA
) {
    # Get length of input vector and apply sanity check
    n <- length(x)
    if (n < 2) {
        stop("Too few values to analyse ", n)
    }

    if (!edge %in% c("weight", "periodic", "none")) {
        stop(
            "Invalid edge method ",
            edge,
            " - use none, weight, or periodic"
        )
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

    # For periodic edge calculate new dm
    if (edge == "periodic") {
        range <- xmax - xmin
        lower <- x < range / 2
        x_periodic <- x
        x_periodic[lower] <- x_periodic[lower] + range
        dm_periodic <- dist(x_periodic)
        dm <- pmin(dm, dm_periodic)
    }

    # Setup output data frame
    output <- data.frame(t = t, K = NA, L = NA)

    # Look at all values for t
    for (ix in 1:nrow(output)) {
        tt <- output$t[ix]
        dm_t <- as.dist(matrix(
            as.numeric(as.matrix(dm) < tt),
            ncol = n
        ))
        if (edge == "weight") {
            w1 <- 1 - pmin(x - xmin - tt, 0) / tt
            w2 <- 1 - pmin(xmax - x - tt, 0) / tt
            w <- pmax(w1, w2)
            # Combine weights to form a matrix
            wm <- matrix(nrow = n, ncol = n)
            for (i in 1:(n - 1)) {
                for (j in i:n) {
                    wm[j, i] <- max(w[i], w[j])
                }
            }
            wm <- as.dist(wm)
            K <- sum(dm_t * wm) / n / lindensity
        } else {
            K <- sum(dm_t) / n / lindensity
        }

        output$K[ix] <- K
        output$L[ix] <- sqrt(K / pi)
    }
    output
}
