#' Rapidly generate a table of critical values of tau or an equivalent test
#' statistic.
#' 
#' @param n A vector of positive integers.
#' @param alpha A vector of probabilities (preferably below .5).
#' @param incl.len Logical; whether to include a column of values of \code{n}
#'   (default to \code{FALSE}; the rows are named by the values of \code{n} in
#'   any case).
#' @param stat A character string; any of several recognized test statistics
#'   (\code{"tau"}, \code{"concord*"}, , \code{"discord*"}, \code{"P"},
#'   \code{"Q"}, \code{"K"}, \code{"k"}, \code{"inv"}, \code{"S"}).
#' @return The number of permutations of \code{n} having \code{k} inversions, as
#'   \code{k} ranges from 0 to \code{n}.
#' @export
#' @examples
#' tau_crit_table(n = 1:5, alpha = 10 ^ (-1:-3), incl.len = TRUE, stat = "inv")
#' tau_crit_table(n = 1:5, alpha = 10 ^ (-1:-3), incl.len = TRUE, stat = "tau")
tau_crit_table <- function(n, alpha, incl.len = TRUE, stat = "tau") {
    if(any(n <= 0 | n %% 1 != 0)) stop("n must be a positive integer")
    if(any(alpha > .5 | alpha <= 0)) stop("alpha must be positive and below .5")
    vec <- 1; m <- 1; max.n <- max(n)
    mat <- matrix(NA, nr = if(min(n) == 1) 1 else 0, nc = length(alpha))
    while(max.n > m) {
        # Update Mahonian vector
        vec <- rowSums(sapply(0:m, function(i) {
            c(rep(0, i), vec, rep(0, m - i))
        }))
        # Update Mahonian index
        m <- m + 1
        # If m has a value from n...
        if(m %in% n) {
            # Compute proportions of total sum at first half of partial sums
            s <- cumsum(vec[1:floor(choose(m, 2) / 2)]) / sum(vec)
            # Identify index (if any) at which proportion drops below alpha
            i <- sapply(alpha,
                        function(a) if(s[1] < a) max(which(s < a)) - 1 else NA)
            # Record these indices (numbers of inversons) in matrix
            mat <- rbind(mat, i)
        }
    }
    stopifnot(length(vec) == choose(max.n, 2) + 1)
    # Data frame
    colnames(mat) <- alpha
    rownames(mat) <- n
    # Transform to desired statistic
    if(stat == "tau") {
        mat <- (choose(n, 2) - 2 * mat) / choose(n, 2)
    } else if(stat == "P" | grepl("^concord", stat, ignore.case = TRUE)) {
        mat <- choose(n, 2) - mat
    } else if(stat == "S") {
        mat <- choose(n, 2) - 2 * mat
    } else if(!(stat %in% c("K", "k", "Q")) &
                  !grepl("^discord|^inv", stat, ignore.case = TRUE)) {
        stop("unknown statistic")
    }
    # Add length column if desired
    if(incl.len) mat <- cbind(n = n, mat)
    return(mat)
}