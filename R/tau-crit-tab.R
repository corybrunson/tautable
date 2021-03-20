#' Tabulate critical values for tau and other test statistics via backward use
#' of correlation test functions.
#'
#' @param n A vector of positive integers.
#' @param alpha A vector of probabilities (preferably below .5).
#' @return For each pair of values from `n` and `alpha`, the critical value of
#'   the desired test statistic for rankings of `1:n`.
#' @param ... Additional parameters passed to [tau_a_crit()].
#' @export
#' @examples
#' tau_crit_tab(n = seq(10, 30, 10), alpha = 10 ^ (-1:-2))
tau_crit_tab <- function(n, alpha, ...) {
    mat <- sapply(n, function(nn) tau_a_crit(nn, alpha, ...))
    df <- as.data.frame(t(mat))
    names(df) <- alpha
    rownames(df) <- n
    return(df)
}
