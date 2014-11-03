#' Tabulate critical values for tau and other test statistics.
#' 
#' @param n A vector of positive integers.
#' @param alpha A vector of probabilities (preferably below .5).
#' @return For each pair of values from \code{n} and \code{alpha}, the critical value of the desired test statistic for rankings of \code{1:n} .
#' @examples
#' tau.crit.tab(n = seq(10, 30, 10), alpha = 10 ^ (-1:-2))
tau.crit.tab <- function(n, alpha, ...) {
    mat <- sapply(n, function(nn) tau.a.crit(nn, alpha, ...))
    df <- as.data.frame(t(mat))
    names(df) <- alpha
    rownames(df) <- n
    return(df)
}
