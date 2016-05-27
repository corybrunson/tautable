#' Find the critical test statistic for a given number of objects and type-I error.
#' 
#' @param n A positive integer.
#' @param alpha A vector of probabilities (preferably below .5).
#' @param stat A character string; the name of the desired element returned by tau.a.p.
#' @return For each entry of alpha, the least extreme test statistic whose p-value is smaller.
#' @export
#' @examples
#' tau.a.crit(5, 10 ^ (-1:-3), stat = 'k')
tau.a.crit <- function(n, alpha, stat = 'tau', ...) {
    tab <- tau.a.p(n, ...)
    i <- which(names(tab) == stat)
    mat <- sapply(alpha, function(a) {
        tab[which(tab$p < a)[1], i]
    })
    return(mat)
}
