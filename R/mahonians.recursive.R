#' Compute the Mahonian numbers recursively (very slow for n > 8).
#' 
#' @param n A positive integer.
#' @return The number of permutations of \code{n} having \code{k} inversions, as \code{k} ranges from 0 to \code{n}.
#' @export
#' @examples
#' mahonians.recursive(1)
#' mahonians.recursive(5)
mahonians.recursive <- function(n) {
    if(n <= 0 | n%%1 != 0) stop('n must be a positive integer')
    if(n == 1) 1 else rowSums(sapply(0:(n - 1), function(i) {
        c(rep(0, i), mahonians.recursive(n - 1), rep(0, n - 1 - i))
    }))
}