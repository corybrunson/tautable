#' Compute the Mahonian numbers recursively (very slow for n > 8).
#'
#' @param n A positive integer.
#' @return The number of permutations of `n` having `k` inversions, as `k`
#'   ranges from 0 to `n`.
#' @export
#' @examples
#' mahonians_recursive(1)
#' mahonians_recursive(5)
mahonians_recursive <- function(n) {
    if(n <= 0 | n%%1 != 0) stop("`n` must be a positive integer")
    if(n == 1) 1 else rowSums(sapply(0:(n - 1), function(i) {
        c(rep(0, i), mahonians_recursive(n - 1), rep(0, n - 1 - i))
    }))
}
