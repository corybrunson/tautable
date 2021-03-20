#' Compute the Mahonian numbers.
#'
#' @param n A positive integer.
#' @return The number of permutations of `n` having `k` inversions, as `k`
#'   ranges from 0 to `n`.
#' @export
#' @examples
#' mahonians(1)
#' mahonians(5)
mahonians <- function(n) {
    if (n <= 0 | n %% 1 != 0) stop("`n` must be a positive integer.")
    vec <- 1; m <- 1
    while (n > m) {
        vec <- rowSums(sapply(0:m, function(i) {
            c(rep(0, i), vec, rep(0, m - i))
        }))
        m <- m + 1
    }
    stopifnot(length(vec) == choose(n, 2) + 1)
    return(vec)
}
