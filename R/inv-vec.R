#' Generate a permutation of a given length having a given number of inversions
#'
#' @param n A positive integer.
#' @param k A nonnegative integer at most `n`.
#' @param allow.dual Logical; whether to allow exchanging large values of `k`
#'   for their duals to improve efficiency.
#' @return A permutation of `n` having `k` (non-arbitrary) inversions.
#' @export
#' @examples
#' inv_vec(n = 5, k = choose(5, 2))
#' inv_vec(n = 5, k = 7)
#' inv_vec(n = 5, k = 7, allow.dual = FALSE)
inv_vec <- function(n, k, allow.dual = TRUE) {
    if(k > choose(n, 2) | k < 0 | any(c(n, k) %% 1 != 0))
      stop("Invalid values of `n` and `k`.")
    backwards <- (k > choose(n, 2) / 2) & allow.dual
    if(backwards) k <- choose(n, 2) - k
    vec <- 1:n; ninv <- 0
    while(ninv < k) {
        i <- vec[1]
        j <- if(k - ninv < n - i) (2 + k - ninv):n else {
            if(i == 1) c() else (2 + n - i):n
        }
        vec <- c(vec[-c(1, j)], i, vec[j])
        ninv <- ninv + {
            if(k - ninv < n - i) (k - ninv) else (n - i)
        }
        stopifnot(ninv == vec_inv(vec))
    }
    if(backwards) vec <- rev(vec)
    return(vec)
}
