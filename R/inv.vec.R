#' Generate a permutation of a given length having a given number of inversions
#' 
#' @param n A positive integer.
#' @param k A nonnegative integer at most \code{n}.
#' @param bw.ok Logical; whether to allow exchanging large values of \code{k} for their duals to improve efficiency.
#' @return A permutation of \code{n} having \code{k} (non-arbitrary) inversions.
#' @export
#' @examples
#' inv.vec(n = 5, k = choose(5, 2))
#' inv.vec(n = 5, k = 7)
#' inv.vec(n = 5, k = 7, bw.ok = FALSE)
inv.vec <- function(n, k, bw.ok = TRUE) {
    if(k > choose(n, 2) | k < 0 | any(c(n, k) %% 1 != 0)) stop('invalid n or k')
    backwards <- (k > choose(n, 2) / 2) & bw.ok
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
        stopifnot(ninv == vec.inv(vec))
    }
    if(backwards) vec <- rev(vec)
    return(vec)
}