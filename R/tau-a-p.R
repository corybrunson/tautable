#' Tabulate tau and p-values for a given length and number of inversions.
#'
#' @importFrom Kendall Kendall
#' @param n A positive integer.
#' @param k A vector of nonnegative integers at most `n`.
#' @param sided Number of sides to a test; either 1 or 0.
#' @param tau.fn A function that computes elements tau and p-value from a pair
#'   of vectors.
#' @param tau.part The name of the tau element of the value of tau.fn.
#' @param p.part The name of the p-value element of the value of tau.fn.
#' @param ... Additional parameters passed to `tau.fn`.
#' @return For each entry of k, the corresponding value of tau and associated
#'   p-value.
#' @export
#' @examples
#' tau_a_p(5, 0:5, sided = 1)
#' tau_a_p(5, 5:10, sided = 1)
#' tau_a_p(5, 0:5, tau.fn = Kendall::Kendall, tau.part = "tau", p.part = "sl")
tau_a_p <- function(
    n, k = floor((choose(n, 2) / 2)):0, sided = 1,
    tau.fn = stats::cor.test, tau.part = "estimate", p.part = "p.value", ...
) {
    # All tests are two-sided; need to adjust the results for one-sided pvals
    stopifnot(sided %in% 1:2)
    mat <- sapply(k, function(kk) {
        test <- tau.fn(x = 1:n, y = inv_vec(n, kk), ...)
        wh.tau.p <- sapply(c(tau.part, p.part), function(part) {
            which(names(test) == part)
        })
        tau.p <- unlist(unname(test[wh.tau.p]))
        return(c(k = kk, K = (tau.p[1] + 1) * n * (n - 1) / 4,
                 S = tau.p[1] * choose(n, 2), tau = tau.p[1],
                 p = tau.p[2] / (3 - sided)))
    })
    return(as.data.frame(t(mat)))
}
