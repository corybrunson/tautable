#' Count the inversions in a vector.
#' 
#' @param vec A numerical vector (theoretically, with no repeated values).
#' @return The number of pairs of entries *i* and *j* of \code{vec} for which *i* > *j*.
#' @examples
#' vec.inv(1:5)
#' vec.inv(5:1)
vec.inv <- function(vec) {
    sum(unlist(sapply(1:(length(vec) - 1), function(i) {
        sapply((i + 1):length(vec), function(j) vec[i] > vec[j])
    })))
}
