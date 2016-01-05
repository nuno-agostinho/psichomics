#' Creates a data frame of parsed junctions filled with NAs
#' 
#' @param nrow Integer: Number of rows
#' 
#' @return A data frame with the junctions coordinate names pre-filled with NAs
#' @export
#' 
#' @examples
#' createFilledJunctions(nrow = 8)
createFilledJunctions <- function(nrow) {
    parsed <- as.data.frame(matrix(NA, nrow = nrow, ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")
    return(parsed)
}


#' Converts a vector of numeric characters to a list of vectors of integers
#'
#' @param i Numeric character vector
#'
#' @return List of numeric vectors
#' @export
numericList <- function(i) {
    return(as.list(as.numeric(i)))
}