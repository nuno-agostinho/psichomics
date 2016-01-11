#' Creates a data frame of parsed junctions filled with NAs
#' 
#' @param nrow Integer: Number of rows
#' @param program Character: Program used to get the junctions
#' @param event.type Character: Event type of the respective events
#' @param chromosome Character: Chromosome of the junctions
#' @param strand Strand: positive ("+") or negative ("-") strand of the event
#' 
#' @return A data frame with the junctions coordinate names pre-filled with NAs
#' @export
#' 
#' @examples
#' createFilledJunctions(nrow = 8)
createFilledJunctions <- function(nrow, program = character(0),
                                  event.type = character(0),
                                  chromosome = character(0),
                                  strand = character(0)) {
    parsed <- as.data.frame(matrix(NA, nrow = nrow, ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")
    
    if (length(program) > 0)    parsed[["Program"]] <- "MISO"
    if (length(event.type) > 0) parsed[["Event.type"]] <- event.type
    if (length(chromosome) > 0) parsed[["Chromosome"]] <- chromosome
    if (length(strand) > 0)     parsed[["Strand"]] <- strand
    return(parsed)
}

#' Converts a data frame of junctions to a list of numeric junctions
#'
#' @param junctions Data frame of junctions
#' @param rows Rows of the junctions
#' @param cols Columns of the junctions
#'
#' @return List of numeric junctions
#' @export
listJunctions <- function(junctions, rows, cols) {
    apply(junctions[rows, cols], 1, function(i) as.list(as.numeric(i)))
}