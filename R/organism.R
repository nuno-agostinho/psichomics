#' Organism representation
#'
#' @slot species Character: species name
#' @slot common.name Character: common name
#' @slot clinical.information Data.frame: clinical information
#' @slot inclusion.levels Data.frame: inclusion levels of a given exon
#' @export
#'
#' @examples
#' # create a new organism object ('Homo sapiens' by default)
#' org <- new("Organism") 
#' 
#' # create a new organism object (this time a mouse)
#' org <- new("Organism", species = "Mus musculus", common.name = "Mouse")
#' 
#' # string representation of the object (displays the species and respective
#' # common name of a organism)
#' as(org, "character")
setClass("Organism",
         representation = representation(
             species = "character",
             common.name = "character",
             clinical.information = "data.frame",
             inclusion.levels = "data.frame"),
         prototype = prototype(
             species = "Homo sapiens",
             common.name = "Human"))

# Checks validity of Organism's attributes
setValidity("Organism",
            function(object) {
                if (length(object@species) != 1) 
                    "'species' should be a single string"
                else if (length(object@common.name) != 1) 
                    "'common.name' should be a single string"
                else if (!nzchar(object@species))
                    "'species' is empty"
                else
                    TRUE
            })

# String representation of the object
setAs("Organism", "character", function(from)
    if (from@common.name != "")
        sprintf("%s (%s)", from@species, from@common.name)
    else
        sprintf("%s", from@species))


#' Gets the inclusion levels of an organism
#'
#' @param object Organism
#'
#' @return Inclusion levels of an organism
#' @export
#'
#' @examples
#' org <- new("Organism",
#'            inclusion.levels = data.frame(1,2,3))
#' inclusion.levels(org)
setGeneric("inclusion.levels", function(object)
    standardGeneric("inclusion.levels"))

#' @rdname inclusion.levels
setMethod("inclusion.levels", "Organism", function(object)
    object@inclusion.levels)

#' Sets the inclusion levels of an organism
#'
#' @param object Organism
#' @param value Data.frame
#' @export
#'
#' @examples
#' org <- new("Organism")
#' inclusion.levels(org) <- data.frame(1,2,3)
setGeneric("inclusion.levels<-", function(object, value)
    standardGeneric("inclusion.levels<-"))

#' @rdname inclusion.levels-set
setReplaceMethod("inclusion.levels", "Organism",
                 function(object, value)
                     initialize(object, inclusion.levels = value))