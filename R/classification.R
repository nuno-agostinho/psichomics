#' Classification representation
#'
#' @slot species Character: species name
#' @slot common.name Character: common name
#' @slot condition Character: condition type of the group
#' @slot date Date: time stamp of the samples
#' @slot clinical Data.frame: clinical information
#' @slot inclusion.levels Data.frame: inclusion levels of a given exon
#' @export
#'
#' @examples
#' # create a new classification object ('Homo sapiens' by default)
#' org <- new("Classification") 
#' 
#' # create a new classification object (this time a mouse)
#' org <- new("Classification", species = "Mus musculus", common.name = "Mouse")
#' 
#' # string representation of the object (displays the species and respective
#' # common name of a classification)
#' as(org, "character")
setClass("Classification",
         representation = representation(
             species = "character",
             common.name = "character",
             condition = "character",
             date = "Date",
             clinical = "data.frame",
             junction.quantification = "data.frame",
             inclusion.levels = "data.frame"),
         prototype = prototype(
             species = "Homo sapiens",
             common.name = "Human",
             ## TODO(NunoA): check if it's possible to circumvent this band-aid
             ## solution that avoids error when date is not filled
             date = structure(list(), class = "Date")))

# Checks validity of Classification's attributes
setValidity("Classification",
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
setAs("Classification", "character", function(from)
    if (from@common.name != "")
        sprintf("%s (%s)", from@species, from@common.name)
    else
        sprintf("%s", from@species))


#' Gets the inclusion levels of a classification
#'
#' @param object Classification
#'
#' @return Inclusion levels of a classification
#' @export
#'
#' @examples
#' org <- new("Classification",
#'            inclusion.levels = data.frame(1,2,3))
#' inclusion.levels(org)
setGeneric("inclusion.levels", function(object)
    standardGeneric("inclusion.levels"))

#' @rdname inclusion.levels
setMethod("inclusion.levels", "Classification", function(object)
    object@inclusion.levels)

#' Sets the inclusion levels of a classification
#'
#' @param object Classification
#' @param value Data.frame: inclusion levels
#' @export
#'
#' @examples
#' org <- new("Classification")
#' inclusion.levels(org) <- data.frame(1,2,3)
setGeneric("inclusion.levels<-", function(object, value)
    standardGeneric("inclusion.levels<-"))

#' @rdname inclusion.levels-set
setReplaceMethod("inclusion.levels", "Classification",
                 function(object, value)
                     initialize(object, inclusion.levels = value))

#' Sets the junction quantification of a classification
#'
#' @param object Classification
#' @param value Data.frame: junction quantification
#' @export
#'
#' @examples
#' org <- new("Classification")
#' junction.quantification(org) <- data.frame(1,2,3)
setGeneric("junction.quantification<-", function(object, value)
    standardGeneric("junction.quantification<-"))

#' @rdname junction.quantification-set
setReplaceMethod("junction.quantification", "Classification",
                 function(object, value)
                     initialize(object, junction.quantification = value))

#' Sets the clinical data of a classification
#'
#' @param object Classification
#' @param value Data.frame: clinical data
#' @export
#'
#' @examples
#' org <- new("Classification")
#' clinical(org) <- data.frame(1,2,3)
setGeneric("clinical<-", function(object, value)
    standardGeneric("clinical<-"))

#' @rdname clinical-set
setReplaceMethod("clinical", "Classification",
                 function(object, value)
                     initialize(object, clinical = value))