setClass("Organism",
         representation = representation(
             species = "character",
             common.name = "character",
             clinical.information = "data.frame",
             inclusion.levels = "data.frame"),
         prototype = prototype(
             species = "Homo sapiens",
             common.name = "Human"))

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

setAs("Organism", "character", function(from)
    if (from@common.name != "")
        sprintf("%s (%s)", from@species, from@common.name)
    else
        sprintf("%s", from@species))

setGeneric("inclusion.levels", function(object)
    standardGeneric("inclusion.levels"))
setMethod("inclusion.levels", "Organism", function(object)
    object@inclusion.levels)

setGeneric("inclusion.levels<-", function(object, value)
    standardGeneric("inclusion.levels<-"))
setReplaceMethod("inclusion.levels", "Organism", function(object, value)
    initialize(object, inclusion.levels = value))

# -----------------------------------------------------------------------
# org <- new("Organism",
#            species = "Mus musculus",
#            common.name = "Mouse")
# as(org, "character")
# inclusion.levels(org)
# inclusion.levels(org) <- data.frame(1,2,3)
# inclusion.levels(org)
# 
# org <- new("Organism",
#            species = "Mus musculus",
#            common.name = "Mouse",
#            inclusion.levels = data.frame(1,2,3))