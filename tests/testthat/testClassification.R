context("Classification class")

test_that("Classification object is created", {
    classif <- new("Classification")
    expect_equal(classif@species, "Homo sapiens")
    expect_equal(classif@common.name, "Human")
    
    classif <- new("Classification",
               species = "Mus musculus",
               common.name = "Mouse")
    expect_equal(classif@species, "Mus musculus")
    expect_equal(classif@common.name, "Mouse")
    
    classif <- new("Classification",
               inclusion.levels = data.frame(1, 2, 3),
               clinical = data.frame("a", "b", "c"))
    expect_equal(classif@inclusion.levels, data.frame(1, 2, 3))
    expect_equal(classif@clinical, data.frame("a", "b", "c"))
})

test_that("Classification object's attributes are validated", {
    # 'species' should be a single string
    expect_error( new("Classification", species = c("Homo sapiens", "Mus musculus")) )
    
    # 'species' should not be empty
    expect_error( new("Classification", species = "" ))
    
    # 'common.name' should be a single string
    expect_error( new("Classification", common.name = c("Humna", "Mouse")) )
})

test_that("Classification has a string representation", {
    classif <- new("Classification")
    expect_equal(as(classif, "character"), "Homo sapiens (Human)")
    
    classif <- new("Classification", common.name = "")
    expect_equal(as(classif, "character"), "Homo sapiens")
})