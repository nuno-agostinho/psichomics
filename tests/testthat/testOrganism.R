context("Organism class")

test_that("Organism object is created", {
    org <- new("Organism")
    expect_equal(org@species, "Homo sapiens")
    expect_equal(org@common.name, "Human")
    
    org <- new("Organism",
               species = "Mus musculus",
               common.name = "Mouse")
    expect_equal(org@species, "Mus musculus")
    expect_equal(org@common.name, "Mouse")
    
    org <- new("Organism",
               inclusion.levels = data.frame(1, 2, 3),
               clinical.information = data.frame("a", "b", "c"))
    expect_equal(org@inclusion.levels, data.frame(1, 2, 3))
    expect_equal(org@clinical.information, data.frame("a", "b", "c"))
})

test_that("Organism object's attributes are validated", {
    # 'species' should be a single string
    expect_error( new("Organism", species = c("Homo sapiens", "Mus musculus")) )
    
    # 'species' should not be empty
    expect_error( new("Organism", species = "" ))
    
    # 'common.name' should be a single string
    expect_error( new("Organism", common.name = c("Humna", "Mouse")) )
})

test_that("Organism has a string representation", {
    org <- new("Organism")
    expect_equal(as(org, "character"), "Homo sapiens (Human)")
    
    org <- new("Organism", common.name = "")
    expect_equal(as(org, "character"), "Homo sapiens")
})