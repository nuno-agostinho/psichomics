context("Test data grouping functions")

df <- data.frame(gender=c("male", "female"),
                 stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))

test_that("Create groups by column", {
    group <- createGroupByAttribute(col="stage", dataset=df)
    expect_is(group, "list")
    expect_equal(names(group), paste("stage", 1:4))
    expect_equivalent(unlist(group), c(1, 3, 5, 7, 8, 2, 6, 4))
})

context("Assign one group per patient")

test_that("Each group will be placed in the respective index", {
    # Sequential order
    groups <- list(1:3, 4:7, 8:10)
    names(groups) <- paste("Stage", 1:3)
    
    expect_equal(groupPerPatient(groups, patients = 10),
                 c(rep("Stage 1", 3), rep("Stage 2", 4), rep("Stage 3", 3)))
    
    # Random order
    groups <- list(c(2, 4), c(1, 5), c(3, 6, 7))
    names(groups) <- c("Alive", "Dead", "Zombie")
    expect_equal(
        groupPerPatient(groups, patients = 7),
        c("Dead", "Alive", "Zombie", "Alive", "Dead", "Zombie", "Zombie"))
})

test_that("Each index can belong to multiple groups", {
    # Example 1
    groups <- list(1:3, 4:7, 8:10, c(1, 8))
    names(groups) <- paste("Stage", 1:4)
    expect_equal(groupPerPatient(groups, patients = 10),
                 c("Stage 1, Stage 4", rep("Stage 1", 2), rep("Stage 2", 4), 
                   "Stage 3, Stage 4", rep("Stage 3", 2)))
    
    # Example 2
    groups <- list(c(2, 4), c(1, 5), c(3, 6, 7), c(3, 6, 4))
    names(groups) <- c("Alive", "Dead", "Zombie", "Possibly")
    expect_equal(
        groupPerPatient(groups, patients = 7),
        c("Dead", "Alive", "Zombie, Possibly", "Alive, Possibly", "Dead", 
          "Zombie, Possibly", "Zombie"))
})

test_that("Non-matching patients are returned as NAs or custom group", {
    # Return non-matching patients as NAs
    groups <- list(c(2, 4), c(1, 6), c(9, 10))
    names(groups) <- c("Alive", "Dead", "Zombie")
    expect_equal(
        groupPerPatient(groups, patients = 10),
        c("Dead", "Alive", NA, "Alive", NA, "Dead", NA, NA, "Zombie", "Zombie"))
    
    # Return non-matching patients as part of a custom group
    groups <- list(c(2, 4), c(1, 6), c(9, 10))
    names(groups) <- c("Alive", "Dead", "Zombie")
    expect_equal(
        groupPerPatient(groups, patients = 10, includeOuterGroup = TRUE, 
                        outerGroupName = "Others"),
        c("Dead", "Alive", "Others", "Alive", "Others", "Dead", "Others", 
          "Others", "Zombie", "Zombie"))
})

test_that("No groups returns a custom string", {
    groups <- list()
    expect_equal(groupPerPatient(groups, patients = 10),
                 rep("Single group", 10))
})

context("Assign one group per sample")

test_that("Each group will be placed in the respective index", {
    # Sequential order
    groups <- list(letters[1:3], letters[4:7], letters[8:10])
    names(groups) <- paste("Stage", 1:3)
    samples <- letters[1:10]
    
    expect_equal(groupPerSample(groups, samples),
                 c(rep("Stage 1", 3), rep("Stage 2", 4), rep("Stage 3", 3)))
    
    # Random order
    groups <- list(letters[c(2, 4)], letters[c(1, 5)], letters[c(3, 6, 7)])
    names(groups) <- c("Alive", "Dead", "Zombie")
    samples <- letters[1:7]
    expect_equal(
        groupPerSample(groups, samples),
        c("Dead", "Alive", "Zombie", "Alive", "Dead", "Zombie", "Zombie"))
})

test_that("Each index can belong to multiple groups", {
    # Example 1
    groups <- list(letters[1:3], letters[4:7], letters[8:10], letters[c(1, 8)])
    names(groups) <- paste("Stage", 1:4)
    samples <- letters[1:10]
    expect_equal(groupPerSample(groups, samples),
                 c("Stage 1, Stage 4", rep("Stage 1", 2), rep("Stage 2", 4), 
                   "Stage 3, Stage 4", rep("Stage 3", 2)))
    
    # Example 2
    groups <- list(letters[c(2, 4)], letters[c(1, 5)], letters[c(3, 6, 7)], 
                   letters[c(3, 6, 4)])
    names(groups) <- c("Alive", "Dead", "Zombie", "Possibly")
    samples <- letters[1:7]
    expect_equal(
        groupPerSample(groups, samples),
        c("Dead", "Alive", "Zombie, Possibly", "Alive, Possibly", "Dead", 
          "Zombie, Possibly", "Zombie"))
})

test_that("Non-matching samples are returned as NAs or custom group", {
    # Return non-matching patients as NAs
    groups <- list(letters[c(2, 4)], letters[c(1, 6)], letters[c(9, 10)])
    names(groups) <- c("Alive", "Dead", "Zombie")
    samples <- letters[1:10]
    expect_equal(
        groupPerSample(groups, samples),
        c("Dead", "Alive", NA, "Alive", NA, "Dead", NA, NA, "Zombie", "Zombie"))
    
    # Return non-matching patients as part of a custom group
    groups <- list(letters[c(2, 4)], letters[c(1, 6)], letters[c(9, 10)])
    names(groups) <- c("Alive", "Dead", "Zombie")
    samples <- letters[1:10]
    expect_equal(
        groupPerSample(groups, samples, includeOuterGroup = TRUE, 
                       outerGroupName = "Others"),
        c("Dead", "Alive", "Others", "Alive", "Others", "Dead", "Others", 
          "Others", "Zombie", "Zombie"))
})

test_that("No groups returns a custom string", {
    groups <- list()
    expect_equal(groupPerSample(groups, letters[1:10]), rep("Single group", 10))
})