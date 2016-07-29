context("Assign one group per patient")

test_that("Each group will be placed in the respective index", {
    # Sequential order
    names <- paste("Stage", 1:3)
    groups <- matrix(ncol=2, byrow=TRUE, dimnames=list(
        names, 
        c("Names", "Rows")),
        c(names[1], list(1:3), 
          names[2], list(4:7),
          names[3], list(8:10)))
    expect_equal(groupPerPatient(groups, patients = 10),
                 c(rep("Stage 1", 3), rep("Stage 2", 4), rep("Stage 3", 3)))
    
    # Random order
    names <- c("Alive", "Dead", "Zombie")
    groups <- matrix(ncol=2, byrow=TRUE, dimnames=list(
        names, 
        c("Names", "Rows")),
        c(names[1], list(c(2, 4)), 
          names[2], list(c(1, 5)),
          names[3], list(c(3, 6, 7))))
    expect_equal(
        groupPerPatient(groups, patients = 7),
        c("Dead", "Alive", "Zombie", "Alive", "Dead", "Zombie", "Zombie"))
})

test_that("Each index can belong to multiple groups", {
    # Example 1
    names <- paste("Stage", 1:4)
    groups <- matrix(ncol=1, byrow=TRUE, dimnames=list(names, "Rows"),
                     c(list(1:3), list(4:7), list(8:10), list(c(1, 8))))
    expect_equal(groupPerPatient(groups, patients = 10),
                 c("Stage 1, Stage 4", rep("Stage 1", 2), rep("Stage 2", 4), 
                   "Stage 3, Stage 4", rep("Stage 3", 2)))
    
    # Example 2
    names <- c("Alive", "Dead", "Zombie", "Possibly")
    groups <- matrix(ncol=1, byrow=TRUE, dimnames=list(names, "Rows"),
                     c(list(c(2, 4)), list(c(1, 5)), list(c(3, 6, 7)),
                       list(c(3, 6, 4))))
    expect_equal(
        groupPerPatient(groups, patients = 7),
        c("Dead", "Alive", "Zombie, Possibly", "Alive, Possibly", "Dead", 
          "Zombie, Possibly", "Zombie"))
})

test_that("Non-matching patients are returned as NAs or custom group", {
    # Return non-matching patients as NAs
    names <- c("Alive", "Dead", "Zombie")
    groups <- matrix(ncol=2, byrow=TRUE, dimnames=list(
        names, 
        c("Names", "Rows")),
        c(names[1], list(c(2, 4)), 
          names[2], list(c(1, 6)),
          names[3], list(c(9, 10))))
    expect_equal(
        groupPerPatient(groups, patients = 10),
        c("Dead", "Alive", NA, "Alive", NA, "Dead", NA, NA, "Zombie", "Zombie"))
    
    # Return non-matching patients as part of a custom group
    names <- c("Alive", "Dead", "Zombie")
    groups <- matrix(ncol=2, byrow=TRUE, dimnames=list(
        names, 
        c("Names", "Rows")),
        c(names[1], list(c(2, 4)), 
          names[2], list(c(1, 6)),
          names[3], list(c(9, 10))))
    expect_equal(
        groupPerPatient(groups, patients = 10, includeOuterGroup = TRUE, 
                        outerGroupName = "Others"),
        c("Dead", "Alive", "Others", "Alive", "Others", "Dead", "Others", 
          "Others", "Zombie", "Zombie"))
})

test_that("No groups returns a custom string", {
    groups <- matrix(ncol=2, byrow=TRUE, 
                     dimnames=list("Single", c("Names", "Rows")),
                     c("Single", list(c(2, 4))))
    groups <- groups[-1, ]
    expect_equal(
        groupPerPatient(groups, patients = 10, allDataName = "No groups"),
        rep("No groups", 10))
})