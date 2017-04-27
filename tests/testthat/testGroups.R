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

context("Test set operations") ################################################

matches  <- list("1"=c(1:2), "2"=c(3:5), "3"=c(6:10), "4"=c(11:15), "5"=c(),
                 "6"=c(16), "7"=c(17, 18, 20))
inverted <- list("1"=1, "2"=1, "3"=2, "4"=2, "5"=2, "6"=3, "7"=3, "8"=3, "9"=3,
                 "10"=3, "11"=4, "12"=4, "13"=4, "14"=4, "15"=4,
                 "16"=6, "17"=7, "18"=7, "19"=NULL, "20"=7)

returnSamples   <- function(patients) unlist(matches[patients], use.names=FALSE)
returnPatients  <- function(samples)  unlist(inverted[samples], use.names=FALSE)

male   <- 1:3
female <- 4:7
maleSamples   <- returnSamples(male)
femaleSamples <- returnSamples(female)

normal <- 1:10
tumour <- 11:20
normalMatch <- returnPatients(normal)
tumourMatch <- returnPatients(tumour)

df <- rbind(
    cbind("male",   "Attr", "gender",      list(male),   list(maleSamples)),
    cbind("female", "Attr", "gender",      list(female), list(femaleSamples)),
    cbind("normal", "Attr", "sample_type", list(normalMatch), list(normal)),
    cbind("tumour", "Attr", "sample_type", list(tumourMatch), list(tumour)))

colnames(df) <- c("Names", "Subset", "Input", "Patients", "Samples")
rownames(df) <- df[ , "Names"]

test_that("Set union", {
    # Test with 2 groups
    selected <- 2:3
    df2 <- setOperation("union", df, selected)
    
    samples  <- Reduce(union, df[selected, "Samples"])
    patients <- Reduce(union, df[selected, "Patients"])
    patients <- unique(c(patients, returnPatients(samples)))
    
    expect_equal(sort(df2[[1, "Samples"]]), sort(samples))
    expect_equal(sort(df2[[1, "Patients"]]), sort(patients))
})

test_that("Set intersect", {
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("intersect", df, selected)
    
    samples  <- Reduce(intersect, df[selected, "Samples"])
    patients <- Reduce(intersect, df[selected, "Samples"])
    patients <- unique(c(patients, returnPatients(samples)))
    
    expect_equal(sort(df2[[1, "Samples"]]), sort(samples))
    expect_equal(sort(df2[[1, "Patients"]]), sort(patients))
})

test_that("Set complement", {
    # Test with 2 groups
    allSamples <- 1:20
    allPatients <- 1:7
    selected <- c(2, 4)
    
    df2 <- setOperation("complement", df, selected, patients=allPatients,
                        samples=allSamples, symbol="U \u005c ")
    
    samples  <- setdiff(allSamples, Reduce(union, df[selected, "Samples"]))
    patients <- setdiff(allPatients, Reduce(union, df[selected, "Patients"]))
    patients <- unique(c(patients, returnPatients(samples)))
    
    expect_equal(sort(df2[[1, "Samples"]]), sort(samples))
    expect_equal(sort(df2[[1, "Patients"]]), sort(patients))
})

test_that("Set subtract", {
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("subtract", df, selected)
    
    samples  <- setdiff(df[[selected[1], "Samples"]],
                        df[[selected[2], "Samples"]])
    patients <- setdiff(df[[selected[1], "Patients"]],
                        df[[selected[2], "Patients"]])
    patients <- unique(c(patients, returnPatients(samples)))
    
    expect_equal(sort(df2[[1, "Samples"]]), sort(samples))
    expect_equal(sort(df2[[1, "Patients"]]), sort(patients))
})

test_that("Set symmetric difference", {
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("symDiff", df, selected)
    
    samples  <- setdiff(Reduce(union, df[selected, "Samples"]),
                        Reduce(intersect, df[selected, "Samples"]))
    patients <- setdiff(Reduce(union, df[selected, "Patients"]),
                        Reduce(intersect, df[selected, "Patients"]))
    patients <- unique(c(patients, returnPatients(samples)))
    
    expect_equal(sort(df2[[1, "Samples"]]), sort(samples))
    expect_equal(sort(df2[[1, "Patients"]]), sort(patients))
})