context("Test data grouping functions")

df <- data.frame(gender=c("male", "female"),
                 stage=paste("stage", c(1, 3, 1, 4, 2, 3, 2, 2)))

test_that("Create groups by column", {
    group <- createGroupByAttribute(col="stage", dataset=df)
    expect_is(group, "list")
    expect_equal(names(group), paste("stage", 1:4))
    expect_equivalent(unlist(group), as.character(c(1, 3, 5, 7, 8, 2, 6, 4)))
})

context("Assign one group per element")

test_that("Each group will be placed in the respective index", {
    groups <- list(1:3, 4:7, 8:10)
    names(groups) <- paste("Stage", 1:3)
    
    expected <- c(rep("Stage 1", 3), rep("Stage 2", 4), rep("Stage 3", 3))
    names(expected) <- unlist(groups)
    
    expect_equal(groupPerElem(groups), expected)
})

test_that("Each index can belong to multiple groups", {
    # Example 1
    groups <- list(1:3, 4:7, 8:10, c(1, 8))
    names(groups) <- paste("Stage", 1:4)
    
    expected <- c("Stage 1, Stage 4", rep("Stage 1", 2), rep("Stage 2", 4), 
                  "Stage 3, Stage 4", rep("Stage 3", 2))
    names(expected) <- unique(unlist(groups))
    
    expect_equal(groupPerElem(groups), expected)
    
    # Example 2
    groups <- list(c(2, 4), c(1, 5), c(3, 6, 7), c(3, 6, 4))
    names(groups) <- c("Alive", "Dead", "Zombie", "Possibly")
    
    expected <- c("Alive", "Alive, Possibly", "Dead", "Dead",
                  "Zombie, Possibly", "Zombie, Possibly", "Zombie")
    names(expected) <- unique(unlist(groups))
    
    expect_equal(groupPerElem(groups), expected)
})

test_that("Non-matching subjects are returned as NAs or custom group", {
    # Return non-matching subjects as NAs
    groups <- list(c(2, 4), c(1, 6), c(9, 10))
    names(groups) <- c("Alive", "Dead", "Zombie")
    
    elem <- 1:10
    expected <- c("Dead", "Alive", NA, "Alive", NA, "Dead", NA, NA, "Zombie", 
                  "Zombie")
    names(expected) <- elem
        
    expect_equal(groupPerElem(groups, elem), expected)
    
    # Return non-matching subjects as part of a custom group
    groups <- list(c(2, 4), c(1, 6), c(9, 10))
    names(groups) <- c("Alive", "Dead", "Zombie")
    
    elem <- 1:10
    expected <- c("Dead", "Alive", "Others", "Alive", "Others", "Dead", 
                  "Others", "Others", "Zombie", "Zombie")
    names(expected) <- elem
    expect_equal(groupPerElem(groups, elem, outerGroupName="Others"), expected)
})

test_that("No groups returns a custom string", {
    groups <- list()
    elem <- 1:10
    expect_equal(groupPerElem(groups, elem), rep("Single group", length(elem)))
})

context("Test set operations") ################################################

# Prepare groups containing both subjects and samples
dummySampleId <- function(i)
    if (length(i) > 0) paste0("sample-", i, "-test")
dummySubjectId <- function(i)
    if (length(i) > 0) paste0("subject-", i)

matches <- list("1"=c(1:2), "2"=c(3:5), "3"=c(6:10), "4"=c(11:15), "5"=c(),
                "6"=c(16), "7"=c(17, 18, 20))
matches <- sapply(matches, dummySampleId)
names(matches) <- dummySubjectId(names(matches))

inverted <- c("1"=1, "2"=1, "3"=2, "4"=2, "5"=2, "6"=3, "7"=3, "8"=3, "9"=3,
              "10"=3, "11"=4, "12"=4, "13"=4, "14"=4, "15"=4, "16"=6, 
              "17"=7, "18"=7, "19"=NULL, "20"=7)
ns <- dummySampleId(names(inverted))
inverted <- dummySubjectId(inverted)
names(inverted) <- ns

returnSamples   <- function(subjects) {
    res <- unique(unname(unlist(matches[subjects], use.names=FALSE)))
    res[!is.na(res)]
}
returnSubjects  <- function(samples) {
    res <- unique(unname(unlist(inverted[samples], use.names=FALSE)))
    res[!is.na(res)]
}

prepareTestGroups <- function(matches, inverted) {
    male   <- dummySubjectId(1:3)
    female <- dummySubjectId(4:7)
    maleSamples   <- returnSamples(male)
    femaleSamples <- returnSamples(female)
    
    normal <- dummySampleId(c(1, 3, 6, 10, 11, 15:18))
    tumour <- dummySampleId(c(2, 4, 5, 7:9, 12:14, 20))
    normalMatch <- returnSubjects(normal)
    tumourMatch <- returnSubjects(tumour)
    
    df <- rbind(
        cbind("male",   "Attr", "gender",      list(male),   list(maleSamples)),
        cbind("female", "Attr", "gender",      list(female), list(femaleSamples)),
        cbind("normal", "Attr", "sample_type", list(normalMatch), list(normal)),
        cbind("tumour", "Attr", "sample_type", list(tumourMatch), list(tumour)))
    
    colnames(df) <- c("Names", "Subset", "Input", "Patients", "Samples")
    rownames(df) <- df[ , "Names"]
    return(df)
}
df <- prepareTestGroups(matches, inverted)

test_that("Set union", {
    # Test with 1 group
    selected <- 2
    df2 <- setOperation("union", df, selected, matches=inverted)
    
    samples  <- Reduce(union, df[selected, "Samples"])
    subjects <- Reduce(union, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 2 groups
    selected <- 2:3
    df2 <- setOperation("union", df, selected, matches=inverted)
    
    samples  <- Reduce(union, df[selected, "Samples"])
    subjects <- Reduce(union, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 4 groups
    selected <- 1:4
    df2 <- setOperation("union", df, selected, matches=inverted)
    
    samples  <- Reduce(union, df[selected, "Samples"])
    subjects <- Reduce(union, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
})

test_that("Set intersect", {
    # Test with 1 group
    selected <- 2
    df2 <- setOperation("intersect", df, selected, matches=inverted)
    
    samples  <- Reduce(intersect, df[selected, "Samples"])
    subjects <- Reduce(intersect, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("intersect", df, selected, matches=inverted)
    
    samples  <- Reduce(intersect, df[selected, "Samples"])
    subjects <- Reduce(intersect, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 4 groups
    selected <- 1:4
    df2 <- setOperation("intersect", df, selected, matches=inverted)
    
    samples  <- Reduce(intersect, df[selected, "Samples"])
    subjects <- Reduce(intersect, df[selected, "Patients"])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
})

test_that("Set complement", {
    allSamples <- dummySampleId(1:20)
    allSubjects <- dummySubjectId(1:7)
    
    # Test with 1 group
    selected <- 3
    df2 <- setOperation("complement", df, selected, first=allSubjects,
                        second=allSamples, symbol="U \u005c ", matches=inverted)
    
    samples  <- setdiff(allSamples, df[[selected, "Samples"]])
    subjects <- setdiff(allSubjects, df[[selected, "Patients"]])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("complement", df, selected, first=allSubjects,
                        second=allSamples, symbol="U \u005c ", matches=inverted)
    
    samples  <- setdiff(allSamples, Reduce(union, df[selected, "Samples"]))
    subjects <- setdiff(allSubjects, Reduce(union, df[selected, "Patients"]))
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 4 groups
    selected <- c(1:4)
    df2 <- setOperation("complement", df, selected, first=allSubjects,
                        second=allSamples, symbol="U \u005c ", matches=inverted)
    
    samples  <- setdiff(allSamples, Reduce(union, df[selected, "Samples"]))
    subjects <- setdiff(allSubjects, Reduce(union, df[selected, "Patients"]))
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
})

test_that("Set subtract", {
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("subtract", df, selected, matches=inverted)
    
    samples  <- setdiff(df[[selected[1], "Samples"]],
                        df[[selected[2], "Samples"]])
    subjects <- setdiff(df[[selected[1], "Patients"]],
                        df[[selected[2], "Patients"]])
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Error if only 1 group is provided
    selected <- 2
    expect_error(setOperation("subtract", df, selected, matches=inverted), 
                 "set subtract requires 2 groups")
    
    # Error if more than 2 groups are provided
    selected <- 1:4
    expect_error(setOperation("subtract", df, selected, matches=inverted), 
                 "set subtract requires 2 groups")
})

test_that("Set symmetric difference", {
    # Test with 1 group
    selected <- 2
    df2 <- setOperation("symDiff", df, selected, matches=inverted)
    
    samples  <- setdiff(Reduce(union, df[selected, "Samples"]),
                        Reduce(intersect, df[selected, "Samples"]))
    subjects <- setdiff(Reduce(union, df[selected, "Patients"]),
                        Reduce(intersect, df[selected, "Patients"]))
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 2 groups
    selected <- c(2, 4)
    df2 <- setOperation("symDiff", df, selected, matches=inverted)
    
    samples  <- setdiff(Reduce(union, df[selected, "Samples"]),
                        Reduce(intersect, df[selected, "Samples"]))
    subjects <- setdiff(Reduce(union, df[selected, "Patients"]),
                        Reduce(intersect, df[selected, "Patients"]))
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
    
    # Test with 4 groups
    selected <- 1:4
    df2 <- setOperation("symDiff", df, selected, matches=inverted)
    
    samples  <- setdiff(Reduce(union, df[selected, "Samples"]),
                        Reduce(intersect, df[selected, "Samples"]))
    subjects <- setdiff(Reduce(union, df[selected, "Patients"]),
                        Reduce(intersect, df[selected, "Patients"]))
    subjects <- unique(c(subjects, returnSubjects(samples)))
    
    expect_equal(df2[[1, "Samples"]], samples)
    expect_equal(df2[[1, "Patients"]], subjects)
})

test_that("Rename groups", {
    # Rename a group
    selected <- 3
    name <- "non-tumour"
    df2 <- setOperation("rename", df, selected, matches=inverted,
                        groupName=name)
    expect_equal(df2[[selected, "Names"]], name)
    
    # Rename multiple groups
    selected <- c(1, 3)
    name <- c("non-female", "non-tumour")
    df2 <- setOperation("rename", df, selected, matches=inverted, 
                        groupName=name)
    expect_equal(as.character(df2[selected, "Names"]), name)
    
    # Properly name a new group after a set operation
    selected <- c(2, 4)
    name <- "(female tumour)"
    df2 <- setOperation("union", df, selected, matches=inverted, groupName=name)
    expect_equal(df2[[1, "Names"]], name)
    
    ### Avoid renaming groups with a previously used name ###
    
    # Rename a group with a previously used name
    selected <- 3
    name <- "tumour"
    df2 <- setOperation("rename", df, selected, matches=inverted, 
                        groupName=name)
    expect_equal(as.character(df2[selected, "Names"]), paste(name, "(1)"))
    
    # Rename multiple groups with the same name
    selected <- c(1, 3)
    name <- "group of interest"
    df2 <- setOperation("rename", df, selected, matches=inverted, 
                        groupName=name)
    expect_equal(as.character(df2[selected, "Names"]),
                 c(name, paste(name, "(1)")))
    
    # Rename multiple groups with a previously used name
    selected <- c(1, 3)
    name <- "tumour"
    df2 <- setOperation("rename", df, selected, matches=inverted, 
                        groupName=name)
    expect_equal(as.character(df2[selected, "Names"]),
                 paste(name, c("(1)", "(2)")))
})

test_that("Remove groups", {
    # Remove a group
    selected <- 2
    name <- as.character(df[selected, "Names"])
    df2 <- setOperation("remove", df, selected, matches=inverted)
    expect_false(name %in% as.character(df2[, "Names"]))
    
    # Remove multiple groups
    selected <- c(2, 4)
    name <- as.character(df[selected, "Names"])
    df2 <- setOperation("remove", df, selected, matches=inverted)
    expect_false( any(name %in% as.character(df2[, "Names"])) )
})
