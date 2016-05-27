context("Auxiliary functions to parse events")

test_that("createJunctionsTemplate creates a template of junctions with NAs", {
    nrow <- 8
    temp <- createJunctionsTemplate(nrow)
    expect_equal(nrow(temp), nrow)
    expect_true(all(is.na(temp)))
    expect_equal(names(temp), c("C1.start", "C1.end", "A1.start", "A1.end",
                                "A2.start", "A2.end", "C2.start", "C2.end"))
    
    nrow <- 2
    temp <- createJunctionsTemplate(nrow, program="MISO",
                                    event.type=c("A5SS", "SE"), 
                                    chromosome=c(2, 4),
                                    strand=c("-", "+"))
    expect_equal(nrow(temp), nrow)
    expect_false(all(is.na(temp)))
    expect_equal(names(temp), c("C1.start", "C1.end", "A1.start", "A1.end",
                                "A2.start", "A2.end", "C2.start", "C2.end",
                                "Program", "Event.type", "Chromosome", "Strand"))
    expect_equal(temp$Program, c("MISO", "MISO"))
    expect_equal(temp$Event.type, c("A5SS", "SE"))
    expect_equal(temp$Chromosome, c(2, 4))
    expect_equal(temp$Strand, c("-", "+"))
})

test_that("Calculate inclusion levels for skipping exon", {
    library(fastmatch)
    
    eventType <- "SE"
    annot <- read.table(text = "1 + 32 35 37 38
                                2 + 32 35 37 38
                                3 + 32 35 37 38")
    names(annot) <- c("Chromosome", "Strand",
                      "C1.end", "A1.start", "A1.end", "C2.start")
    junctionQuant <- read.table(text = "10 10 10 10 10 10
                                        10 10 10 10 10 10
                                        10 10 10 10 10 10
                                        27 20 90 24 14 35
                                        10 18 13 12 10 21
                                        30 24 92 26 13 29
                                        27 20 90 24 14 35
                                        90 98 93 92 90 91
                                        30 24 92 26 13 29")
    names(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
    rownames(junctionQuant) <- c("chr1:32:+,chr1:35:+",
                                 "chr1:32:+,chr1:38:+",
                                 "chr1:37:+,chr1:38:+",
                                 "chr2:32:+,chr2:35:+",
                                 "chr2:32:+,chr2:38:+",
                                 "chr2:37:+,chr2:38:+",
                                 "chr3:32:+,chr3:35:+",
                                 "chr3:32:+,chr3:38:+",
                                 "chr3:37:+,chr3:38:+")
    psi <- calculateInclusionLevels(eventType, junctionQuant, annot)
    
    expect_true(all(psi[1, ] == 0.5)) # Same reads for all junctions
    expect_true(all(psi[2, ] > 0.5)) # More reads for inclusion isoform
    expect_true(all(psi[3, ] < 0.5)) # More reads for exclusive isoform
})

test_that("Calculate inclusion levels for skipping exon with minimum reads", {
    library(fastmatch)
    
    eventType <- "SE"
    minReads <- 15
    annot <- read.table(text = "1 + 32 35 37 38
                                2 + 32 35 37 38
                                3 + 32 35 37 38")
    names(annot) <- c("Chromosome", "Strand",
                      "C1.end", "A1.start", "A1.end", "C2.start")
    junctionQuant <- read.table(text = "10 10 10 10 10 10
                                        10 10 10 10 10 10
                                        10 10 10 10 10 10
                                        27 20 90 24 14 35
                                        10 18 13 12 10 21
                                        30 24 92 26 13 29
                                        27 20 90 24 14 35
                                        90 98 93 92 90 91
                                        30 24 92 26 13 29")
    names(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
    rownames(junctionQuant) <- c("chr1:32:+,chr1:35:+",
                                 "chr1:32:+,chr1:38:+",
                                 "chr1:37:+,chr1:38:+",
                                 "chr2:32:+,chr2:35:+",
                                 "chr2:32:+,chr2:38:+",
                                 "chr2:37:+,chr2:38:+",
                                 "chr3:32:+,chr3:35:+",
                                 "chr3:32:+,chr3:38:+",
                                 "chr3:37:+,chr3:38:+")
    psi <- calculateInclusionLevels(eventType, junctionQuant, annot, minReads)
    
    expect_equal(nrow(psi), 2) # Discard first event based on few reads
    expect_true(is.na(psi[1, 5]))
    expect_true(all(psi[1, -5] > 0.5)) # More reads for inclusion isoform
    expect_true(all(psi[2, ] < 0.5)) # More reads for exclusive isoform
})

test_that("Calculate inclusion levels for mutually exclusive exons", {
    library(fastmatch)
    
    eventType <- "MX"
    annot <- read.table(text = "1 + 32 35 37 38 40 42
                                2 + 32 35 37 38 40 42
                                3 + 32 35 37 38 40 42")
    names(annot) <- c("Chromosome", "Strand",
                      "C1.end", "A1.start", "A1.end",
                      "A2.start", "A2.end", "C2.start")
    junctionQuant <- read.table(text = "10 10 10 10 10 10
                                10 10 10 10 10 10
                                10 10 10 10 10 10
                                10 10 10 10 10 10
                                27 20 90 24 14 35
                                30 24 92 26 13 29
                                10 18 13 12 10 21
                                10 14 12 16 13 19
                                10 18 13 12 10 21
                                10 14 12 16 13 19
                                27 20 90 24 14 35
                                30 24 92 26 13 29")
    names(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
    rownames(junctionQuant) <- c("chr1:32:+,chr1:35:+",
                                 "chr1:37:+,chr1:42:+",
                                 "chr1:32:+,chr1:38:+",
                                 "chr1:40:+,chr1:42:+",
                                 "chr2:32:+,chr2:35:+",
                                 "chr2:37:+,chr2:42:+",
                                 "chr2:32:+,chr2:38:+",
                                 "chr2:40:+,chr2:42:+",
                                 "chr3:32:+,chr3:35:+",
                                 "chr3:37:+,chr3:42:+",
                                 "chr3:32:+,chr3:38:+",
                                 "chr3:40:+,chr3:42:+")
    psi <- calculateInclusionLevels(eventType, junctionQuant, annot)
    
    expect_true(all(psi[1, ] == 0.5)) # Same reads for all junctions
    expect_true(all(psi[2, ] > 0.5)) # More reads for inclusive isoform
    expect_true(all(psi[3, ] < 0.5)) # More reads for exclusive isoform
})

test_that("Calculate inclusion levels for alternative 5' splice site", {
    library(fastmatch)
    
    eventType <- "A5"
    annot <- read.table(text = "1 + 32 35 37
                                2 + 32 35 37
                                3 + 32 35 37")
    names(annot) <- c("Chromosome", "Strand", "C1.end", "A1.end", "C2.start")
    junctionQuant <- read.table(text = "10 10 10 10 10 10
                                10 10 10 10 10 10
                                10 14 12 16 13 19
                                27 20 90 24 14 35
                                30 24 92 26 13 29
                                10 18 13 12 10 21")
    names(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
    rownames(junctionQuant) <- c("chr1:32:+,chr1:37:+",
                                 "chr1:35:+,chr1:37:+",
                                 "chr2:32:+,chr2:37:+",
                                 "chr2:35:+,chr2:37:+",
                                 "chr3:32:+,chr3:37:+",
                                 "chr3:35:+,chr3:37:+")
    psi <- calculateInclusionLevels(eventType, junctionQuant, annot)
    
    expect_true(all(psi[1, ] == 0.5)) # Same reads for all junctions
    expect_true(all(psi[2, ] > 0.5)) # More reads for inclusive isoform
    expect_true(all(psi[3, ] < 0.5)) # More reads for exclusive isoform
})

test_that("Calculate inclusion levels for alternative 3' splice site", {
    library(fastmatch)
    
    eventType <- "A3"
    annot <- read.table(text = "1 + 32 35 37
                                2 + 32 35 37
                                3 + 32 35 37")
    names(annot) <- c("Chromosome", "Strand", "C1.end", "A1.start", "C2.start")
    junctionQuant <- read.table(text = "10 10 10 10 10 10
                                        10 10 10 10 10 10
                                        27 20 90 24 14 35
                                        10 14 12 16 13 19
                                        10 18 13 12 10 21
                                        30 24 92 26 13 29")
    names(junctionQuant) <- c(paste("Normal", 1:3), paste("Cancer", 1:3))
    rownames(junctionQuant) <- c("chr1:32:+,chr1:35:+",
                                 "chr1:32:+,chr1:37:+",
                                 "chr2:32:+,chr2:35:+",
                                 "chr2:32:+,chr2:37:+",
                                 "chr3:32:+,chr3:35:+",
                                 "chr3:32:+,chr3:37:+")
    psi <- calculateInclusionLevels(eventType, junctionQuant, annot)
    
    expect_true(all(psi[1, ] == 0.5)) # Same reads for all junctions
    expect_true(all(psi[2, ] > 0.5)) # More reads for inclusive isoform
    expect_true(all(psi[3, ] < 0.5)) # More reads for exclusive isoform
})