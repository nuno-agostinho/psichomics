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