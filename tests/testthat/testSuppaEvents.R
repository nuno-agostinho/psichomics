context("Parse SUPPA splicing events")

test_that("parseSuppaAnnotation parses annotation from SUPPA", {
    folder <- "extdata/eventsAnnotSample/suppa_output/suppaEvents"
    suppaOutput <- system.file(folder, package="psichomics")
    
    suppa <- parseSuppaAnnotation(suppaOutput)
    expect_is(suppa, "ASevents")
    expect_equal(length(suppa), 14)
    expect_equal(unique(suppa$Program), "SUPPA")
    expect_equal(unique(suppa$Chromosome), "20")
    expect_equal(unique(suppa$Strand), "+")
})

test_that("parseSuppaEvent parses multiple exon skipping event IDs at once", {
    # Load all types of events to test
    events <- c(
        "ENSG00000131002.7;SE:chrY:21751498-21753666:21753845-21755285:+",
        "ENSG00000131002.7;SE:chrY:21759551-21760438:21760525-21761625:+",
        "ENSG00000131002.7;SE:chrY:21729837-21731271:21731345-21749096:+",
        "ENSG00000147761.4;SE:chrY:9544678-9544925:9545180-9546154:-")
    expect_silent(parsed <- parseSuppaEvent(events))
    expect_is(parsed, "data.frame")
    # number of elements in list is the same as number of events
    expect_equal(nrow(parsed), length(events))
    expect_equal(parsed$Program, rep("SUPPA", 4))
    expect_equal(parsed$Event.type, rep("SE", 4))
    expect_equal(parsed$C1.end, c("21751498", "21759551", "21729837",
                                  "9546154"))
})

test_that("parseSuppaEvent parses multiple alt. 3' SS event IDs at once", {
    # Load all types of events to test
    events <- c(
        "ENSG00000260117;A3:HG531_PATCH:115086707-115089263:115086707-115089266:+",
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153768260-153768394:153768260-153768553:+",
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153769180-153769533:153769180-153769638:+",
        "ENSG00000131591;A3:1:1019466-1019733:1019391-1019733:-")
    expect_silent(parsed <- parseSuppaEvent(events))
    expect_is(parsed, "data.frame")
    # number of elements in list is the same as number of events
    expect_equal(nrow(parsed), length(events))
    expect_equal(parsed$Program, rep("SUPPA", 4))
    expect_equal(parsed$Event.type, rep("A3SS", 4))
    expect_equal(parsed$C1.end, c("115086707", "153768260", "153769180",
                                  "1019733"))
})

test_that("parseSuppaSE parses a exon skipping event junctions", {
    # Strand plus
    junctions <- read.table(text = "169768099 169770024 169770112 169771762")
    parsed <- parseSuppaSE(junctions, "+")
    expect_equal(parsed$C1.end,   169768099)
    expect_equal(parsed$A1.start, 169770024)
    expect_equal(parsed$A1.end,   169770112)
    expect_equal(parsed$C2.start, 169771762)
    
    # Strand minus
    junctions <- read.table(text = "1192510 1192588 1192690 1198726")
    parsed <- parseSuppaSE(junctions, "-")
    expect_equal(parsed$C1.end,   1198726)
    expect_equal(parsed$A1.start, 1192690)
    expect_equal(parsed$A1.end,   1192588)
    expect_equal(parsed$C2.start, 1192510)
})

test_that("parseSuppaMXE parses a mutually exclusive exon event junctions", {
    # Strand plus
    junctions <- read.table(
        text = "202060671 202068453 202068489 202073793 202060671 202072798 202072906 202073793")
    parsed <- parseSuppaMXE(junctions, "+")
    expect_equal(parsed$C1.end,   202060671)
    expect_equal(parsed$A1.start, 202068453)
    expect_equal(parsed$A1.end,   202068489)
    expect_equal(parsed$A2.start, 202072798)
    expect_equal(parsed$A2.end,   202072906)
    expect_equal(parsed$C2.start, 202073793)
    
    # Strand minus
    junctions <- read.table(
        text = "220236296 220239415 220239459 220240644 220236296 220240372 220240467 220240644")
    parsed <- parseSuppaMXE(junctions, "-")
    expect_equal(parsed$C1.end,   220240644)
    expect_equal(parsed$A1.start, 220240467)
    expect_equal(parsed$A1.end,   220240372)
    expect_equal(parsed$A2.start, 220239459)
    expect_equal(parsed$A2.end,   220239415)
    expect_equal(parsed$C2.start, 220236296)
})

test_that("parseSuppaRI parses a retained intron event junctions", {
    # Strand plus
    junctions <- read.table(text = "196709749 196709922 196711005 196711181")
    parsed <- parseSuppaRI(junctions, "+")
    expect_equal(parsed$C1.start, 196709749)
    expect_equal(parsed$C1.end,   196709922)
    expect_equal(parsed$C2.start, 196711005)
    expect_equal(parsed$C2.end,   196711181)
    
    # Strand minus
    junctions <- read.table(text = "1326146 1326955 1328170 1328183")
    parsed <- parseSuppaRI(junctions, "-")
    expect_equal(parsed$C1.start, 1328183)
    expect_equal(parsed$C1.end,   1328170)
    expect_equal(parsed$C2.start, 1326955)
    expect_equal(parsed$C2.end,   1326146)
})

test_that("parseSuppaA3SS parses an alt. 3' splice site event junctions", {
    # Strand plus
    junctions <- read.table(text = "169772450 169773216 169772450 169773253")
    parsed <- parseSuppaA3SS(junctions, "+")
    expect_equal(parsed$C1.end,   169772450)
    expect_equal(parsed$A1.start, 169773216)
    expect_equal(parsed$A2.start, 169773253)
    
    # Strand minus
    junctions <- read.table(text = "1019466 1019733 1019391 1019733")
    parsed <- parseSuppaA3SS(junctions, "-")
    expect_equal(parsed$C1.end,   1019733)
    expect_equal(parsed$A1.start, 1019391)
    expect_equal(parsed$A2.start, 1019466)
})

test_that("parseSuppaA5SS parses an alt. 5' splice site event junctions", {
    # Strand plus
    junctions <- read.table(text = "50193276 50197008 50192997 50197008")
    parsed <- parseSuppaA5SS(junctions, "+")
    expect_equal(parsed$A2.end,   50192997)
    expect_equal(parsed$A1.end,   50193276)
    expect_equal(parsed$C2.start, 50197008)
    
    # Strand minus
    junctions <- read.table(text = "29543197 29547229 29543197 29547350")
    parsed <- parseSuppaA5SS(junctions, "-")
    expect_equal(parsed$A2.end,   29547229)
    expect_equal(parsed$A1.end,   29547350)
    expect_equal(parsed$C2.start, 29543197)
})

test_that("parseSuppaAFE parses an alt. first exon event junctions", {
    # Strand plus
    junctions <- read.table(
        text = "169763871 169764046 169767998 169764550 169765124 169767998")
    parsed <- parseSuppaAFE(junctions, "+")
    expect_equal(parsed$A2.start, 169764550)
    expect_equal(parsed$A2.end,   169765124)
    expect_equal(parsed$A1.start, 169763871)
    expect_equal(parsed$A1.end,   169764046)
    expect_equal(parsed$C2.start, 169767998)
    
    # Strand minus
    junctions <- read.table(
        text = "2341890 2343830 2344010 2341890 2345036 2345236")
    parsed <- parseSuppaAFE(junctions, "-")
    expect_equal(parsed$A2.start, 2344010)
    expect_equal(parsed$A2.end,   2343830)
    expect_equal(parsed$A1.start, 2345236)
    expect_equal(parsed$A1.end,   2345036)
    expect_equal(parsed$C2.start, 2341890)
})

test_that("parseSuppaALE parses an alt. last exon event junctions", {
    # Strand plus
    junctions <- read.table(
        text = "24790610 24792494 24792800 24790610 24795476 24795797")
    parsed <- parseSuppaALE(junctions, "+")
    expect_equal(parsed$C1.end,   24790610)
    expect_equal(parsed$A1.start, 24792494)
    expect_equal(parsed$A1.end,   24792800)
    expect_equal(parsed$A2.start, 24795476)
    expect_equal(parsed$A2.end,   24795797)
    
    # Strand minus
    junctions <- read.table(
        text = "31478894 31479307 31504938 31496482 31496997 31504938")
    parsed <- parseSuppaALE(junctions, "-")
    expect_equal(parsed$C1.end,   31504938)
    expect_equal(parsed$A1.start, 31496997)
    expect_equal(parsed$A1.end,   31496482)
    expect_equal(parsed$A2.start, 31479307)
    expect_equal(parsed$A2.end,   31478894)
})
