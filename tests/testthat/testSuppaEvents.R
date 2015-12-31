context("Parse SUPPA splicing events")

test_that("parseSuppaEvent parses multiple skipping exon event IDs at once", {
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
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153768260-153768465:153768260-153768553:-")
    expect_silent(parsed <- parseSuppaEvent(events))
    expect_is(parsed, "data.frame")
    # number of elements in list is the same as number of events
    expect_equal(nrow(parsed), length(events))
    expect_equal(parsed$Program, rep("SUPPA", 4))
    expect_equal(parsed$Event.type, rep("A3SS", 4))
    expect_equal(parsed$C1.end, c("115086707", "153768260", "153769180",
                                  "153768553"))
})

test_that("parseSuppaSE parses a skipping exon event junctions", {
    junctions <- read.table(text = "169768099 169770024 169770112 169771762")
    parsed <- parseSuppaSE(junctions, "+")
    expect_equal(parsed$C1.end,   169768099)
    expect_equal(parsed$A1.start, 169770024)
    expect_equal(parsed$A1.end,   169770112)
    expect_equal(parsed$C2.start, 169771762)
})

test_that("parseSuppaMXE parses a mutually exclusive exon event junctions", {
    junctions <- read.table(text = "202060671 202068453 202068489 202073793 202060671 202072798 202072906 202073793")
    parsed <- parseSuppaMXE(junctions, "+")
    expect_equal(parsed$C1.end,   202060671)
    expect_equal(parsed$A1.start, 202068453)
    expect_equal(parsed$A1.end,   202068489)
    expect_equal(parsed$A2.start, 202072798)
    expect_equal(parsed$A2.end,   202072906)
    expect_equal(parsed$C2.start, 202073793)
})

test_that("parseSuppaRI parses an intron retention event junctions", {
    junctions <- read.table(text = "196709749 196709922 196711005 196711181")
    parsed <- parseSuppaRI(junctions, "+")
    expect_equal(parsed$C1.start, 196709749)
    expect_equal(parsed$C1.end,   196709922)
    expect_equal(parsed$C2.start, 196711005)
    expect_equal(parsed$C2.end,   196711181)
})

test_that("parseSuppaA3SS parses an alt. 3' splice site event junctions", {
    junctions <- read.table(text = "169772450 169773216 169772450 169773253")
    parsed <- parseSuppaA3SS(junctions, "+")
    expect_equal(parsed$C1.end,   169772450)
    expect_equivalent(parsed$C2.start, list(c(169773216, 169773253)))
})

test_that("parseSuppaA5SS parses an alt. 5' splice site event junctions", {
    junctions <- read.table(text = "50193276 50197008 50192997 50197008")
    parsed <- parseSuppaA5SS(junctions, "+")
    expect_equivalent(parsed$C1.end, list(c(50193276, 50192997)))
    expect_equal(parsed$C2.start, 50197008)
})

test_that("parseSuppaAFE parses an alt. first exon event junctions", {
    junctions <- read.table(text = "169763871 169764046 169767998 169764550 169765124 169767998")
    parsed <- parseSuppaAFE(junctions, "+")
    expect_equal(parsed$C1.start, 169763871)
    expect_equal(parsed$C1.end,   169764046)
    expect_equal(parsed$A1.start, 169764550)
    expect_equal(parsed$A1.end,   169765124)
    expect_equal(parsed$C2.start, 169767998)
})

test_that("parseSuppaALE parses an alt. last exon event junctions", {
    junctions <- read.table(text = "24790610 24792494 24792800 24790610 24795476 24795797")
    parsed <- parseSuppaALE(junctions, "+")
    expect_equal(parsed$C1.end,   24790610)
    expect_equal(parsed$A1.start, 24795476)
    expect_equal(parsed$A1.end,   24795797)
    expect_equal(parsed$C2.start, 24792494)
    expect_equal(parsed$C2.end,   24792800)
})