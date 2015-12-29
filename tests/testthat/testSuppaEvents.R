context("Parse SUPPA splicing events")

test_that("parseSuppaEventID parses multiple event IDs at once", {
    # Load all types of events to test
    events <- c(
        "ENSG00000260117;A3:HG531_PATCH:115086707-115089263:115086707-115089266:+ ",
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153768260-153768394:153768260-153768553:+",
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153769180-153769533:153769180-153769638:+",
        "ENSG00000263163;A3:HSCHR1_1_CTG31:153768260-153768465:153768260-153768553:+")
    expect_silent(parsed <- parseSuppaEvent(events))
    expect_is(parsed, "data.frame")
    # number of elements in list is the same as number of events
    expect_equal(nrow(parsed), length(events))
    expect_equal(parsed$Program, rep("SUPPA", 4))
    expect_equal(parsed$Event.type, rep("A3SS", 4))
    expect_equal(parsed$C1.end, c("115089266", "153768553", "153769638",
                                  "153768553"))
})

test_that("parseSuppaJunctions parses alt. 3' splice site (+ strand)", {
    junctions <- read.table(text="169772450 169773216 169772450 169773253")
    parsed <- parseSuppaJunctions("A3SS", "+", junctions)
    expect_equal(parsed$C1.end, 169772450)
    expect_equivalent(parsed$C2.start, list(c(169773216, 169773253)))
})

test_that("parseSuppaJunctions parses alt. 3' splice site (- strand)", {
    junctions <- read.table(text="49557492 49557642 49557470 49557642")
    parsed <- parseSuppaJunctions("A3SS", "-", junctions)
    expect_equal(parsed$C1.end,     49557642)
    expect_equivalent(parsed$C2.start, list(c(49557492, 49557470)))
})

test_that("parseSuppaJunctions parses alt. 5' splice site (+ strand)", {
    junctions <- read.table(text="50193276 50197008 50192997 50197008")
    parsed <- parseSuppaJunctions("A5SS", "+", junctions)
    expect_equivalent(parsed$C1.end, list(c(50193276, 50192997)))
    expect_equal(parsed$C2.start, 50197008)
})

test_that("parseSuppaJunctions parses alt. 5' splice site (- strand)", { 
    junctions <- read.table(text="99890743 99891188 99890743 99891605")
    parsed <- parseSuppaJunctions("A5SS", "-", junctions)
    expect_equivalent(parsed$C1.end, list(c(99891188, 99891605)))
    expect_equal(parsed$C2.start, 99890743)
})

test_that("parseSuppaJunctions parses alt. first exon (+ strand)", {
    junctions <- read.table(
        text="169763871 169764046 169767998 169764550 169765124 169767998")
    parsed <- parseSuppaJunctions("AFE", "+", junctions)
    expect_equal(parsed$C1.start, 169763871)
    expect_equal(parsed$C1.end,   169764046)
    expect_equal(parsed$A1.start, 169764550)
    expect_equal(parsed$A1.end,   169765124)
    expect_equal(parsed$C2.start, 169767998)
})

test_that("parseSuppaJunctions parses alt. first exon (- strand)", {
    junctions <- read.table(
        text="169858031 169862929 169863076 169858031 169863148 169863408")
    parsed <- parseSuppaJunctions("AFE", "-", junctions)
    expect_equal(parsed$C1.start, 169863408)
    expect_equal(parsed$C1.end,   169863148)
    expect_equal(parsed$A1.start, 169863076)
    expect_equal(parsed$A1.end,   169862929)
    expect_equal(parsed$C2.start, 169858031)
})

test_that("parseSuppaJunctions parses alt. last exon (+ strand)", {
    junctions <- read.table(
        text="24790610 24792494 24792800 24790610 24795476 24795797")
    parsed <- parseSuppaJunctions("ALE", "+", junctions)
    expect_equal(parsed$C1.end,   24790610)
    expect_equal(parsed$A1.start, 24795476)
    expect_equal(parsed$A1.end,   24795797)
    expect_equal(parsed$C2.start, 24792494)
    expect_equal(parsed$C2.end,   24792800)
})

test_that("parseSuppaJunctions parses alt. last exon (- strand)", {
    junctions <- read.table(
        text="64037473 64037809 64051654 64044233 64044515 64051654")
    parsed <- parseSuppaJunctions("ALE", "-", junctions)
    expect_equal(parsed$C1.end,   64051654)
    expect_equal(parsed$A1.start, 64037809)
    expect_equal(parsed$A1.end,   64037473)
    expect_equal(parsed$C2.start, 64044515)
    expect_equal(parsed$C2.end,   64044233)
})

test_that("parseSuppaJunctions parses skipping exon (+ strand)", {
    junctions <- read.table(
        text="169768099 169770024 169770112 169771762")
    parsed <- parseSuppaJunctions("SE", "+", junctions)
    expect_equal(parsed$C1.end,   169768099)
    expect_equal(parsed$A1.start, 169770024)
    expect_equal(parsed$A1.end,   169770112)
    expect_equal(parsed$C2.start, 169771762)
})

test_that("parseSuppaJunctions parses skipping exon (- strand)", {
    junctions <- read.table(
        text="49557470 49557642 49557746 49558568")
    parsed <- parseSuppaJunctions("SE", "-", junctions)
    expect_equal(parsed$C1.end,   49558568)
    expect_equal(parsed$A1.start, 49557746)
    expect_equal(parsed$A1.end,   49557642)
    expect_equal(parsed$C2.start, 49557470)
})

test_that("parseSuppaJunctions parses mutually excl. exon (+ strand)", {
    junctions <- read.table(
        text="202060671 202068453 202068489 202073793 202060671 202072798 202072906 202073793")
    parsed <- parseSuppaJunctions("MXE", "+", junctions)
    expect_equal(parsed$C1.end,   202060671)
    expect_equal(parsed$A1.start, 202068453)
    expect_equal(parsed$A1.end,   202068489)
    expect_equal(parsed$A2.start, 202072798)
    expect_equal(parsed$A2.end,   202072906)
    expect_equal(parsed$C2.start, 202073793)
})

test_that("parseSuppaJunctions parses mutually excl. exon (- strand)", {
    junctions <- read.table(
        text="49557470 49557666 49557746 49562274 49557470 49558568 49558663 49562274")
    parsed <- parseSuppaJunctions("MXE", "-", junctions)
    expect_equal(parsed$C1.end,   49562274)
    expect_equal(parsed$A1.start, 49557746)
    expect_equal(parsed$A1.end,   49557666)
    expect_equal(parsed$A2.start, 49558663)
    expect_equal(parsed$A2.end,   49558568)
    expect_equal(parsed$C2.start, 49557470)
})

test_that("parseSuppaJunctions parses retained intron (+ strand)", {
    junctions <- read.table(text="196709749 196709922 196711005 196711181")
    parsed <- parseSuppaJunctions("RI", "+", junctions)
    expect_equal(parsed$C1.start, 196709749)
    expect_equal(parsed$C1.end,   196709922)
    expect_equal(parsed$C2.start, 196711005)
    expect_equal(parsed$C2.end,   196711181)
})

test_that("parseSuppaJunctions parses retained intron (- strand)", {
    junctions <- read.table(text="1038930 1039052 1039217 1039310")
    parsed <- parseSuppaJunctions("RI", "-", junctions)
    expect_equal(parsed$C1.start, 1039310)
    expect_equal(parsed$C1.end,   1039217)
    expect_equal(parsed$C2.start, 1039052)
    expect_equal(parsed$C2.end,   1038930)
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