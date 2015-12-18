context("Parse SUPPA splicing events")
## parseSuppaEvent is tested by testing parseSuppaEventID

test_that("parseSuppaJunctions parses alt. 3' splice site (+ strand)", {
    junctions <- c("169772450-169773216", "169772450-169773253")
    parsed <- parseSuppaJunctions("A3SS", "+", junctions)
    expect_equal(parsed$"C1 end",     169772450)
    expect_equal(parsed$"C2 start", c(169773216, 169773253))
})

test_that("parseSuppaJunctions parses alt. 3' splice site (- strand)", {
    junctions <- c("49557492-49557642", "49557470-49557642")
    parsed <- parseSuppaJunctions("A3SS", "-", junctions)
    expect_equal(parsed$"C1 end",     49557642)
    expect_equal(parsed$"C2 start", c(49557492, 49557470))
})

test_that("parseSuppaJunctions parses alt. 5' splice site (+ strand)", {
    junctions <- c("50193276-50197008", "50192997-50197008")
    parsed <- parseSuppaJunctions("A5SS", "+", junctions)
    expect_equal(parsed$"C1 end", c(50193276, 50192997))
    expect_equal(parsed$"C2 start", 50197008)
})

test_that("parseSuppaJunctions parses alt. 5' splice site (- strand)", { 
    junctions <- c("99890743-99891188", "99890743-99891605")
    parsed <- parseSuppaJunctions("A5SS", "-", junctions)
    expect_equal(parsed$"C1 end", c(99891188, 99891605))
    expect_equal(parsed$"C2 start", 99890743)
})

test_that("parseSuppaJunctions parses alt. first exon (+ strand)", {
    junctions <- c("169763871", "169764046-169767998", "169764550",
                   "169765124-169767998")
    parsed <- parseSuppaJunctions("AFE", "+", junctions)
    expect_equal(parsed$`C1 start`, 169763871)
    expect_equal(parsed$`C1 end`,   169764046)
    expect_equal(parsed$`A1 start`, 169764550)
    expect_equal(parsed$`A1 end`,   169765124)
    expect_equal(parsed$`C2 start`, 169767998)
})

test_that("parseSuppaJunctions parses alt. first exon (- strand)", {
    junctions <- c("169858031-169862929", "169863076", "169858031-169863148",
                   "169863408")
    parsed <- parseSuppaJunctions("AFE", "-", junctions)
    expect_equal(parsed$`C1 start`, 169863408)
    expect_equal(parsed$`C1 end`,   169863148)
    expect_equal(parsed$`A1 start`, 169863076)
    expect_equal(parsed$`A1 end`,   169862929)
    expect_equal(parsed$`C2 start`, 169858031)
})

test_that("parseSuppaJunctions parses alt. last exon (+ strand)", {
    junctions <- c("24790610-24792494", "24792800", "24790610-24795476",
                   "24795797")
    parsed <- parseSuppaJunctions("ALE", "+", junctions)
    expect_equal(parsed$`C1 end`,   24790610)
    expect_equal(parsed$`A1 start`, 24795476)
    expect_equal(parsed$`A1 end`,   24795797)
    expect_equal(parsed$`C2 start`, 24792494)
    expect_equal(parsed$`C2 end`,   24792800)
})

test_that("parseSuppaJunctions parses alt. last exon (- strand)", {
    junctionsAL <- c("64037473", "64037809-64051654", "64044233",
                     "64044515-64051654")
    parsed <- parseSuppaJunctions("ALE", "-", junctionsAL)
    expect_equal(parsed$`C1 end`,   64051654)
    expect_equal(parsed$`A1 start`, 64037809)
    expect_equal(parsed$`A1 end`,   64037473)
    expect_equal(parsed$`C2 start`, 64044515)
    expect_equal(parsed$`C2 end`,   64044233)
})

test_that("parseSuppaJunctions parses skipping exon (+ strand)", {
    junctionsSE <- c("169768099-169770024", "169770112-169771762")
    parsed <- parseSuppaJunctions("SE", "+", junctionsSE)
    expect_equal(parsed$`C1 end`,   169768099)
    expect_equal(parsed$`A1 start`, 169770024)
    expect_equal(parsed$`A1 end`,   169770112)
    expect_equal(parsed$`C2 start`, 169771762)
})

test_that("parseSuppaJunctions parses skipping exon (- strand)", {
    junctionsSE <- c("49557470-49557642", "49557746-49558568")
    parsed <- parseSuppaJunctions("SE", "-", junctionsSE)
    expect_equal(parsed$`C1 end`,   49558568)
    expect_equal(parsed$`A1 start`, 49557746)
    expect_equal(parsed$`A1 end`,   49557642)
    expect_equal(parsed$`C2 start`, 49557470)
})

test_that("parseSuppaJunctions parses mutually excl. exon (+ strand)", {
    junctions <- c("202060671-202068453", "202068489-202073793",
                   "202060671-202072798", "202072906-202073793")
    parsed <- parseSuppaJunctions("MXE", "+", junctions)
    expect_equal(parsed$`C1 end`,   202060671)
    expect_equal(parsed$`A1 start`, 202068453)
    expect_equal(parsed$`A1 end`,   202068489)
    expect_equal(parsed$`A2 start`, 202072798)
    expect_equal(parsed$`A2 end`,   202072906)
    expect_equal(parsed$`C2 start`, 202073793)
})

test_that("parseSuppaJunctions parses mutually excl. exon (- strand)", {
    junctions <- c("49557470-49557666", "49557746-49562274",
                   "49557470-49558568", "49558663-49562274")
    parsed <- parseSuppaJunctions("MXE", "-", junctions)
    expect_equal(parsed$`C1 end`,   49562274)
    expect_equal(parsed$`A1 start`, 49557746)
    expect_equal(parsed$`A1 end`,   49557666)
    expect_equal(parsed$`A2 start`, 49558663)
    expect_equal(parsed$`A2 end`,   49558568)
    expect_equal(parsed$`C2 start`, 49557470)
})

test_that("parseSuppaJunctions parses retained intron (+ strand)", {
    junctions <- c("196709749", "196709922-196711005", "196711181")
    parsed <- parseSuppaJunctions("RI", "+", junctions)
    expect_equal(parsed$`C1 start`, 196709749)
    expect_equal(parsed$`C1 end`,   196709922)
    expect_equal(parsed$`C2 start`, 196711005)
    expect_equal(parsed$`C2 end`,   196711181)
})

test_that("parseSuppaJunctions parses retained intron (- strand)", {
    junctions <- c("1038930", "1039052-1039217", "1039310")
    parsed <- parseSuppaJunctions("RI", "-", junctions)
    expect_equal(parsed$`C1 start`, 1039310)
    expect_equal(parsed$`C1 end`,   1039217)
    expect_equal(parsed$`C2 start`, 1039052)
    expect_equal(parsed$`C2 end`,   1038930)
})

test_that("parseSuppaEventID parses multiple event IDs at once", {
    # Load all types of events to test
    eventID_A3 <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
    eventID_A5 <- "ENSG00000000003;A5:X:99890743-99891188:99890743-99891605:-"
    eventID_AF <- "ENSG00000000457;AF:1:169858031-169862929:169863076:169858031-169863148:169863408:-"
    eventID_AL <- "ENSG00000001461;AL:1:24790610-24792494:24792800:24790610-24795476:24795797:+"
    eventID_SE <- "ENSG00000000419;SE:20:49557470-49557642:49557746-49558568:-"
    eventID_MX <- "ENSG00000000419;MX:20:49557470-49557666:49557746-49562274:49557470-49558568:49558663-49562274:-"
    eventID_RI <- "ENSG00000000971;RI:1:196709749:196709922-196711005:196711181:+"
    events <- c(eventID_A3, eventID_A5, eventID_AF, eventID_AL, eventID_SE,
                eventID_MX, eventID_RI)
    
    expect_silent(parsed <- parseSuppaEventID(events))
    expect_is(parsed, "list")
    # number of elements in list is the same as number of events
    expect_equal(length(parsed), length(events))
})

test_that("parseSuppaSE parses a skipping exon event junctions", {
    junctions <- c(169768099, 169770024, 169770112, 169771762)
    parsed <- parseSuppaSE(junctions, "+")
    expect_equal(parsed$`C1 end`,   169768099)
    expect_equal(parsed$`A1 start`, 169770024)
    expect_equal(parsed$`A1 end`,   169770112)
    expect_equal(parsed$`C2 start`, 169771762)
})

test_that("parseSuppaMXE parses a mutually exclusive exon event junctions", {
    junctions <- c(202060671, 202068453, 202068489, 202073793,
                   202060671, 202072798, 202072906, 202073793)
    parsed <- parseSuppaMXE(junctions, "+")
    expect_equal(parsed$`C1 end`,   202060671)
    expect_equal(parsed$`A1 start`, 202068453)
    expect_equal(parsed$`A1 end`,   202068489)
    expect_equal(parsed$`A2 start`, 202072798)
    expect_equal(parsed$`A2 end`,   202072906)
    expect_equal(parsed$`C2 start`, 202073793)
})

test_that("parseSuppaRI parses an intron retention event junctions", {
    junctions <- c(196709749, 196709922, 196711005, 196711181)
    parsed <- parseSuppaRI(junctions, "+")
    expect_equal(parsed$`C1 start`, 196709749)
    expect_equal(parsed$`C1 end`,   196709922)
    expect_equal(parsed$`C2 start`, 196711005)
    expect_equal(parsed$`C2 end`,   196711181)
})

test_that("parseSuppaA3SS parses an alt. 3' splice site event junctions", {
    junctions <- c(169772450, 169773216, 169772450, 169773253)
    parsed <- parseSuppaA3SS(junctions, "+")
    expect_equal(parsed$`C1 end`,   169772450)
    expect_equal(parsed$`C2 start`, c(169773216, 169773253))
})

test_that("parseSuppaA5SS parses an alt. 5' splice site event junctions", {
    junctions <- c(50193276, 50197008, 50192997, 50197008)
    parsed <- parseSuppaA5SS(junctions, "+")
    expect_equal(parsed$`C1 end`,   c(50193276, 50192997))
    expect_equal(parsed$`C2 start`, 50197008)
})

test_that("parseSuppaAFE parses an alt. first exon event junctions", {
    junctions <- c(169763871, 169764046, 169767998, 169764550, 169765124, 
                   169767998)
    parsed <- parseSuppaAFE(junctions, "+")
    expect_equal(parsed$`C1 start`, 169763871)
    expect_equal(parsed$`C1 end`,   169764046)
    expect_equal(parsed$`A1 start`, 169764550)
    expect_equal(parsed$`A1 end`,   169765124)
    expect_equal(parsed$`C2 start`, 169767998)
})

test_that("parseSuppaALE parses an alt. last exon event junctions", {
    junctions <- c(24790610, 24792494, 24792800, 24790610, 24795476, 24795797)
    parsed <- parseSuppaALE(junctions, "+")
    expect_equal(parsed$`C1 end`,   24790610)
    expect_equal(parsed$`A1 start`, 24795476)
    expect_equal(parsed$`A1 end`,   24795797)
    expect_equal(parsed$`C2 start`, 24792494)
    expect_equal(parsed$`C2 end`,   24792800)
})