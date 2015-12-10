context("Parse SUPPA splicing events")
## parseSuppaEvent is tested by testing parseSuppaEventID

test_that("parseSuppaJunctions parses alt. 3' splice site (positive strand)", {
  junctionsA3 <- c("169772450-169773216", "169772450-169773253")
  res <- parseSuppaJunctions("A3", "+", junctionsA3)
  expect_equal(res$"C1 end",     169772450)
  expect_equal(res$"C2 start", c(169773216, 169773253))
})

test_that("parseSuppaJunctions parses alt. 3' splice site (negative strand)", {
  junctionsA3 <- c("49557492-49557642", "49557470-49557642")
  res <- parseSuppaJunctions("A3", "-", junctionsA3)
  expect_equal(res$"C1 end",     49557642)
  expect_equal(res$"C2 start", c(49557492, 49557470))
})

test_that("parseSuppaJunctions parses alt. 5' splice site (positive strand)", {
  junctionsA5 <- c("50193276-50197008", "50192997-50197008")
  res <- parseSuppaJunctions("A5", "+", junctionsA5)
  expect_equal(res$"C1 end", c(50193276, 50192997))
  expect_equal(res$"C2 start", 50197008)
})

test_that("parseSuppaJunctions parses alt. 5' splice site (negative strand)", { 
  junctionsA5 <- c("99890743-99891188", "99890743-99891605")
  res <- parseSuppaJunctions("A5", "-", junctionsA5)
  expect_equal(res$"C1 end", c(99891188, 99891605))
  expect_equal(res$"C2 start", 99890743)
})

test_that("parseSuppaJunctions parses alt. first exon (positive strand)", {
  junctionsAF <- c("169763871", "169764046-169767998", "169764550",
                   "169765124-169767998")
  res <- parseSuppaJunctions("AF", "+", junctionsAF)
  expect_equal(res$`C1 start`, 169763871)
  expect_equal(res$`C1 end`,   169764046)
  expect_equal(res$`A1 start`, 169764550)
  expect_equal(res$`A1 end`,   169765124)
  expect_equal(res$`C2 start`, 169767998)
})

test_that("parseSuppaJunctions parses alt. first exon (negative strand)", {
  junctionsAF <- c("169858031-169862929", "169863076", "169858031-169863148",
                   "169863408")
  res <- parseSuppaJunctions("AF", "-", junctionsAF)
  expect_equal(res$`C1 start`, 169863408)
  expect_equal(res$`C1 end`,   169863148)
  expect_equal(res$`A1 start`, 169863076)
  expect_equal(res$`A1 end`,   169862929)
  expect_equal(res$`C2 start`, 169858031)
})

test_that("parseSuppaJunctions parses alt. last exon (positive strand)", {
  junctionsAL <- c("24790610-24792494", "24792800", "24790610-24795476",
                   "24795797")
  res <- parseSuppaJunctions("AL", "+", junctionsAL)
  expect_equal(res$`C1 end`,   24790610)
  expect_equal(res$`A1 start`, 24795476)
  expect_equal(res$`A1 end`,   24795797)
  expect_equal(res$`C2 start`, 24792494)
  expect_equal(res$`C2 end`,   24792800)
})

test_that("parseSuppaJunctions parses alt. last exon (negative strand)", {
  junctionsAL <- c("64037473", "64037809-64051654", "64044233",
                   "64044515-64051654")
  res <- parseSuppaJunctions("AL", "-", junctionsAL)
  expect_equal(res$`C1 end`,   64051654)
  expect_equal(res$`A1 start`, 64037809)
  expect_equal(res$`A1 end`,   64037473)
  expect_equal(res$`C2 start`, 64044515)
  expect_equal(res$`C2 end`,   64044233)
})

test_that("parseSuppaJunctions parses skipping exon (positive strand)", {
  junctionsSE <- c("169768099-169770024", "169770112-169771762")
  res <- parseSuppaJunctions("SE", "+", junctionsSE)
  expect_equal(res$`C1 end`,   169768099)
  expect_equal(res$`A1 start`, 169770024)
  expect_equal(res$`A1 end`,   169770112)
  expect_equal(res$`C2 start`, 169771762)
})

test_that("parseSuppaJunctions parses skipping exon (negative strand)", {
  junctionsSE <- c("49557470-49557642", "49557746-49558568")
  res <- parseSuppaJunctions("SE", "-", junctionsSE)
  expect_equal(res$`C1 end`,   49558568)
  expect_equal(res$`A1 start`, 49557746)
  expect_equal(res$`A1 end`,   49557642)
  expect_equal(res$`C2 start`, 49557470)
})

test_that("parseSuppaJunctions parses mutually excl. exon (positive strand)", {
  junctionsMX <- c("202060671-202068453", "202068489-202073793",
                   "202060671-202072798", "202072906-202073793")
  res <- parseSuppaJunctions("MX", "+", junctionsMX)
  expect_equal(res$`C1 end`,   202060671)
  expect_equal(res$`A1 start`, 202068453)
  expect_equal(res$`A1 end`,   202068489)
  expect_equal(res$`A2 start`, 202072798)
  expect_equal(res$`A2 end`,   202072906)
  expect_equal(res$`C2 start`, 202073793)
})

test_that("parseSuppaJunctions parses mutually excl. exon (negative strand)", {
  junctionsMX <- c("49557470-49557666", "49557746-49562274",
                   "49557470-49558568", "49558663-49562274")
  res <- parseSuppaJunctions("MX", "-", junctionsMX)
  expect_equal(res$`C1 end`,   49562274)
  expect_equal(res$`A1 start`, 49557746)
  expect_equal(res$`A1 end`,   49557666)
  expect_equal(res$`A2 start`, 49558663)
  expect_equal(res$`A2 end`,   49558568)
  expect_equal(res$`C2 start`, 49557470)
})

test_that("parseSuppaJunctions parses retained intron (positive strand)", {
  junctionsRI <- c("196709749", "196709922-196711005", "196711181")
  res <- parseSuppaJunctions("RI", "+", junctionsRI)
  expect_equal(res$`C1 start`, 196709749)
  expect_equal(res$`C1 end`,   196709922)
  expect_equal(res$`C2 start`, 196711005)
  expect_equal(res$`C2 end`,   196711181)
})

test_that("parseSuppaJunctions parses retained intron (negative strand)", {
  junctionsRI <- c("1038930", "1039052-1039217", "1039310")
  res <- parseSuppaJunctions("RI", "-", junctionsRI)
  expect_equal(res$`C1 start`, 1039310)
  expect_equal(res$`C1 end`,   1039217)
  expect_equal(res$`C2 start`, 1039052)
  expect_equal(res$`C2 end`,   1038930)
})

test_that("parseSuppaEventID parses multiple event types at once", {
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