context("Parse SUPPA splicing events")
# source("R/readPSIfiles.R")

# Load all types of events to test
eventA3 <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
eventA5 <- "ENSG00000000003;A5:X:99890743-99891188:99890743-99891605:-"
eventAF <- "ENSG00000000457;AF:1:169858031-169862929:169863076:169858031-169863148:169863408:-"
eventAL <- "ENSG00000001461;AL:1:24790610-24792494:24792800:24790610-24795476:24795797:+"
eventSE <- "ENSG00000000419;SE:20:49557470-49557642:49557746-49558568:-"
eventMX <- "ENSG00000000419;MX:20:49557470-49557666:49557746-49562274:49557470-49558568:49558663-49562274:-"
eventRI <- "ENSG00000000971;RI:1:196709749:196709922-196711005:196711181:+"

# Load types of junctions
junctionsA3 <- c("49557492-49557642", "49557470-49557642")
junctionsA5 <- c("99890743-99891188", "99890743-99891605")
junctionsAF <- c("169858031-169862929", "169863076", "169858031-169863148",
                 "169863408")
junctionsAL <- c("24790610-24792494", "24792800", "24790610-24795476",
                 "24795797")
junctionsSE <- c("49557470-49557642", "49557746-49558568")
junctionsMX <- c("49557470-49557666", "49557746-49562274", "49557470-49558568", 
                 "49558663-49562274")
junctionsRI <- c("196709749", "196709922-196711005", "196711181")

test_that("parseSuppaJunctions correctly parses junctions by event type", {
  a <- parseSuppaJunctions("A3", "-", junctionsA3)
  expect_equal(a$"C1 end", "49557642")
  expect_equal(a$"C2 start", c("49557470", "49557492"))
  
  ## TODO: continue the tests
  parseSuppaJunctions("A5", "-", junctionsA5)
  parseSuppaJunctions("AF", "-", junctionsAF)
  parseSuppaJunctions("AL", "+", junctionsAL)
  parseSuppaJunctions("SE", "-", junctionsSE)
  parseSuppaJunctions("MX", "-", junctionsMX)
  parseSuppaJunctions("RI", "+", junctionsRI)
})

parsed_eventA3 <- parseSuppaEvent(eventA3)
parsed_eventA5 <- parseSuppaEvent(eventA5)
parsed_eventAF <- parseSuppaEvent(eventAF)
parsed_eventAL <- parseSuppaEvent(eventAL)
parsed_eventSE <- parseSuppaEvent(eventSE)
parsed_eventMX <- parseSuppaEvent(eventMX)
parsed_eventRI <- parseSuppaEvent(eventRI)