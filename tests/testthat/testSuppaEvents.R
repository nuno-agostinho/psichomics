# source("R/readPSIfiles.R")
eventA3 <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
eventA5 <- "ENSG00000000003;A5:X:99890743-99891188:99890743-99891605:-"
eventAF <- "ENSG00000000457;AF:1:169858031-169862929:169863076:
169858031-169863148:169863408:-"
eventAL <- "ENSG00000001461;AL:1:24790610-24792494:24792800:24790610-24795476:24795797:+"
eventSE <- "ENSG00000000419;SE:20:49557470-49557642:49557746-49558568:-"
eventMX <-
  "ENSG00000000419;MX:20:49557470-49557666:49557746-49562274:49557470-
49558568:49558663-49562274:-"
eventRI <- rownames(suppa_RI.tab)[1]

# parsed_eventA3 <- parseSuppaEventE(eventA3)
# parsed_eventA5 <- parseSuppaEventE(eventA5)
# parsed_eventAF <- parseSuppaEventE(eventAF)
# parsed_eventAL <- parseSuppaEventE(eventAL)
# parsed_eventSE <- parseSuppaEventE(eventSE)
# parsed_eventMX <- parseSuppaEventE(eventMX)
# parsed_eventRI <- parseSuppaEventE(eventRI)

a <- lapply(events, parseSuppaEventE)