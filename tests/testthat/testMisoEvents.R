context("Parse MISO splicing events")
library(fastmatch)

annotation <- read.table(text="
    chr19 AFE gene 50015886 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE mRNA 50028313 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1.A;Parent=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1.A;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE exon 50028313 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1.A.0;Parent=2217@uc002poi.1@uc002poe.1.A;Name=2217@uc002poi.1@uc002poe.1.A.0;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE mRNA 50015886 50016730  .  +  . ID=2217@uc002poi.1@uc002poe.1.B;Parent=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1.B;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE exon 50015886 50016007  .  +  . ID=2217@uc002poi.1@uc002poe.1.B.0;Parent=2217@uc002poi.1@uc002poe.1.B;Name=2217@uc002poi.1@uc002poe.1.B.0;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE exon 50016644 50016730  .  +  . ID=2217@uc002poi.1@uc002poe.1.B.1;Parent=2217@uc002poi.1@uc002poe.1.B;Name=2217@uc002poi.1@uc002poe.1.B.1;gid=2217@uc002poi.1@uc002poe.1
    chr19 AFE gene 50016492 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1;gid=2217@uc002poh.1@uc002pog.1
    chr19 AFE mRNA 50017468 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1.A;Parent=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1.A;gid=2217@uc002poh.1@uc002pog.1
    chr19 AFE exon 50017468 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1.A.0;Parent=2217@uc002poh.1@uc002pog.1.A;Name=2217@uc002poh.1@uc002pog.1.A.0;gid=2217@uc002poh.1@uc002pog.1
    chr19 AFE mRNA 50016492 50016730  .  +  . ID=2217@uc002poh.1@uc002pog.1.B;Parent=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1.B;gid=2217@uc002poh.1@uc002pog.1
    chr19 AFE exon 50016492 50016730  .  +  . ID=2217@uc002poh.1@uc002pog.1.B.0;Parent=2217@uc002poh.1@uc002pog.1.B;Name=2217@uc002poh.1@uc002pog.1.B.0;gid=2217@uc002poh.1@uc002pog.1")

event <- c(1, NA, 7)
next_event <- c(6, NA, NA)

test_that("parseMultipleMisoEvents parses two different events", {
    parsed <- parseMultipleMisoEvents(annotation)
    expect_is(parsed, "list")
    expect_equal(length(parsed), 2)
})

test_that("parseMultipleMisoEvents parses one event", {
    parsed <- parseMultipleMisoEvents(annotation[1:6, ])
    expect_is(parsed, "list")
    expect_equal(length(parsed), 1)
})

test_that("getDataRows retrieves rows of a data frame", {
    ret <- getDataRows(1, annotation, event, next_event)
    expect_is(ret, "data.frame")
    expect_equal(ret, annotation[1:6, 1:8])
    
    ret <- getDataRows(3, annotation, event, next_event)
    expect_is(ret, "data.frame")
    expect_equal(ret, annotation[7:11, 1:8])
})

test_that("getDataRows returns NA if value is outside data frame", {
    ret <- getDataRows(2, annotation, event, next_event)
    expect_true(is.na(ret))
})

test_that("parseMisoEventID matches events ID (returns NA if not possible)", {
    eventID <- c("2217@uc002poi.1@uc002poe.1", "57705@uc009xob.1@uc001jgy.2", 
                 "2217@uc002poh.1@uc002pog.1")
    columnID <- 9
    
    events <- parseMisoEventID(eventID, annotation, columnID)
    expect_equal(length(events), 3)
    expect_equal(events[[1]], annotation[1:6, 1:8])
    expect_true(is.na(events[[2]]))
    expect_equal(events[[3]], annotation[7:11, 1:8])
})

test_that("parseMisoEvent parses alternative splicing event", {
    event <- read.table(text = "
    chr1 SE gene 16854	18061	. - .
    chr1 SE mRNA 16854 18061 . - .
    chr1 SE exon 16854 17055 . - .
    chr1 SE exon 17233 17742 . - .
    chr1 SE exon 17915 18061 . - .
    chr1 SE mRNA 16854 18061 . - .
    chr1 SE exon 16854 17955 . - .
    chr1 SE exon 17915 18061 . - .")
    parsed <- parseMisoEvent(event)
    expect_equal(parsed$Program, "MISO")
    expect_equal(parsed$Chromosome, "chr1")
    expect_equal(parsed$`Event type`, "SE")
    expect_equal(parsed$Strand, "-")
    # the rest of the returned value is tested through the following tests
})

test_that("parseMisoEventSE parses exon skipping junctions (+ strand)", {
    event <- read.table(text = "
                      chr1 SE gene 1370903 1378262  .  +  .
                      chr1 SE mRNA 1370903 1378262  .  +  .
                      chr1 SE exon 1370903 1371201  .  +  .
                      chr1 SE exon 1372702 1372864  .  +  .
                      chr1 SE exon 1374461 1378262  .  +  .
                      chr1 SE mRNA 1370903 1378262  .  +  .
                      chr1 SE exon 1370903 1371201  .  +  .
                      chr1 SE exon 1374461 1378262  .  +  .")
    parsed <- parseMisoSE(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 1370903)
    expect_equal(parsed$`C1 end`,   1371201)
    expect_equal(parsed$`A1 start`, 1372702)
    expect_equal(parsed$`A1 end`,   1372864)
    expect_equal(parsed$`C2 start`, 1374461)
    expect_equal(parsed$`C2 end`,   1378262)
})

test_that("parseMisoEventSE parses exon skipping junctions (- strand)", {
    event <- read.table(text = "
    chr1 SE gene 16854 18061 . - .
    chr1 SE mRNA 16854 18061 . - .
    chr1 SE exon 16854 17055 . - .
    chr1 SE exon 17233 17742 . - .
    chr1 SE exon 17915 18061 . - .
    chr1 SE mRNA 16854 18061 . - .
    chr1 SE exon 16854 17955 . - .
    chr1 SE exon 17915 18061 . - .")
    parsed <- parseMisoSE(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 18061)
    expect_equal(parsed$`C1 end`,   17915)
    expect_equal(parsed$`A1 start`, 17742)
    expect_equal(parsed$`A1 end`,   17233)
    expect_equal(parsed$`C2 start`, 17055)
    expect_equal(parsed$`C2 end`,   16854)
})

test_that("parseMisoEventSE doesn't parse unrecognized event", {
    event <- read.table(text = "
                      chr1 SE gene 1370903 1378262  .  +  .
                      chr1 SE mRNA 1370903 1378262  .  +  .
                      chr1 SE exon 1370903 1371201  .  +  .
                      chr1 SE exon 1372702 1372864  .  +  .
                      chr1 SE exon 1374461 1378262  .  +  .
                      chr7 SE exon 1374461 1378262  .  +  .
                      chr1 SE mRNA 1370903 1378262  .  +  .
                      chr1 SE exon 1370903 1371201  .  +  .
                      chr1 SE exon 1374461 1378262  .  +  .")
    parsed <- parseMisoSE(event, strand = "+", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    #expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventMXE parses mutually exc. exons junctions (+ strand)", {
    event <- read.table(text = "
                      chr1 MXE gene 764383 788090 . + .
                      chr1 MXE mRNA 764383 788090 . + .
                      chr1 MXE exon 764383 764484 . + .
                      chr1 MXE exon 776580 776753 . + .
                      chr1 MXE exon 787307 788090 . + .
                      chr1 MXE mRNA 764383 788090 . + .
                      chr1 MXE exon 764383 764484 . + .
                      chr1 MXE exon 783034 783186 . + .
                      chr1 MXE exon 787307 788090 . + .")
    parsed <- parseMisoMXE(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 764383)
    expect_equal(parsed$`C1 end`,   764484)
    expect_equal(parsed$`A1 start`, 776580)
    expect_equal(parsed$`A1 end`,   776753)
    expect_equal(parsed$`A2 start`, 783034)
    expect_equal(parsed$`A2 end`,   783186)
    expect_equal(parsed$`C2 start`, 787307)
    expect_equal(parsed$`C2 end`,   788090)
})

test_that("parseMisoEventMXE parses mutually exc. exons junctions (- strand)", {
    event <- read.table(text = "
                      chr1 MXE gene 1027371 1051736  .  -  .
                      chr1 MXE mRNA 1027371 1051736  .  -  .
                      chr1 MXE exon 1027371 1027483  .  -  .
                      chr1 MXE exon 1041336 1041429  .  -  .
                      chr1 MXE exon 1051440 1051736  .  -  .
                      chr1 MXE mRNA 1027371 1051736  .  -  .
                      chr1 MXE exon 1027371 1027483  .  -  .
                      chr1 MXE exon 1050402 1050455  .  -  .
                      chr1 MXE exon 1051440 1051736  .  -  .")
    parsed <- parseMisoMXE(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 1051736)
    expect_equal(parsed$`C1 end`,   1051440)
    expect_equal(parsed$`A1 start`, 1041429)
    expect_equal(parsed$`A1 end`,   1041336)
    expect_equal(parsed$`A2 start`, 1050455)
    expect_equal(parsed$`A2 end`,   1050402)
    expect_equal(parsed$`C2 start`, 1027483)
    expect_equal(parsed$`C2 end`,   1027371)
})

test_that("parseMisoEventMXE doesn't parse unrecognized event", {
    event <- read.table(text = "
                      chr1 MXE gene 1027371 1051736  .  -  .
                      chr1 MXE mRNA 1027371 1051736  .  -  .
                      chr1 MXE exon 1027371 1027483  .  -  .
                      chr1 MXE exon 1041336 1041429  .  -  .
                      chr1 MXE exon 1051440 1051736  .  -  .
                      chr1 MXE mRNA 1027371 1051736  .  -  .
                      chr1 MXE exon 1027371 1027483  .  -  .
                      chr1 MXE exon 1050402 1050455  .  -  .
                      chr6 MXE gene 51440 51736  .  +  .")
    parsed <- parseMisoSE(event, strand = "+", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    #expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventRI parses intron retention junctions (+ strand)", {
    event <- read.table(text = "
                      chr1 RI gene 1223053 1223417  .  +  .
                      chr1 RI mRNA 1222888 1223216  .  +  .
                      chr1 RI exon 1222888 1223216  .  +  .
                      chr1 RI mRNA 1222888 1223216  .  +  .
                      chr1 RI exon 1222888 1222976  .  +  .
                      chr1 RI exon 1223053 1223216  .  +  .")
    parsed <- parseMisoRI(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 1222888)
    expect_equal(parsed$`C1 end`,   1222976)
    expect_equal(parsed$`C2 start`, 1223053)
    expect_equal(parsed$`C2 end`,   1223216)
})

test_that("parseMisoEventRI parses intron retention junctions (- strand)", {
    event <- read.table(text = "
                      chr1 RI gene 17233 17742 . - .
                      chr1 RI mRNA 17233 17742 . - .
                      chr1 RI exon 17233 17742 . - .
                      chr1 RI mRNA 17233 17742 . - .
                      chr1 RI exon 17233 17364 . - .
                      chr1 RI exon 17601 17742 . - .")
    parsed <- parseMisoRI(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 17742)
    expect_equal(parsed$`C1 end`,   17601)
    expect_equal(parsed$`C2 start`, 17364)
    expect_equal(parsed$`C2 end`,   17233)
})

test_that("parseMisoEventRI doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr1 RI gene 17233 17742 . - .
                      chr1 RI mRNA 17233 17742 . - .
                      chr1 RI exon 17233 17742 . - .
                      chr10 RI gene 17233 17742 . - .
                      chr10 RI mRNA 17233 17364 . - .
                      chr10 RI exon 17601 17742 . - .")
    parsed <- parseMisoRI(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventA5SS parses alt. 5' SS junctions (+ strand)", {
    event <- read.table(text = "
                      chr1 A5SS gene 874655 876686  .  +  .
                      chr1 A5SS mRNA 874655 876686  .  +  .
                      chr1 A5SS exon 874655 874840  .  +  .
                      chr1 A5SS exon 876524 876686  .  +  .
                      chr1 A5SS mRNA 874655 876686  .  +  .
                      chr1 A5SS exon 874655 874792  .  +  .
                      chr1 A5SS exon 876524 876686  .  +  .")
    parsed <- parseMisoA5SS(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 874655)
    expect_equal(parsed$`C1 end`, c(874840, 874792))
    expect_equal(parsed$`C2 start`, 876524)
    expect_equal(parsed$`C2 end`,   876686)
})

test_that("parseMisoEventA5SS parses alt. 5' SS junctions (- strand)", {
    event <- read.table(text = "
                      chr1 A5SS gene 17233 17742 . - .
                      chr1 A5SS mRNA 17233 17742 . - .
                      chr1 A5SS exon 17233 17368 . - .
                      chr1 A5SS exon 17526 17742 . - .
                      chr1 A5SS mRNA 17233 17742 . - .
                      chr1 A5SS exon 17233 17368 . - .
                      chr1 A5SS exon 17606 17742 . - .")
    parsed <- parseMisoA5SS(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 17742)
    expect_equal(parsed$`C1 end`, c(17526, 17606))
    expect_equal(parsed$`C2 start`, 17368)
    expect_equal(parsed$`C2 end`,   17233)
})

test_that("parseMisoEventA5SS doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr1 A5SS gene 17233 17742 . - .
                      chr1 A5SS mRNA 17233 17742 . - .
                      chr1 A5SS exon 17233 17368 . - .
                      chr1 A5SS mRNA 17526 17742 . - .
                      chr1 A5SS mRNA 17233 17742 . - .")
    parsed <- parseMisoA5SS(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventA3SS parses alt. 3' SS junctions (+ strand)", {
    event <- read.table(text = "
                      chr1 A3SS gene 898084 898633  .  +  .
                      chr1 A3SS mRNA 898084 898633  .  +  .
                      chr1 A3SS exon 898084 898297  .  +  .
                      chr1 A3SS exon 898412 898633  .  +  .
                      chr1 A3SS mRNA 898084 898633  .  +  .
                      chr1 A3SS exon 898084 898297  .  +  .
                      chr1 A3SS exon 898489 898633  .  +  .")
    parsed <- parseMisoA3SS(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 898084)
    expect_equal(parsed$`C1 end`,   898297)
    expect_equal(parsed$`C2 start`, c(898412, 898489))
    expect_equal(parsed$`C2 end`,   898633)
})

test_that("parseMisoEventA3SS parses alt. 3' SS junctions (- strand)", {
    event <- read.table(text = "
                      chr1 A3SS gene 15796 16765 . - .
                      chr1 A3SS mRNA 15796 16765 . - .
                      chr1 A3SS exon 15796 15947 . - .
                      chr1 A3SS exon 16607 16765 . - .
                      chr1 A3SS mRNA 15796 16765 . - .
                      chr1 A3SS exon 15796 15942 . - .
                      chr1 A3SS exon 16607 16765 . - .")
    parsed <- parseMisoA3SS(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 16765)
    expect_equal(parsed$`C1 end`,   16607)
    expect_equal(parsed$`C2 start`, c(15947, 15942))
    expect_equal(parsed$`C2 end`,   15796)
})

test_that("parseMisoEventA3SS doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr1 A3SS gene 17233 17742 . - .
                      chr1 A3SS mRNA 17233 17742 . - .
                      chr1 A3SS exon 17233 17368 . - .
                      chr1 A3SS mRNA 17526 17742 . - .
                      chr1 A3SS mRNA 17233 17742 . - .")
    parsed <- parseMisoA3SS(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})


test_that("parseMisoEventTandemUTR parses tandem UTR junctions (+ strand)", {
    event <- read.table(text = "
                      chr12 TandemUTR gene  28702106  28703099  .  +  .
                      chr12 TandemUTR mRNA  28702106  28703099  .  +  .
                      chr12 TandemUTR exon  28702106  28703099  .  +  .
                      chr12 TandemUTR mRNA  28702106  28702181  .  +  .
                      chr12 TandemUTR exon  28702106  28702181  .  +  .")
    parsed <- parseMisoTandemUTR(event, strand = "+", parsed = list())
    expect_equal(parsed$`C2 start`, 28702106)
    expect_equal(parsed$`C2 end`, c(28703099, 28702181))
})

test_that("parseMisoEventTandemUTR parses tandem UTR junctions (- strand)", {
    event <- read.table(text = "
                      chr19 TandemUTR gene  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10663759  10664625  .  -  .
                      chr19 TandemUTR exon  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .")
    parsed <- parseMisoTandemUTR(event, strand = "-", parsed = list())
    expect_equal(parsed$`C2 start`, 10664625)
    expect_equal(parsed$`C2 end`, c(10663759, 10664223))
})

test_that("parseMisoEventTandemUTR doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr1 TandemUTR gene 17233 17742 . - .
                      chr1 TandemUTR mRNA 17233 17742 . - .
                      chr1 TandemUTR exon 17233 17368 . - .
                      chr1 TandemUTR mRNA 17526 17742 . - .")
    parsed <- parseMisoTandemUTR(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventAFE parses alt. first exon junctions (+ strand)", {
    # Last exons of each mRNA group
    event <- read.table(text = "
                      chr17 AFE gene   4843303   4848515  .  +  .
                      chr17 AFE mRNA   4846721   4848515  .  +  .
                      chr17 AFE exon   4846721   4846814  .  +  .
                      chr17 AFE exon   4847853   4848515  .  +  .
                      chr17 AFE mRNA   4843303   4844288  .  +  .
                      chr17 AFE exon   4843303   4843525  .  +  .
                      chr17 AFE exon   4843782   4844021  .  +  .
                      chr17 AFE exon   4844168   4844288  .  +  .")
    parsed <- parseMisoAFE(event, strand = "+", parsed = list())
    expect_equal(parsed$`C1 start`, 4847853)
    expect_equal(parsed$`C1 end`,   4848515)
    expect_equal(parsed$`A1 start`, 4844168)
    expect_equal(parsed$`A1 end`,   4844288)
})

test_that("parseMisoEventAFE parses alt. first exon junctions (- strand)", {
    # Last exons of each mRNA group
    event <- read.table(text = "
                      chr6 AFE gene 38561740 38607924  .  -  .
                      chr6 AFE mRNA 38561740 38563843  .  -  .
                      chr6 AFE exon 38561740 38562103  .  -  .
                      chr6 AFE exon 38563546 38563843  .  -  .
                      chr6 AFE mRNA 38565686 38607924  .  -  .
                      chr6 AFE exon 38565686 38565897  .  -  .
                      chr6 AFE exon 38607576 38607924  .  -  .")
    parsed <- parseMisoAFE(event, strand = "-", parsed = list())
    expect_equal(parsed$`C1 start`, 38562103)
    expect_equal(parsed$`C1 end`,   38561740)
    expect_equal(parsed$`A1 start`, 38565897)
    expect_equal(parsed$`A1 end`,   38565686)
})

test_that("parseMisoEventAFE doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr6 AFE gene 30620579 30822593  .  -  .
                      chr6 AFE mRNA 30822190 30822593  .  -  .
                      chr6 AFE mRNA 30620579 30620982  .  -  .")
    parsed <- parseMisoALE(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})

test_that("parseMisoEventALE parses alt. last exon junctions (+ strand)", {
    # First exons of each mRNA group
    event <- read.table(text = "
                      chr19 ALE gene 7830524 7834491  .  +  .
                      chr19 ALE mRNA 7830524 7834491  .  +  .
                      chr19 ALE exon 7830524 7831693  .  +  .
                      chr19 ALE exon 7832402 7832514  .  +  .
                      chr19 ALE exon 7833724 7834491  .  +  .
                      chr19 ALE mRNA 7833724 7834491  .  +  .
                      chr19 ALE exon 7833724 7833929  .  +  .
                      chr19 ALE exon 7834259 7834491  .  +  .")
    parsed <- parseMisoALE(event, strand = "+", parsed = list())
    expect_equal(parsed$`A1 start`, 7830524)
    expect_equal(parsed$`A1 end`,   7831693)
    expect_equal(parsed$`C2 start`, 7833724)
    expect_equal(parsed$`C2 end`,   7833929)
})

test_that("parseMisoEventALE parses alt. last exon junctions (- strand)", {
    # First exons of each mRNA group
    event <- read.table(text = "
                      chr17 ALE gene 26931076 26939102  .  -  .
                      chr17 ALE mRNA 26931076 26939102  .  -  .
                      chr17 ALE exon 26931076 26932233  .  -  .
                      chr17 ALE exon 26938786 26938813  .  -  .
                      chr17 ALE exon 26939062 26939102  .  -  .
                      chr17 ALE mRNA 26934982 26938674  .  -  .
                      chr17 ALE exon 26934982 26935547  .  -  .
                      chr17 ALE exon 26938562 26938674  .  -  .")
    parsed <- parseMisoALE(event, strand = "-", parsed = list())
    expect_equal(parsed$`A1 start`, 26939102)
    expect_equal(parsed$`A1 end`,   26939062)
    expect_equal(parsed$`C2 start`, 26938674)
    expect_equal(parsed$`C2 end`,   26938562)
})

test_that("parseMisoEventALE doesn't parse unrecognized events", {
    event <- read.table(text = "
                      chr6 ALE gene 30620579 30822593  .  -  .
                      chr6 ALE mRNA 30822190 30822593  .  -  .
                      chr6 ALE mRNA 30620579 30620982  .  -  .")
    parsed <- parseMisoALE(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
    
    event <- read.table(text = "
                      chr6 ALE gene 30620579 30822593  .  -  .")
    parsed <- parseMisoALE(event, strand = "-", parsed = list())
    expect_equal(parsed$`MISO condition`, "unrecognized event")
    # expect_equal(parsed$`MISO event`, event)
})

test_that("remove_wrong_mRNA removes mRNAs from other chromosomes", {
    event <- read.table(text = "
                      chr6 ALE gene 30620579 30822593  .  +  .
                      chr7 ALE mRNA 30620579 30620982  .  +  .
                      chr7 ALE exon 30620579 30620982  .  +  .
                      chr6 ALE mRNA 30822190 30822593  .  +  .
                      chr6 ALE exon 30822190 30822593  .  +  .")
    new <- remove_wrong_mRNA(event)
    expect_is(new, "data.frame")
    expect_equal(nrow(new), 3)
    expect_equal(new, event[c(1, 4, 5), ])
})

test_that("remove_wrong_mRNA removes mRNAs outside the event boundary", {
    event <- read.table(text = "
                      chr6 ALE gene 30620579 30822593  .  +  .
                      chr6 ALE mRNA 30822190 30822593  .  +  .
                      chr6 ALE exon 30822190 30822593  .  +  .
                      chr6 ALE mRNA 40620579 40620982  .  +  .
                      chr6 ALE exon 40620579 40620982  .  +  .
                      chr6 ALE mRNA 20620579 20620982  .  +  .
                      chr6 ALE exon 20620579 20620982  .  +  .")
    new <- remove_wrong_mRNA(event)
    expect_is(new, "data.frame")
    expect_equal(nrow(new), 3)
    expect_equal(new, event[1:3, ])
})

test_that("list_mRNA creates a list with mRNAs and respective exons", {
    event <- read.table(text = "
                      chr19 TandemUTR gene  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10663759  10664625  .  -  .
                      chr19 TandemUTR exon  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .")
    mRNA <- list_mRNA(event)
    expect_is(mRNA, "list")
    expect_equal(length(mRNA), 4)
    expect_equal(mRNA[[1]], event[2:3, ])
    expect_equal(mRNA[[2]], event[4:5, ])
    expect_equal(mRNA[[3]], event[6:7, ])
    expect_equal(mRNA[[4]], event[8:9, ])
})

test_that("remove_duplicated_mRNA removes duplicated mRNAs", {
    event <- read.table(text = "
                      chr19 TandemUTR gene  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10663759  10664625  .  -  .
                      chr19 TandemUTR exon  10663759  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .
                      chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                      chr19 TandemUTR exon  10664223  10664625  .  -  .")
    mRNA <- list_mRNA(event)
    new <- remove_duplicated_mRNA(mRNA)
    expect_is(new, "list")
    expect_equal(length(new), 2)
    expect_equal(new, mRNA[1:2])
})

## TODO: test parsing multiple events at once