context("Parse MISO splicing events")

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

test_that("getsRowData retrieves rows of a data frame", {
  ret <- getRowsData(1, annotation, event, next_event)
  expect_is(ret, "data.frame")
  expect_equal(ret, annotation[1:6, 1:8])
  
  ret <- getRowsData(3, annotation, event, next_event)
  expect_is(ret, "data.frame")
  expect_equal(ret, annotation[7:11, 1:8])
})

test_that("getsRowData returns NA if value is outside data frame", {
  ret <- getRowsData(2, annotation, event, next_event)
  expect_true(is.na(ret))
})

test_that("matchMisoEventID matches events ID (returns NA if not possible)", {
  eventID <- c("2217@uc002poi.1@uc002poe.1", "57705@uc009xob.1@uc001jgy.2", 
               "2217@uc002poh.1@uc002pog.1")
  columnID <- 9
    
  events <- matchMisoEventID(eventID, annotation, columnID)
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

test_that("parseMisoEventSE parses exon skipping junctions", {
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

test_that("parseMisoEventMXE parses mutually exclusive exons junctions", {
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

test_that("parseMisoEventRI parses intron retention junctions", {
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

test_that("parseMisoEventA5SS parses alt. 5' splice site junctions", {
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

test_that("parseMisoEventA3SS parses alt. 3' splice site junctions", {
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

test_that("parseMisoEventALE parses alt. last exon junctions", {
  event <- read.table(text = "
                    chr12 AFE gene 57916659 57920171  .  +  .
                    chr12 AFE mRNA 57919131 57920171  .  +  .
                    chr12 AFE exon 57919131 57920171  .  +  .
                    chr12 AFE mRNA 57916659 57918199  .  +  .
                    chr12 AFE exon 57916659 57916794  .  +  .
                    chr12 AFE exon 57917812 57917875  .  +  .
                    chr12 AFE exon 57918063 57918199  .  +  .")
  parsed <- parseMisoAFE(event, strand = "+", parsed = list())
  expect_equal(parsed$`C1 start`, 57919131)
  expect_equal(parsed$`C1 end`,   57920171)
  expect_equal(parsed$`A1 start`, 57918063)
  expect_equal(parsed$`A1 end`,   57918199)
})

test_that("parseMisoEventALE parses alt. last exon junctions", {
  event <- read.table(text = "
                      chr6 ALE gene 30620579 30822593  .  +  .
                      chr6 ALE mRNA 30822190 30822593  .  +  .
                      chr6 ALE exon 30822190 30822593  .  +  .
                      chr6 ALE mRNA 30620579 30620982  .  +  .
                      chr6 ALE exon 30620579 30620982  .  +  .")
  parsed <- parseMisoALE(event, strand = "+", parsed = list())
  expect_equal(parsed$`A1 start`, 30822190)
  expect_equal(parsed$`A1 end`,   30822593)
  expect_equal(parsed$`C2 start`, 30620579)
  expect_equal(parsed$`C2 end`,   30620982)
})

test_that("parseMisoEventTandemUTR parses tandem UTR junctions", {
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