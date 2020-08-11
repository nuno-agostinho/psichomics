context("Parse MISO splicing events")

library(fastmatch)

test_that("parseMisoAnnotation parses annotation from MISO", {
    folder <- "extdata/eventsAnnotSample/miso_annotation"
    misoOutput <- system.file(folder, package="psichomics")
    
    miso <- parseMisoAnnotation(misoOutput)
    expect_is(miso, "ASevents")
    expect_equal(length(miso), 13)
    expect_equal(unique(miso$Program), "MISO")
    expect_equal(unique(miso$Strand), c("-", "+"))
})

test_that("parseMisoEventID parses an event identifier from MISO", {
    eventID <- c("114785@uc001sok.1@uc001soj.1", "114784@uc001bxm.1@uc001bxn.1")
    # the annotation is one of the GFF3 files needed to run MISO
    gff3 <- system.file("extdata", "miso_AS_annot_example.gff3", 
                        package="psichomics")
    annotation <- read.delim(gff3, header=FALSE, comment.char="#")
    IDcolumn <- 9
    event <- parseMisoEventID(eventID, annotation, IDcolumn)
    
    expect_is(event, "list")
    expect_length(event, 2)
    expect_is(event[[1]], "data.frame")
    expect_equal(unique(event[[1]][[2]]), "AFE")
    expect_equivalent(event[[1]][[3]][[1]], "gene")
    expect_is(event[[2]], "data.frame")
    expect_equal(unique(event[[2]][[2]]), "AFE")
    expect_equivalent(event[[2]][[3]][[1]], "gene")
})

test_that("getValidEvents returns valid events depending on the validator", {
    event <- read.table(text = "
                        chr1 SE gene 17233 18061  .  -  .
                        chr1 SE jfkdjkf 00000 30000  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17526 17742  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE gene 17233 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17606 17742  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE gene 17915 18061  .  -  .
                        chr1 SE gene 17915 18061  .  -  .
                        chr1 SE mRNA 17915 18061  .  -  .
                        ")
    validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
    res <- getValidEvents(event, validator)
    expect_is(res, "data.frame")
    expect_equal(nrow(res), 8)
    expect_equal(as.character(res[[3]]), validator)
    expect_equal(res, event[10:17, ])
})

test_that("getValidEvents returns one valid event", {
    event <- read.table(text = "
                        chr1 SE gene 17233 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17606 17742  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        ")
    validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
    res <- getValidEvents(event, validator)
    expect_is(res, "data.frame")
    expect_equal(nrow(res), 8)
    expect_equal(as.character(res[[3]]), validator)
    expect_equal(res, event)
})

test_that("getValidEvents returns NULL when there are no valid events", {
    event <- read.table(text = "
                        chr1 SE gene 17233 18061  .  -  .
                        chr1 SE dfkdsfjkkfd 00000 30000  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17526 17742  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        chr1 SE gene 17233 18061  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17606 17742  .  -  .
                        chr1 SE ssdfdsfs 110 367  .  -  .
                        chr1 SE mRNA 17233 18061  .  -  .
                        chr1 SE exon 17233 17368  .  -  .
                        chr1 SE exon 17915 18061  .  -  .
                        ")
    validator <- c("gene", "mRNA", rep("exon", 3), "mRNA", rep("exon", 2))
    res <- getValidEvents(event, validator)
    expect_null(res)
})

# annotation <- read.table(text="
#     chr19 AFE gene 50015886 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE mRNA 50028313 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1.A;Parent=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1.A;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE exon 50028313 50028830  .  +  . ID=2217@uc002poi.1@uc002poe.1.A.0;Parent=2217@uc002poi.1@uc002poe.1.A;Name=2217@uc002poi.1@uc002poe.1.A.0;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE mRNA 50015886 50016730  .  +  . ID=2217@uc002poi.1@uc002poe.1.B;Parent=2217@uc002poi.1@uc002poe.1;Name=2217@uc002poi.1@uc002poe.1.B;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE exon 50015886 50016007  .  +  . ID=2217@uc002poi.1@uc002poe.1.B.0;Parent=2217@uc002poi.1@uc002poe.1.B;Name=2217@uc002poi.1@uc002poe.1.B.0;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE exon 50016644 50016730  .  +  . ID=2217@uc002poi.1@uc002poe.1.B.1;Parent=2217@uc002poi.1@uc002poe.1.B;Name=2217@uc002poi.1@uc002poe.1.B.1;gid=2217@uc002poi.1@uc002poe.1
#     chr19 AFE gene 50016492 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1;gid=2217@uc002poh.1@uc002pog.1
#     chr19 AFE mRNA 50017468 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1.A;Parent=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1.A;gid=2217@uc002poh.1@uc002pog.1
#     chr19 AFE exon 50017468 50017743  .  +  . ID=2217@uc002poh.1@uc002pog.1.A.0;Parent=2217@uc002poh.1@uc002pog.1.A;Name=2217@uc002poh.1@uc002pog.1.A.0;gid=2217@uc002poh.1@uc002pog.1
#     chr19 AFE mRNA 50016492 50016730  .  +  . ID=2217@uc002poh.1@uc002pog.1.B;Parent=2217@uc002poh.1@uc002pog.1;Name=2217@uc002poh.1@uc002pog.1.B;gid=2217@uc002poh.1@uc002pog.1
#     chr19 AFE exon 50016492 50016730  .  +  . ID=2217@uc002poh.1@uc002pog.1.B.0;Parent=2217@uc002poh.1@uc002pog.1.B;Name=2217@uc002poh.1@uc002pog.1.B.0;gid=2217@uc002poh.1@uc002pog.1")
# 
# event <- c(1, NA, 7)
# next_event <- c(6, NA, NA)
# 
# test_that("getDataRows retrieves rows of a data frame", {
#     ret <- getDataRows(1, annotation, event, next_event)
#     expect_is(ret, "data.frame")
#     expect_equal(ret, annotation[1:6, 1:8])
#     
#     ret <- getDataRows(3, annotation, event, next_event)
#     expect_is(ret, "data.frame")
#     expect_equal(ret, annotation[7:11, 1:8])
# })
# 
# test_that("getDataRows returns NA if value is outside data frame", {
#     ret <- getDataRows(2, annotation, event, next_event)
#     expect_true(is.na(ret))
# })
# 
# test_that("parseMisoID matches events ID (returns NA if not possible)", {
#     eventID <- c("2217@uc002poi.1@uc002poe.1", "57705@uc009xob.1@uc001jgy.2", 
#                  "2217@uc002poh.1@uc002pog.1")
#     columnID <- 9
#     
#     events <- parseMisoID(eventID, annotation, columnID)
#     expect_equal(length(events), 3)
#     expect_equal(events[[1]], annotation[1:6, 1:8])
#     expect_true(is.na(events[[2]]))
#     expect_equal(events[[3]], annotation[7:11, 1:8])
# })

test_that("parseMisoEvent parses alternative splicing events", {
    event <- read.table(text = "
                        chr1 SE gene 16854 18061 . - .
                        chr1 SE mRNA 16854 18061 . - .
                        chr1 SE exon 16854 17055 . - .
                        chr1 SE exon 17233 17742 . - .
                        chr1 SE exon 17915 18061 . - .
                        chr1 SE mRNA 16854 18061 . - .
                        chr1 SE exon 16854 17955 . - .
                        chr1 SE exon 17915 18061 . - .")
    parsed <- parseMisoEvent(event)
    parsed2 <- parseMisoSE(event)
    expect_equal(parsed, parsed2)
})

test_that("parseMisoSE parses exon skipping junctions", {
    event <- read.table(text = "
                        chr1 SE gene 1370903 1378262  .  +  .
                        chr1 SE mRNA 1370903 1378262  .  +  .
                        chr1 SE exon 1370903 1371201  .  +  .
                        chr1 SE exon 1372702 1372864  .  +  .
                        chr1 SE exon 1374461 1378262  .  +  .
                        chr1 SE mRNA 1370903 1378262  .  +  .
                        chr1 SE exon 1370903 1371201  .  +  .
                        chr1 SE exon 1374461 1378262  .  +  .
                        chr1 SE gene 16854 18061 . - .
                        chr1 SE mRNA 16854 18061 . - .
                        chr1 SE exon 16854 17055 . - .
                        chr1 SE exon 17233 17742 . - .
                        chr1 SE exon 17915 18061 . - .
                        chr1 SE mRNA 16854 18061 . - .
                        chr1 SE exon 16854 17955 . - .
                        chr1 SE exon 17915 18061 . - .")
    parsed <- parseMisoSE(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, rep("chr1", 2))
    expect_equal(parsed$Event.type, rep("SE", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$C1.start[[1]], 1370903)
    expect_equal(parsed$C1.end[[1]],   1371201)
    expect_equal(parsed$A1.start[[1]], 1372702)
    expect_equal(parsed$A1.end[[1]],   1372864)
    expect_equal(parsed$C2.start[[1]], 1374461)
    expect_equal(parsed$C2.end[[1]],   1378262)
    # Minus strand
    expect_equal(parsed$C1.start[[2]], 18061)
    expect_equal(parsed$C1.end[[2]],   17915)
    expect_equal(parsed$A1.start[[2]], 17742)
    expect_equal(parsed$A1.end[[2]],   17233)
    expect_equal(parsed$C2.start[[2]], 17055)
    expect_equal(parsed$C2.end[[2]],   16854)
})

test_that("parseMisoSE doesn't parse unrecognized event", {
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
    parsed <- parseMisoSE(event)
    expect_null(parsed)
})

test_that("parseMisoMXE parses mutually exc. exons junctions", {
    event <- read.table(text = "
                        chr1 MXE gene 764383 788090 . + .
                        chr1 MXE mRNA 764383 788090 . + .
                        chr1 MXE exon 764383 764484 . + .
                        chr1 MXE exon 776580 776753 . + .
                        chr1 MXE exon 787307 788090 . + .
                        chr1 MXE mRNA 764383 788090 . + .
                        chr1 MXE exon 764383 764484 . + .
                        chr1 MXE exon 783034 783186 . + .
                        chr1 MXE exon 787307 788090 . + .
                        chr1 MXE gene 1027371 1051736  .  -  .
                        chr1 MXE mRNA 1027371 1051736  .  -  .
                        chr1 MXE exon 1027371 1027483  .  -  .
                        chr1 MXE exon 1041336 1041429  .  -  .
                        chr1 MXE exon 1051440 1051736  .  -  .
                        chr1 MXE mRNA 1027371 1051736  .  -  .
                        chr1 MXE exon 1027371 1027483  .  -  .
                        chr1 MXE exon 1050402 1050455  .  -  .
                        chr1 MXE exon 1051440 1051736  .  -  .")
    parsed <- parseMisoMXE(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, rep("chr1", 2))
    expect_equal(parsed$Event.type, rep("MXE", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$C1.start[[1]], 764383)
    expect_equal(parsed$C1.end[[1]],   764484)
    expect_equal(parsed$A1.start[[1]], 776580)
    expect_equal(parsed$A1.end[[1]],   776753)
    expect_equal(parsed$A2.start[[1]], 783034)
    expect_equal(parsed$A2.end[[1]],   783186)
    expect_equal(parsed$C2.start[[1]], 787307)
    expect_equal(parsed$C2.end[[1]],   788090)
    # Minus strand
    expect_equal(parsed$C1.start[[2]], 1051736)
    expect_equal(parsed$C1.end[[2]],   1051440)
    expect_equal(parsed$A1.start[[2]], 1041429)
    expect_equal(parsed$A1.end[[2]],   1041336)
    expect_equal(parsed$A2.start[[2]], 1050455)
    expect_equal(parsed$A2.end[[2]],   1050402)
    expect_equal(parsed$C2.start[[2]], 1027483)
    expect_equal(parsed$C2.end[[2]],   1027371)
})

test_that("parseMisoMXE doesn't parse unrecognized event", {
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
    parsed <- parseMisoMXE(event)
    expect_null(parsed)
})

test_that("parseMisoRI parses retained intron junctions", {
    event <- read.table(text = "
                        chr1 RI gene 1223053 1223417  .  +  .
                        chr1 RI mRNA 1222888 1223216  .  +  .
                        chr1 RI exon 1222888 1223216  .  +  .
                        chr1 RI mRNA 1222888 1223216  .  +  .
                        chr1 RI exon 1222888 1222976  .  +  .
                        chr1 RI exon 1223053 1223216  .  +  .
                        chr1 RI gene 17233 17742 . - .
                        chr1 RI mRNA 17233 17742 . - .
                        chr1 RI exon 17233 17742 . - .
                        chr1 RI mRNA 17233 17742 . - .
                        chr1 RI exon 17233 17364 . - .
                        chr1 RI exon 17601 17742 . - .")
    parsed <- parseMisoRI(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, rep("chr1", 2))
    expect_equal(parsed$Event.type, rep("RI", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$C1.start[[1]], 1222888)
    expect_equal(parsed$C1.end[[1]],   1222976)
    expect_equal(parsed$C2.start[[1]], 1223053)
    expect_equal(parsed$C2.end[[1]],   1223216)
    # Minus strand
    expect_equal(parsed$C1.start[[2]], 17742)
    expect_equal(parsed$C1.end[[2]],   17601)
    expect_equal(parsed$C2.start[[2]], 17364)
    expect_equal(parsed$C2.end[[2]],   17233)
})

test_that("parseMisoRI doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr1 RI gene 17233 17742 . - .
                        chr1 RI mRNA 17233 17742 . - .
                        chr1 RI exon 17233 17742 . - .
                        chr10 RI gene 17233 17742 . - .
                        chr10 RI mRNA 17233 17364 . - .
                        chr10 RI exon 17601 17742 . - .")
    parsed <- parseMisoRI(event)
    expect_null(parsed)
})

test_that("parseMisoA5SS parses alt. 5' SS junctions", {
    event <- read.table(text = "
                        chr1 A5SS gene 874655 876686  .  +  .
                        chr1 A5SS mRNA 874655 876686  .  +  .
                        chr1 A5SS exon 874655 874840  .  +  .
                        chr1 A5SS exon 876524 876686  .  +  .
                        chr1 A5SS mRNA 874655 876686  .  +  .
                        chr1 A5SS exon 874655 874792  .  +  .
                        chr1 A5SS exon 876524 876686  .  +  .
                        chr1 A5SS gene 17233 17742 . - .
                        chr1 A5SS mRNA 17233 17742 . - .
                        chr1 A5SS exon 17233 17368 . - .
                        chr1 A5SS exon 17526 17742 . - .
                        chr1 A5SS mRNA 17233 17742 . - .
                        chr1 A5SS exon 17233 17368 . - .
                        chr1 A5SS exon 17606 17742 . - .")
    parsed <- parseMisoA5SS(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, rep("chr1", 2))
    expect_equal(parsed$Event.type, rep("A5SS", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$A2.start[[1]], 874655)
    expect_equal(parsed$A2.end[[1]],   874792)
    expect_equal(parsed$A1.end[[1]],   874840)
    expect_equal(parsed$C2.start[[1]], 876524)
    expect_equal(parsed$C2.end[[1]],   876686)
    # Minus strand
    expect_equal(parsed$A2.start[[2]], 17742)
    expect_equal(parsed$A2.end[[2]],   17606)
    expect_equal(parsed$A1.end[[2]],   17526)
    expect_equal(parsed$C2.start[[2]], 17368)
    expect_equal(parsed$C2.end[[2]],   17233)
})

test_that("parseMisoA5SS doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr1 A5SS gene 17233 17742 . - .
                        chr1 A5SS mRNA 17233 17742 . - .
                        chr1 A5SS exon 17233 17368 . - .
                        chr1 A5SS mRNA 17526 17742 . - .
                        chr1 A5SS mRNA 17233 17742 . - .")
    parsed <- parseMisoA5SS(event)
    expect_null(parsed)
})

test_that("parseMisoA3SS parses alt. 3' SS junctions", {
    event <- read.table(text = "
                        chr1 A3SS gene 898084 898633  .  +  .
                        chr1 A3SS mRNA 898084 898633  .  +  .
                        chr1 A3SS exon 898084 898297  .  +  .
                        chr1 A3SS exon 898412 898633  .  +  .
                        chr1 A3SS mRNA 898084 898633  .  +  .
                        chr1 A3SS exon 898084 898297  .  +  .
                        chr1 A3SS exon 898489 898633  .  +  .
                        chr1 A3SS gene 15796 16765 . - .
                        chr1 A3SS mRNA 15796 16765 . - .
                        chr1 A3SS exon 15796 15947 . - .
                        chr1 A3SS exon 16607 16765 . - .
                        chr1 A3SS mRNA 15796 16765 . - .
                        chr1 A3SS exon 15796 15942 . - .
                        chr1 A3SS exon 16607 16765 . - .")
    parsed <- parseMisoA3SS(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, rep("chr1", 2))
    expect_equal(parsed$Event.type, rep("A3SS", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$C1.start[[1]], 898084)
    expect_equal(parsed$C1.end[[1]],   898297)
    expect_equal(parsed$A1.start[[1]], 898412)
    expect_equal(parsed$A2.start[[1]], 898489)
    expect_equal(parsed$A2.end[[1]],   898633)
    # Minus strand
    expect_equal(parsed$C1.start[[2]], 16765)
    expect_equal(parsed$C1.end[[2]],   16607)
    expect_equal(parsed$A1.start[[2]], 15947)
    expect_equal(parsed$A2.start[[2]], 15942)
    expect_equal(parsed$A2.end[[2]],   15796)
})

test_that("parseMisoA3SS doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr1 A3SS gene 17233 17742 . - .
                        chr1 A3SS mRNA 17233 17742 . - .
                        chr1 A3SS exon 17233 17368 . - .
                        chr1 A3SS mRNA 17526 17742 . - .
                        chr1 A3SS mRNA 17233 17742 . - .")
    parsed <- parseMisoA3SS(event)
    expect_null(parsed)
})

test_that("parseMisoTandemUTR parses tandem UTR junctions", {
    event <- read.table(text = "
                        chr12 TandemUTR gene  28702106  28703099  .  +  .
                        chr12 TandemUTR mRNA  28702106  28703099  .  +  .
                        chr12 TandemUTR exon  28702106  28703099  .  +  .
                        chr12 TandemUTR mRNA  28702106  28702181  .  +  .
                        chr12 TandemUTR exon  28702106  28702181  .  +  .
                        chr19 TandemUTR gene  10663759  10664625  .  -  .
                        chr19 TandemUTR mRNA  10663759  10664625  .  -  .
                        chr19 TandemUTR exon  10663759  10664625  .  -  .
                        chr19 TandemUTR mRNA  10664223  10664625  .  -  .
                        chr19 TandemUTR exon  10664223  10664625  .  -  .")
    parsed <- parseMisoTandemUTR(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, c("chr12", "chr19"))
    expect_equal(parsed$Event.type, rep("TandemUTR", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$A2.start[[1]], 28702106)
    expect_equal(parsed$A2.end[[1]], 28703099)
    expect_equal(parsed$A1.end[[1]], 28702181)
    # Minus strand
    expect_equal(parsed$A2.start[[2]], 10664625)
    expect_equal(parsed$A2.end[[2]], 10663759)
    expect_equal(parsed$A1.end[[2]], 10664223)
})

test_that("parseMisoTandemUTR doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr1 TandemUTR gene 17233 17742 . - .
                        chr1 TandemUTR mRNA 17233 17742 . - .
                        chr1 TandemUTR exon 17233 17368 . - .
                        chr1 TandemUTR mRNA 17526 17742 . - .")
    parsed <- parseMisoTandemUTR(event)
    expect_null(parsed)
})

test_that("parseMisoAFE parses alt. first exon junctions", {
    # Last exons of each mRNA group
    event <- read.table(text = "
                        chr17 AFE gene   4843303   4848515  .  +  .
                        chr17 AFE mRNA   4846721   4848515  .  +  .
                        chr17 AFE exon   4846721   4846814  .  +  .
                        chr17 AFE exon   4847853   4848515  .  +  .
                        chr17 AFE mRNA   4843303   4844288  .  +  .
                        chr17 AFE exon   4843303   4843525  .  +  .
                        chr17 AFE exon   4843782   4844021  .  +  .
                        chr17 AFE exon   4844168   4844288  .  +  .
                        chr6 AFE gene 38561740 38607924  .  -  .
                        chr6 AFE mRNA 38561740 38563843  .  -  .
                        chr6 AFE exon 38561740 38562103  .  -  .
                        chr6 AFE exon 38563546 38563843  .  -  .
                        chr6 AFE mRNA 38565686 38607924  .  -  .
                        chr6 AFE exon 38565686 38565897  .  -  .
                        chr6 AFE exon 38607576 38607924  .  -  .")
    parsed <- parseMisoAFE(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, c("chr17", "chr6"))
    expect_equal(parsed$Event.type, rep("AFE", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$A2.start[[1]], 4844168)
    expect_equal(parsed$A2.end[[1]],   4844288)
    expect_equal(parsed$A1.start[[1]], 4847853)
    expect_equal(parsed$A1.end[[1]],   4848515)
    # Minus strand
    expect_equal(parsed$A2.start[[2]], 38565897)
    expect_equal(parsed$A2.end[[2]],   38565686)
    expect_equal(parsed$A1.start[[2]], 38562103)
    expect_equal(parsed$A1.end[[2]],   38561740)
})

test_that("parseMisoAFE parses a single alt. first exon event", {
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
    parsed <- parseMisoAFE(event)
    expect_equal(parsed$Program, "MISO")
    expect_equal(parsed$Chromosome, "chr17")
    expect_equal(parsed$Event.type, "AFE")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$A2.start[[1]], 4844168)
    expect_equal(parsed$A2.end[[1]],   4844288)
    expect_equal(parsed$A1.start[[1]], 4847853)
    expect_equal(parsed$A1.end[[1]],   4848515)
})

test_that("parseMisoAFE doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr6 AFE gene 30620579 30822593  .  -  .
                        chr6 AFE mRNA 30822190 30822593  .  -  .
                        chr6 AFE mRNA 30620579 30620982  .  -  .")
    parsed <- parseMisoAFE(event)
    expect_null(parsed)
    
    event <- read.table(text = "chr6 AFE gene 30620579 30822593  .  -  .")
    parsed <- parseMisoAFE(event)
    expect_null(parsed)
})

test_that("parseMisoALE parses alt. last exon events", {
    # First exons of each mRNA group
    event <- read.table(text = "
                        chr19 ALE gene 7830524 7834491  .  +  .
                        chr19 ALE mRNA 7830524 7834491  .  +  .
                        chr19 ALE exon 7830524 7831693  .  +  .
                        chr19 ALE exon 7832402 7832514  .  +  .
                        chr19 ALE exon 7833724 7834491  .  +  .
                        chr19 ALE mRNA 7833724 7834491  .  +  .
                        chr19 ALE exon 7833724 7833929  .  +  .
                        chr19 ALE exon 7834259 7834491  .  +  .
                        chr17 ALE gene 26931076 26939102  .  -  .
                        chr17 ALE mRNA 26931076 26939102  .  -  .
                        chr17 ALE exon 26931076 26932233  .  -  .
                        chr17 ALE exon 26938786 26938813  .  -  .
                        chr17 ALE exon 26939062 26939102  .  -  .
                        chr17 ALE mRNA 26934982 26938674  .  -  .
                        chr17 ALE exon 26934982 26935547  .  -  .
                        chr17 ALE exon 26938562 26938674  .  -  .")
    parsed <- parseMisoALE(event)
    expect_equal(parsed$Program, rep("MISO", 2))
    expect_equal(parsed$Chromosome, c("chr19", "chr17"))
    expect_equal(parsed$Event.type, rep("ALE", 2))
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$A1.start[[1]], 7830524)
    expect_equal(parsed$A1.end[[1]],   7831693)
    expect_equal(parsed$A2.start[[1]], 7833724)
    expect_equal(parsed$A2.end[[1]],   7833929)
    # Minus strand
    expect_equal(parsed$A1.start[[2]], 26939102)
    expect_equal(parsed$A1.end[[2]],   26939062)
    expect_equal(parsed$A2.start[[2]], 26938674)
    expect_equal(parsed$A2.end[[2]],   26938562)
})

test_that("parseMisoALE parses a single alt. last exon event", {
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
    parsed <- parseMisoALE(event)
    expect_equal(parsed$Program, "MISO")
    expect_equal(parsed$Chromosome, "chr19")
    expect_equal(parsed$Event.type, "ALE")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$A1.start[[1]], 7830524)
    expect_equal(parsed$A1.end[[1]],   7831693)
    expect_equal(parsed$A2.start[[1]], 7833724)
    expect_equal(parsed$A2.end[[1]],   7833929)
})

test_that("parseMisoALE doesn't parse unrecognized events", {
    event <- read.table(text = "
                        chr6 ALE gene 30620579 30822593  .  -  .
                        chr6 ALE mRNA 30822190 30822593  .  -  .
                        chr6 ALE mRNA 30620579 30620982  .  -  .")
    parsed <- parseMisoALE(event)
    expect_null(parsed)
    
    event <- read.table(text = "
                        chr6 ALE gene 30620579 30822593  .  -  .")
    parsed <- parseMisoALE(event)
    expect_null(parsed)
})

# test_that("remove_wrong_mRNA removes mRNAs from other chromosomes", {
#     event <- read.table(text = "
#                         chr6 ALE gene 30620579 30822593  .  +  .
#                         chr7 ALE mRNA 30620579 30620982  .  +  .
#                         chr7 ALE exon 30620579 30620982  .  +  .
#                         chr6 ALE mRNA 30822190 30822593  .  +  .
#                         chr6 ALE exon 30822190 30822593  .  +  .")
#     new <- remove_wrong_mRNA(event)
#     expect_is(new, "data.frame")
#     expect_equal(nrow(new), 3)
#     expect_equal(new, event[c(1, 4, 5), ])
# })
# 
# test_that("remove_wrong_mRNA removes mRNAs outside the event boundary", {
#     event <- read.table(text = "
#                         chr6 ALE gene 30620579 30822593  .  +  .
#                         chr6 ALE mRNA 30822190 30822593  .  +  .
#                         chr6 ALE exon 30822190 30822593  .  +  .
#                         chr6 ALE mRNA 40620579 40620982  .  +  .
#                         chr6 ALE exon 40620579 40620982  .  +  .
#                         chr6 ALE mRNA 20620579 20620982  .  +  .
#                         chr6 ALE exon 20620579 20620982  .  +  .")
#     new <- remove_wrong_mRNA(event)
#     expect_is(new, "data.frame")
#     expect_equal(nrow(new), 3)
#     expect_equal(new, event[1:3, ])
# })
# 
# test_that("list_mRNA creates a list with mRNAs and respective exons", {
#     event <- read.table(text = "
#                         chr19 TandemUTR gene  10663759  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10663759  10664625  .  -  .
#                         chr19 TandemUTR exon  10663759  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .")
#     mRNA <- list_mRNA(event)
#     expect_is(mRNA, "list")
#     expect_equal(length(mRNA), 4)
#     expect_equal(mRNA[[1]], event[2:3, ])
#     expect_equal(mRNA[[2]], event[4:5, ])
#     expect_equal(mRNA[[3]], event[6:7, ])
#     expect_equal(mRNA[[4]], event[8:9, ])
# })
# 
# test_that("remove_duplicated_mRNA removes duplicated mRNAs", {
#     event <- read.table(text = "
#                         chr19 TandemUTR gene  10663759  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10663759  10664625  .  -  .
#                         chr19 TandemUTR exon  10663759  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .
#                         chr19 TandemUTR mRNA  10664223  10664625  .  -  .
#                         chr19 TandemUTR exon  10664223  10664625  .  -  .")
#     mRNA <- list_mRNA(event)
#     new <- remove_duplicated_mRNA(mRNA)
#     expect_is(new, "list")
#     expect_equal(length(new), 2)
#     expect_equal(new, mRNA[1:2])
# })
