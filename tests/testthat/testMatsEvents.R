context("Parse MATS splicing events")

test_that("parseMatsAnnotation parses annotation from rMATS", {
    folder <- "extdata/eventsAnnotSample/mats_output/ASEvents"
    matsOutput <- system.file(folder, package="psichomics")
    mats <- parseMatsAnnotation(matsOutput)
    expect_equal(nrow(mats), 83)
    expect_is(mats, "ASevents")
    expect_equal(length(mats), 15)
    expect_equal(unique(mats$Program), "MATS")
    expect_equal(unique(mats$Strand), c("-", "+"))
    
    # Do not parse novel events
    mats <- parseMatsAnnotation(matsOutput, novelEvents=FALSE)
    expect_equal(nrow(mats), 63)
    expect_is(mats, "ASevents")
    expect_equal(length(mats), 15)
    expect_equal(unique(mats$Program), "MATS")
    expect_equal(unique(mats$Strand), c("-", "+"))
})

test_that("parseMatsEvent parses alt. 3' splice site events", {
    event <- read.table(text = "
        3658 ENSG00000067715 SYT1 chr12 + 79685787 79685910 79685796 79685910 79679566 79679751 3658 252 102 73 16 58 56 0.0342916452301 0.274333161841 0.705 0.815 -0.11
        1234 ENSG00000076108 BAZ2A chr12 - 57000030 57000179 57000030 57000096 57000416 57000517 1234 0 18 2 0 112 56 0.00337999081157 0.0540798529851 0 1 -1
    ")
    parsed <- parseMatsEvent(event, "A3SS")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "C1.start"], 79679566)
    expect_equal(parsed[1, "C1.end"],   79679751)
    expect_equal(parsed[1, "A1.start"], 79685787)
    expect_equal(parsed[1, "A2.start"], 79685796)
    expect_equal(parsed[1, "A2.end"],  79685910)
    expect_equal(parsed[1, "P.value"], 0.0342916452301)
    expect_equal(parsed[1, "FDR"],     0.274333161841)
    expect_equal(parsed[1, "Inclusion.level.A"], 0.705)
    expect_equal(parsed[1, "Inclusion.level.B"], 0.815)
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "C1.start"], 57000517)
    expect_equal(parsed[2, "C1.end"],   57000416)
    expect_equal(parsed[2, "A1.start"], 57000179)
    expect_equal(parsed[2, "A2.start"], 57000096)
    expect_equal(parsed[2, "A2.end"],   57000030)
    expect_equal(parsed[2, "P.value"], 0.00337999081157)
    expect_equal(parsed[2, "FDR"],     0.0540798529851)
    expect_equal(parsed[2, "Inclusion.level.A"], 0)
    expect_equal(parsed[2, "Inclusion.level.B"], 1)
})

test_that("parseMatsEvent parses alt. 3' splice site event annotation", {
    event <- read.table(text = "
        1234 ENSG00000076108 BAZ2A chr12 - 57000030 57000179 57000030 57000096 57000416 57000517
    ")
    parsed <- parseMatsEvent(event, "A3SS")
    expect_is(parsed, "data.frame")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$C1.start, 57000517)
    expect_equal(parsed$C1.end,   57000416)
    expect_equal(parsed$A1.start, 57000179)
    expect_equal(parsed$A2.start, 57000096)
    expect_equal(parsed$A2.end,   57000030)
})

test_that("parseMatsEvent parses alt. 5' splice events", {
    event <- read.table(text = "
        1366 ENSG00000172465 TCEAL1 chrX + 102884421 102884501 102884421 102884489 102884812 102885881 1366 2 2 1 0 61 56 0.489673743698 0.498132777769 0.479 1 -0.521
        4064 ENSG00000077782 FGFR1 chr8 - 38285863 38285953 38285869 38285953 38285438 38285611 4064 2 1 1 3 56 56 0.344797594489 0.498132777769 0.667 0.25 0.417
    ")
    parsed <- parseMatsEvent(event, "A5SS")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "A2.start"], 102884421)
    expect_equal(parsed[1, "A2.end"],   102884489)
    expect_equal(parsed[1, "A1.end"],   102884501)
    expect_equal(parsed[1, "C2.start"], 102884812)
    expect_equal(parsed[1, "C2.end"],   102885881)
    expect_equal(parsed[1, "P.value"], 0.489673743698)
    expect_equal(parsed[1, "FDR"],     0.498132777769)
    expect_equal(parsed[1, "Inclusion.level.A"], 0.479)
    expect_equal(parsed[1, "Inclusion.level.B"], 1)
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "A2.start"], 38285953)
    expect_equal(parsed[2, "A2.end"], 38285869)
    expect_equal(parsed[2, "A1.end"], 38285863)
    expect_equal(parsed[2, "C2.start"], 38285611)
    expect_equal(parsed[2, "C2.end"],   38285438)
    expect_equal(parsed[2, "P.value"], 0.344797594489)
    expect_equal(parsed[2, "FDR"],     0.498132777769)
    expect_equal(parsed[2, "Inclusion.level.A"], 0.667)
    expect_equal(parsed[2, "Inclusion.level.B"], 0.25)
})

test_that("parseMatsEvent parses alt. 5' splice event annotation", {
    event <- read.table(text = "
        4064 ENSG00000077782 FGFR1 chr8 - 38285863 38285953 38285869 38285953 38285438 38285611
    ")
    parsed <- parseMatsEvent(event, "A5SS")
    expect_is(parsed, "data.frame")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$A2.start, 38285953)
    expect_equal(parsed$A2.end,   38285869)
    expect_equal(parsed$A1.end,   38285863)
    expect_equal(parsed$C2.start, 38285611)
    expect_equal(parsed$C2.end,   38285438)
    expect_null(parsed$P.value)
    expect_null(parsed$FDR)
    expect_null(parsed$Inclusion.level.A)
    expect_null(parsed$Inclusion.level.B)
})

test_that("parseMatsEvent parses exon skipping events", {
    event <- read.table(text = "
        4626 ENSG00000151422 FER chr5 + 108168470 108168644 108133824 108134090 108171408 108171508 4626 16 0 0 4 112 56 0.000164083368228 0.0164083368228 1 0 1
        16170 ENSG00000151914 DST chr6 - 56463273 56463507 56462537 56462804 56464866 56465019 16170 53 7 64 2 112 56 0.062949258326 0.796580100354 0.791 0.941 -0.15
    ")
    event <- as.data.frame(rbind(event, deparse.level=FALSE))
    parsed <- parseMatsEvent(event, "SE")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "C1.start"], 108133824)
    expect_equal(parsed[1, "C1.end"],   108134090)
    expect_equal(parsed[1, "A1.start"], 108168470)
    expect_equal(parsed[1, "A1.end"],   108168644)
    expect_equal(parsed[1, "C2.start"], 108171408)
    expect_equal(parsed[1, "C2.end"],   108171508)
    expect_equal(parsed[1, "P.value"], 0.000164083368228)
    expect_equal(parsed[1, "FDR"],     0.0164083368228)
    expect_equal(parsed[1, "Inclusion.level.A"], 1)
    expect_equal(parsed[1, "Inclusion.level.B"], 0)
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "C1.start"], 56465019)
    expect_equal(parsed[2, "C1.end"],   56464866)
    expect_equal(parsed[2, "A1.start"], 56463507)
    expect_equal(parsed[2, "A1.end"],   56463273)
    expect_equal(parsed[2, "C2.start"], 56462804)
    expect_equal(parsed[2, "C2.end"],   56462537)
    expect_equal(parsed[2, "P.value"], 0.062949258326)
    expect_equal(parsed[2, "FDR"],     0.796580100354)
    expect_equal(parsed[2, "Inclusion.level.A"], 0.791)
    expect_equal(parsed[2, "Inclusion.level.B"], 0.941)
})

test_that("parseMatsEvent parses exon skipping event annotation", {
    event <- read.table(text = "
        16170 ENSG00000151914 DST chr6 - 56463273 56463507 56462537 56462804 56464866 56465019
    ")
    parsed <- parseMatsEvent(event, "SE")
    expect_is(parsed, "data.frame")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$C1.start, 56465019)
    expect_equal(parsed$C1.end,   56464866)
    expect_equal(parsed$A1.start, 56463507)
    expect_equal(parsed$A1.end,   56463273)
    expect_equal(parsed$C2.start, 56462804)
    expect_equal(parsed$C2.end,   56462537)
    expect_null(parsed$P.value)
    expect_null(parsed$FDR)
    expect_null(parsed$Inclusion.level.A)
    expect_null(parsed$Inclusion.level.B)
})

test_that("parseMatsEvent parses intron retention events", {
    event <- read.table(text = "
        2287 ENSG00000011295 TTC19 chr17 + 15929853 15932100 15929853 15930016 15930687 15932100 2287 0 4 1 0 112 56 0.122321609511 0.842461398441 0 1 -1
        1324 ENSG00000121851 POLR3GL chr1 - 145456990 145457309 145456990 145457104 145457235 145457309 1324 1 12 0 5 112 56 0.588614424316 0.920102552754 0.04 0 0.04
    ")
    parsed <- parseMatsEvent(event, "RI")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "C1.start"], 15929853)
    expect_equal(parsed[1, "C1.end"],   15930016)
    expect_equal(parsed[1, "C2.start"], 15930687)
    expect_equal(parsed[1, "C2.end"],   15932100)
    expect_equal(parsed[1, "P.value"], 0.122321609511)
    expect_equal(parsed[1, "FDR"], 0.842461398441)
    expect_equal(parsed[1, "Inclusion.level.A"], 0)
    expect_equal(parsed[1, "Inclusion.level.B"], 1)
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "C1.start"], 145457309)
    expect_equal(parsed[2, "C1.end"],   145457235)
    expect_equal(parsed[2, "C2.start"], 145457104)
    expect_equal(parsed[2, "C2.end"],   145456990)
    expect_equal(parsed[2, "P.value"], 0.588614424316)
    expect_equal(parsed[2, "FDR"], 0.920102552754)
    expect_equal(parsed[2, "Inclusion.level.A"], 0.04)
    expect_equal(parsed[2, "Inclusion.level.B"], 0)
})

test_that("parseMatsEvent parses intron retention event annotation", {
    event <- read.table(text = "
        1324 ENSG00000121851 POLR3GL chr1 - 145456990 145457309 145456990 145457104 145457235 145457309
    ")
    parsed <- parseMatsEvent(event, "RI")
    expect_is(parsed, "data.frame")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$C1.start, 145457309)
    expect_equal(parsed$C1.end,   145457235)
    expect_equal(parsed$C2.start, 145457104)
    expect_equal(parsed$C2.end,   145456990)
    expect_null(parsed$P.value)
    expect_null(parsed$FDR)
    expect_null(parsed$Inclusion.level.A)
    expect_null(parsed$Inclusion.level.B)
})

test_that("parseMatsEvent parses mutually exc. exons events", {
    event <- read.table(text = "
        217 ENSG00000120251 GRIA2 chr4 + 158282161 158282276 158282689 158282804 158281047 158281295 158283950 158284199 217 10 0 7 1 112 112 0.321107675238 0.559168057922 1 0.875 0.125
        1388 ENSG00000138468 SENP7 chr3 - 101117704 101117899 101136436 101136634 101090851 101090970 101177798 101177896 1388 0 5 1 3 112 112 0.327510084993 0.559168057922 0 0.25 -0.25
    ")
    parsed <- parseMatsEvent(event, "MXE")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "C1.start"], 158281047)
    expect_equal(parsed[1, "C1.end"],   158281295)
    expect_equal(parsed[1, "A1.start"], 158282161)
    expect_equal(parsed[1, "A1.end"],   158282276)
    expect_equal(parsed[1, "A2.start"], 158282689)
    expect_equal(parsed[1, "A2.end"],   158282804)
    expect_equal(parsed[1, "C2.start"], 158283950)
    expect_equal(parsed[1, "C2.end"],   158284199)
    expect_equal(parsed[1, "P.value"], 0.321107675238)
    expect_equal(parsed[1, "FDR"], 0.559168057922)
    expect_equal(parsed[1, "Inclusion.level.A"], 1)
    expect_equal(parsed[1, "Inclusion.level.B"], 0.875)
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "C1.start"], 101177896)
    expect_equal(parsed[2, "C1.end"],   101177798)
    expect_equal(parsed[2, "A1.start"], 101117899)
    expect_equal(parsed[2, "A1.end"],   101117704)
    expect_equal(parsed[2, "A2.start"], 101136634)
    expect_equal(parsed[2, "A2.end"],   101136436)
    expect_equal(parsed[2, "C2.start"], 101090970)
    expect_equal(parsed[2, "C2.end"],   101090851)
    expect_equal(parsed[2, "P.value"], 0.327510084993)
    expect_equal(parsed[2, "FDR"], 0.559168057922)
    expect_equal(parsed[2, "Inclusion.level.A"], 0)
    expect_equal(parsed[2, "Inclusion.level.B"], 0.25)
})

test_that("parseMatsEvent parses mutually exc. exons event annotation", {
    event <- read.table(text = "
        1388 ENSG00000138468 SENP7 chr3 - 101117704 101117899 101136436 101136634 101090851 101090970 101177798 101177896
    ")
    parsed <- parseMatsEvent(event, "MXE")
    expect_is(parsed, "data.frame")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$C1.start, 101177896)
    expect_equal(parsed$C1.end,   101177798)
    expect_equal(parsed$A1.start, 101117899)
    expect_equal(parsed$A1.end,   101117704)
    expect_equal(parsed$A2.start, 101136634)
    expect_equal(parsed$A2.end,   101136436)
    expect_equal(parsed$C2.start, 101090970)
    expect_equal(parsed$C2.end,   101090851)
    expect_null(parsed$P.value)
    expect_null(parsed$FDR)
    expect_null(parsed$Inclusion.level.A)
    expect_null(parsed$Inclusion.level.B)
})

test_that("parseMatsEvent parses alt. first exon event annotation", {
    event <- read.table(text = "
        54 ENSG00000072958 AP1M1 chr19 + 16308723 16308879 16308967 16309119 16314269 16314426
        1 ENSG00000007372 PAX6 chr11 - 31839356 31839483 31822086 31822194 31816177 31816336
    ")
    parsed <- parseMatsEvent(event, "AFE")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "A2.start"], 16308967)
    expect_equal(parsed[1, "A2.end"],   16309119)
    expect_equal(parsed[1, "A1.start"], 16308723)
    expect_equal(parsed[1, "A1.end"],   16308879)
    expect_equal(parsed[1, "C2.start"], 16314269)
    expect_equal(parsed[1, "C2.end"],   16314426)
    expect_null(parsed[1, "P.value"])
    expect_null(parsed[1, "FDR"])
    expect_null(parsed[1, "Inclusion.level.A"])
    expect_null(parsed[1, "Inclusion.level.B"])
    
    # Minus strand
    expect_equal(parsed[2, "Strand"], "-")
    expect_equal(parsed[2, "A2.start"], 31822194)
    expect_equal(parsed[2, "A2.end"],   31822086)
    expect_equal(parsed[2, "A1.start"], 31839483)
    expect_equal(parsed[2, "A1.end"],   31839356)
    expect_equal(parsed[2, "C2.start"], 31816336)
    expect_equal(parsed[2, "C2.end"],   31816177)
    expect_null(parsed[2, "P.value"])
    expect_null(parsed[2, "FDR"])
    expect_null(parsed[2, "Inclusion.level.A"])
    expect_null(parsed[2, "Inclusion.level.B"])
})

test_that("parseMatsEvent parses alt. last exon event annotation", {
    event <- read.table(text = "
		14 ENSG00000153093 ACOXL chr2 + 111858645 111858828 111851063 111851921 111850441 111850543
        32 ENSG00000223745 RP4-717I23.3 chr1 - 93770667 93771027 93776725 93777165 93790191 93790286
    ")
    parsed <- parseMatsEvent(event, "ALE")
    expect_is(parsed, "data.frame")
    
    # Plus strand
    expect_equal(parsed[1, "Strand"], "+")
    expect_equal(parsed[1, "C1.start"], 111850441)
    expect_equal(parsed[1, "C1.end"],   111850543)
    expect_equal(parsed[1, "A1.start"], 111851063)
    expect_equal(parsed[1, "A1.end"],   111851921)
    expect_equal(parsed[1, "A2.start"], 111858645)
    expect_equal(parsed[1, "A2.end"],   111858828)
    expect_null(parsed[1, "P.value"])
    expect_null(parsed[1, "FDR"])
    expect_null(parsed[1, "Inclusion.level.A"])
    expect_null(parsed[1, "Inclusion.level.B"])
    
    # Minus strand
    expect_equal(parsed[2, "C1.start"], 93790286)
    expect_equal(parsed[2, "C1.end"],   93790191)
    expect_equal(parsed[2, "A1.start"], 93777165)
    expect_equal(parsed[2, "A1.end"],   93776725)
    expect_equal(parsed[2, "A2.start"], 93771027)
    expect_equal(parsed[2, "A2.end"],   93770667)
    expect_null(parsed[2, "P.value"])
    expect_null(parsed[2, "FDR"])
    expect_null(parsed[2, "Inclusion.level.A"])
    expect_null(parsed[2, "Inclusion.level.B"])
})

test_that("parseMatsSE parses a exon skipping event's junctions", {
    junctions <- read.table(
        text = "79685787 79685910 79685796 79685910 79679566 79679751")
    parsed <- parseMatsSE(junctions, strand = "+")
    expect_equal(parsed$C1.start, 79685796)
    expect_equal(parsed$C1.end,   79685910)
    expect_equal(parsed$A1.start, 79685787)
    expect_equal(parsed$A1.end,   79685910)
    expect_equal(parsed$C2.start, 79679566)
    expect_equal(parsed$C2.end,   79679751)
})

test_that("parseMatsMXE parses a mutually exclusive exon event's junctions", {
    junctions <- read.table(
        text= "158282161 158282276 158282689 158282804 158281047 158281295 158283950 158284199")
    parsed <- parseMatsMXE(junctions, strand = "+")
    expect_equal(parsed$C1.start, 158281047)
    expect_equal(parsed$C1.end,   158281295)
    expect_equal(parsed$A1.start, 158282161)
    expect_equal(parsed$A1.end,   158282276)
    expect_equal(parsed$A2.start, 158282689)
    expect_equal(parsed$A2.end,   158282804)
    expect_equal(parsed$C2.start, 158283950)
    expect_equal(parsed$C2.end,   158284199)
})

test_that("parseMatsRI parses an intron retention event's junctions", {
    junctions <- read.table(
        text = "15929853 15932100 15929853 15930016 15930687 15932100")
    parsed <- parseMatsRI(junctions, strand = "+")
    expect_equal(parsed$C1.start, 15929853)
    expect_equal(parsed$C1.end,   15930016)
    expect_equal(parsed$C2.start, 15930687)
    expect_equal(parsed$C2.end,   15932100)
})

test_that("parseMatsA3SS parses an alt. 3' splice site event's junctions", {
    junctions <- read.table(
        text = "79685787 79685910 79685796 79685910 79679566 79679751")
    parsed <- parseMatsA3SS(junctions, strand = "+")
    expect_equal(parsed$C1.start, 79679566)
    expect_equal(parsed$C1.end,   79679751)
    expect_equal(parsed$A1.start, 79685787)
    expect_equal(parsed$A2.start, 79685796)
    expect_equal(parsed$A2.end,   79685910)
})

test_that("parseMatsA3SS parses an alt. 5' splice site event's junctions", {
    junctions <- read.table(
        text="102884421 102884501 102884421 102884489 102884812 102885881")
    parsed <- parseMatsA5SS(junctions, strand = "+")
    expect_equal(parsed$A2.start, 102884421)
    expect_equal(parsed$A2.end,   102884489)
    expect_equal(parsed$A1.end,   102884501)
    expect_equal(parsed$C2.start, 102884812)
    expect_equal(parsed$C2.end,   102885881)
})

test_that("parseMatsAFE parses an alt. first exon event's junctions", {
    junctions <- read.table(
        text = "16308723 16308879 16308967 16309119 16314269 16314426")
    parsed <- parseMatsAFE(junctions, strand = "+")
    expect_equal(parsed$A2.start, 16308967)
    expect_equal(parsed$A2.end,   16309119)
    expect_equal(parsed$A1.start, 16308723)
    expect_equal(parsed$A1.end,   16308879)
    expect_equal(parsed$C2.start, 16314269)
    expect_equal(parsed$C2.end,   16314426)
})

test_that("parseMatsALE parses an alt. last exon event's junctions", {
    junctions <- read.table(
        text = "111858645 111858828 111851063 111851921 111850441 111850543")
    parsed <- parseMatsALE(junctions, strand = "+")
    expect_equal(parsed$C1.start, 111850441)
    expect_equal(parsed$C1.end,   111850543)
    expect_equal(parsed$A1.start, 111851063)
    expect_equal(parsed$A1.end,   111851921)
    expect_equal(parsed$A2.start, 111858645)
    expect_equal(parsed$A2.end,   111858828)
})
