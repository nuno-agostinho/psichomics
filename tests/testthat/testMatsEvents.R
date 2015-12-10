context("Parse MATS splicing events")

## parseMatsEvent tests also cover the function parseMatsJunctions
test_that("parseMatsEvent parses alt. 3' splice site event (positive strand)", {
  mats_A3SS <- c(ID = "3658", GeneID = "ENSG00000067715", geneSymbol = "SYT1", 
                 chr = "chr12", strand = "+", longExonStart_0base = "79685787",
                 longExonEnd = "79685910", shortES = "79685796", 
                 shortEE = "79685910", flankingES = "79679566", 
                 flankingEE = "79679751", ID.1 = "3658", IJC_SAMPLE_1 = "252", 
                 SJC_SAMPLE_1 = "102", IJC_SAMPLE_2 = "73", SJC_SAMPLE_2 = "16", 
                 IncFormLen = "58", SkipFormLen = "56",
                 PValue = "0.0342916452301", FDR = "0.274333161841",
                 IncLevel1 = "0.705", IncLevel2 = "0.815", 
                 IncLevelDifference = "-0.11")
  parsed <- parseMatsEvent(mats_A3SS, "A3SS")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`,   79679566)
  expect_equal(parsed$`C1 end`,     79679751)
  expect_equal(parsed$`C2 start`, c(79685787, 79685796))
  expect_equal(parsed$`C2 end`,     79685910)
})

test_that("parseMatsEvent parses alt. 3' splice site event (negative strand)", {
  mats_A3SS <- c(ID = "1234", GeneID = "ENSG00000076108", geneSymbol = "BAZ2A", 
                 chr = "chr12", strand = "-", longExonStart_0base = "57000030", 
                 longExonEnd = "57000179", shortES = "57000030", 
                 shortEE = "57000096", flankingES = "57000416", 
                 flankingEE = "57000517", ID.1 = "1234", IJC_SAMPLE_1 = "0", 
                 SJC_SAMPLE_1 = "18", IJC_SAMPLE_2 = "2", SJC_SAMPLE_2 = "0", 
                 IncFormLen = "112", SkipFormLen = "56",
                 PValue = "0.00337999081157", FDR = "0.0540798529851",
                 IncLevel1 = "0", IncLevel2 = "1", IncLevelDifference = "-1")
  parsed <- parseMatsEvent(mats_A3SS, "A3SS")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 57000517)
  expect_equal(parsed$`C1 end`,   57000416)
  expect_equal(parsed$`C2 start`, c(57000179, 57000096))
  expect_equal(parsed$`C2 end`,   57000030)
})

test_that("parseMatsEvent parses alt. 5' splice event (positive strand)", {
  mats_A5SS <- c(ID = "1366", GeneID = "ENSG00000172465", geneSymbol = "TCEAL1", 
                 chr = "chrX", strand = "+", longExonStart_0base = "102884421", 
                 longExonEnd = "102884501", shortES = "102884421", shortEE = 
                   "102884489", flankingES = "102884812", 
                 flankingEE = "102885881", ID.1 = "1366", IJC_SAMPLE_1 = "2", 
                 SJC_SAMPLE_1 = "2", IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "0", 
                 IncFormLen = "61", SkipFormLen = "56", PValue = "0.489673743698", 
                 FDR = "0.498132777769", IncLevel1 = "0.479", IncLevel2 = "1", 
                 IncLevelDifference = "-0.521")
  parsed <- parseMatsEvent(mats_A5SS, "A5SS")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 102884421)
  expect_equal(parsed$`C1 end`, c(102884501, 102884489))
  expect_equal(parsed$`C2 start`, 102884812)
  expect_equal(parsed$`C2 end`,   102885881)
})

test_that("parseMatsEvent parses alt. 5' splice event (negative strand)", {
  mats_A5SS <- c(ID = "4064", GeneID = "ENSG00000077782", geneSymbol = "FGFR1", 
                 chr = "chr8", strand = "-", longExonStart_0base = "38285863", 
                 longExonEnd = "38285953", shortES = "38285869",
                 shortEE = "38285953", flankingES = "38285438", 
                 flankingEE = "38285611", ID.1 = "4064", IJC_SAMPLE_1 = "2",
                 SJC_SAMPLE_1 = "1", IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "3",
                 IncFormLen = "56", SkipFormLen = "56", 
                 PValue = "0.344797594489", FDR = "0.498132777769",
                 IncLevel1 = "0.667", IncLevel2 = "0.25",
                 IncLevelDifference = "0.417")
  parsed <- parseMatsEvent(mats_A5SS, "A5SS")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 38285953)
  expect_equal(parsed$`C1 end`, c(38285863, 38285869))
  expect_equal(parsed$`C2 start`, 38285611)
  expect_equal(parsed$`C2 end`,   38285438)
})

test_that("parseMatsEvent parses skipping exon event (positive strand)", {
  mats_SE <- c(ID = "4626", GeneID = "ENSG00000151422", geneSymbol = "FER", 
               chr = "chr5", strand = "+", exonStart_0base = "108168470", 
               exonEnd = "108168644", upstreamES = "108133824", 
               upstreamEE = "108134090", downstreamES = "108171408", 
               downstreamEE = "108171508", ID.1 = "4626", IJC_SAMPLE_1 = "16",
               SJC_SAMPLE_1 = "0", IJC_SAMPLE_2 = "0", SJC_SAMPLE_2 = "4",
               IncFormLen = "112", SkipFormLen = "56", 
               PValue = "0.000164083368228", FDR = "0.0164083368228",
               IncLevel1 = "1", IncLevel2 = "0", IncLevelDifference = "1")
  parsed <- parseMatsEvent(mats_SE, "SE")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 108133824)
  expect_equal(parsed$`C1 end`,   108134090)
  expect_equal(parsed$`A1 start`, 108168470)
  expect_equal(parsed$`A1 end`,   108168644)
  expect_equal(parsed$`C2 start`, 108171408)
  expect_equal(parsed$`C2 end`,   108171508)
})

test_that("parseMatsEvent parses skipping exon event (negative strand)", {
  mats_SE <- c(ID = "16170", GeneID = "ENSG00000151914", geneSymbol = "DST",
               chr = "chr6", strand = "-", exonStart_0base = "56463273", 
               exonEnd = "56463507", upstreamES = "56462537", 
               upstreamEE = "56462804", downstreamES = "56464866", 
               downstreamEE = "56465019", ID.1 = "16170", IJC_SAMPLE_1 = "53",
               SJC_SAMPLE_1 = "7", IJC_SAMPLE_2 = "64", SJC_SAMPLE_2 = "2",
               IncFormLen = "112", SkipFormLen = "56", 
               PValue = "0.062949258326", FDR = "0.796580100354", 
               IncLevel1 = "0.791", IncLevel2 = "0.941",
               IncLevelDifference = "-0.15")
  parsed <- parseMatsEvent(mats_SE, "SE")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 56465019)
  expect_equal(parsed$`C1 end`,   56464866)
  expect_equal(parsed$`A1 start`, 56463507)
  expect_equal(parsed$`A1 end`,   56463273)
  expect_equal(parsed$`C2 start`, 56462804)
  expect_equal(parsed$`C2 end`,   56462537)
})

test_that("parseMatsEvent parses intron retention event (positive strand)", {
  mats_RI <- c(ID = "2287", GeneID = "ENSG00000011295", geneSymbol = "TTC19",
               chr = "chr17", strand = "+", riExonStart_0base = "15929853",
               riExonEnd = "15932100", upstreamES = "15929853",
               upstreamEE = "15930016", downstreamES = "15930687",
               downstreamEE = "15932100", ID.1 = "2287", IJC_SAMPLE_1 = "0",
               SJC_SAMPLE_1 = "4", IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "0",
               IncFormLen = "112", SkipFormLen = "56",
               PValue = "0.122321609511", FDR = "0.842461398441",
               IncLevel1 = "0", IncLevel2 = "1", IncLevelDifference = "-1")
  parsed <- parseMatsEvent(mats_RI, "RI")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 15929853)
  expect_equal(parsed$`C1 end`,   15930016)
  expect_equal(parsed$`C2 start`, 15930687)
  expect_equal(parsed$`C2 end`,   15932100)
})

test_that("parseMatsEvent parses intron retention event (negative strand)", {
  mats_RI <- c(ID = "1324", GeneID = "ENSG00000121851", geneSymbol = "POLR3GL", 
               chr = "chr1", strand = "-", riExonStart_0base = "145456990", 
               riExonEnd = "145457309", upstreamES = "145456990",
               upstreamEE = "145457104", downstreamES = "145457235",
               downstreamEE = "145457309", ID.1 = "1324", IJC_SAMPLE_1 = "1",
               SJC_SAMPLE_1 = "12", IJC_SAMPLE_2 = "0", SJC_SAMPLE_2 = "5", 
               IncFormLen = "112", SkipFormLen = "56", 
               PValue = "0.588614424316", FDR = "0.920102552754", 
               IncLevel1 = "0.04", IncLevel2 = "0", IncLevelDifference = "0.04")
  parsed <- parseMatsEvent(mats_RI, "RI")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 145457309)
  expect_equal(parsed$`C1 end`,   145457235)
  expect_equal(parsed$`C2 start`, 145457104)
  expect_equal(parsed$`C2 end`,   145456990)
})

test_that("parseMatsEvent parses mutually exc. exons event (positive strand)", {
  mats_MXE <- c(ID = "217", GeneID = "ENSG00000120251", geneSymbol = "GRIA2", 
                chr = "chr4", strand = "+", X1stExonStart_0base = "158282161", 
                X1stExonEnd = "158282276", X2ndExonStart_0base = "158282689", 
                X2ndExonEnd = "158282804", upstreamES = "158281047", 
                upstreamEE = "158281295", downstreamES = "158283950", 
                downstreamEE = "158284199", ID.1 = "217", IJC_SAMPLE_1 = "10",
                SJC_SAMPLE_1 = "0", IJC_SAMPLE_2 = "7", SJC_SAMPLE_2 = "1",
                IncFormLen = "112", SkipFormLen = "112", 
                PValue = "0.321107675238", FDR = "0.559168057922", 
                IncLevel1 = "1", IncLevel2 = "0.875", 
                IncLevelDifference = "0.125")
  parsed <- parseMatsEvent(mats_MXE, "MXE")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 158281047)
  expect_equal(parsed$`C1 end`,   158281295)
  expect_equal(parsed$`A1 start`, 158282161)
  expect_equal(parsed$`A1 end`,   158282276)
  expect_equal(parsed$`A2 start`, 158282689)
  expect_equal(parsed$`A2 end`,   158282804)
  expect_equal(parsed$`C2 start`, 158283950)
  expect_equal(parsed$`C2 end`,   158284199)
})

test_that("parseMatsEvent parses mutually exc. exons event (negative strand)", {
  mats_MXE <- c(ID = "1388", GeneID = "ENSG00000138468", geneSymbol = "SENP7", 
                chr = "chr3", strand = "-", X1stExonStart_0base = "101117704", 
                X1stExonEnd = "101117899", X2ndExonStart_0base = "101136436", 
                X2ndExonEnd = "101136634", upstreamES = "101090851", 
                upstreamEE = "101090970", downstreamES = "101177798", 
                downstreamEE = "101177896", ID.1 = "1388", IJC_SAMPLE_1 = "0", 
                SJC_SAMPLE_1 = "5", IJC_SAMPLE_2 = "1", SJC_SAMPLE_2 = "3", 
                IncFormLen = "112", SkipFormLen = "112", 
                PValue = "0.327510084993", FDR = "0.559168057922", 
                IncLevel1 = "0", IncLevel2 = "0.25",
                IncLevelDifference = "-0.25")
  parsed <- parseMatsEvent(mats_MXE, "MXE")
  expect_is(parsed, "list")
  expect_equal(parsed$`C1 start`, 101177896)
  expect_equal(parsed$`C1 end`,   101177798)
  expect_equal(parsed$`A1 start`, 101117899)
  expect_equal(parsed$`A1 end`,   101117704)
  expect_equal(parsed$`A2 start`, 101136634)
  expect_equal(parsed$`A2 end`,   101136436)
  expect_equal(parsed$`C2 start`, 101090970)
  expect_equal(parsed$`C2 end`,   101090851)
})