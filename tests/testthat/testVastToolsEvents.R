context("Parse VAST-TOOLS splicing events")

test_that("parseMultipleVastToolsEvents parses multiple events", {
    events <- read.table(text = "
        NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2 0 N,N,N,Bn,S@0,0 0 N,N,N,Bn,S@0,0
        SCYL3   HsaEX0056691    chr1:169839396-169839498        103     chr1:169842837,169839396-169839498,169838180+169838269  S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
        LAP3    HsaEX0035325    chr4:17587574-17587861  288     chr4:17586759,17587574-17587861,17590442        S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
    ")
    parsed <- parseMultipleVastToolsEvents(events)
    expect_is(parsed, "list")
    expect_equal(length(parsed), 3)
    expect_equal(parsed[[1]]$Program, "VAST-TOOLS")
    expect_equal(parsed[[1]]$`C1 end`, 41040823)
    expect_equal(parsed[[2]]$`C1 end`, 169842837)
    expect_equal(parsed[[3]]$`C1 end`, 17586759)
})

test_that("parseVastToolsEvent parses an alternative splicing event", {
    event <- c("NFYA", "HsaEX0042823", "chr6:41046768-41046903", "136", 
               "chr6:41040823,41046768-41046903,41051785", "C2", "0", 
               "N,N,N,Bn,S@0,0", "0", "N,N,N,Bn,S@0,0")
    parsed <- parseVastToolsEvent(event)
    expect_equal(parsed$Program, "VAST-TOOLS")
    expect_equal(parsed$`Gene symbol`, "NFYA")
    expect_equal(parsed$`Event ID`, "HsaEX0042823")
    expect_equal(parsed$`Event type`, "SE")
    expect_equal(parsed$`Inclusion level A`, 0)
    expect_equal(parsed$`Inclusion level B`, 0)
    expect_equal(parsed$Chromosome, "chr6")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   41040823)
    expect_equal(parsed$`A1 start`, 41046768)
    expect_equal(parsed$`A1 end`,   41046903)
    expect_equal(parsed$`C2 start`, 41051785)
})

test_that("parseVastToolsEvent parses an event annotation (event type: SE)", {
    event <- c("NFYA", "HsaEX0042823", "chr6:41046768-41046903", "136", 
               "chr6:41040823,41046768-41046903,41051785", "C2")
    parsed <- parseVastToolsEvent(event)
    expect_equal(parsed$Program, "VAST-TOOLS")
    expect_equal(parsed$`Gene symbol`, "NFYA")
    expect_equal(parsed$`Event ID`, "HsaEX0042823")
    expect_equal(parsed$`Event type`, "SE")
    expect_null(parsed$`Inclusion level A`)
    expect_null(parsed$`Inclusion level B`)
    expect_equal(parsed$Chromosome, "chr6")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   41040823)
    expect_equal(parsed$`A1 start`, 41046768)
    expect_equal(parsed$`A1 end`,   41046903)
    expect_equal(parsed$`C2 start`, 41051785)
})

test_that("parseVastToolsJunctions parses S junctions (+ strand)", {
    coord <- "chr1:169768099,169770024-169770112,169771762"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr1")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   169768099)
    expect_equal(parsed$`A1 start`, 169770024)
    expect_equal(parsed$`A1 end`,   169770112)
    expect_equal(parsed$`C2 start`, 169771762)
})

test_that("parseVastToolsJunctions parses S junctions (- strand)", {
    coord <- "chrX:99887482,99885756-99885863,99884983"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chrX")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   99887482)
    expect_equal(parsed$`A1 start`, 99885863)
    expect_equal(parsed$`A1 end`,   99885756)
    expect_equal(parsed$`C2 start`, 99884983)
})

test_that("parseVastToolsJunctions parses C1 junctions (+ strand)", {
    coord <- "chr2:202014558,202025155+202025436-202025665,202028558"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr2")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   202014558)
    expect_equal(parsed$`A1 start`, c(202025155, 202025436))
    expect_equal(parsed$`A1 end`,   202025665)
    expect_equal(parsed$`C2 start`, 202028558)
})

test_that("parseVastToolsJunctions parses C1 junctions (- strand)", {
    coord <- "chr20:49558568,49557642+49557666-49557746,49557470+49557492+49557495"	
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr20")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   49558568)
    expect_equal(parsed$`A1 start`, 49557746)
    expect_equal(parsed$`A1 end`,   c(49557642, 49557666))
    expect_equal(parsed$`C2 start`, c(49557470, 49557492, 49557495))
})

test_that("parseVastToolsJunctions parses C2 junctions (+ strand)", {
    coord <- "chr6:41040823,41046768-41046903,41051785"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr6")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   41040823)
    expect_equal(parsed$`A1 start`, 41046768)
    expect_equal(parsed$`A1 end`,   41046903)
    expect_equal(parsed$`C2 start`, 41051785)
})

test_that("parseVastToolsJunctions parses C2 junctions (- strand)", {
    coord <- "chr7:8257935,8196746-8196844,8183602"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr7")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   8257935)
    expect_equal(parsed$`A1 start`, 8196844)
    expect_equal(parsed$`A1 end`,   8196746)
    expect_equal(parsed$`C2 start`, 8183602)
})

test_that("parseVastToolsJunctions parses C3 junctions (+ strand) ", {
    coord <- "chr2:202048031,202050494-202050847,202057707"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr2")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   202048031)
    expect_equal(parsed$`A1 start`, 202050494)
    expect_equal(parsed$`A1 end`,   202050847)
    expect_equal(parsed$`C2 start`, 202057707)
})

test_that("parseVastToolsJunctions parses C3 junctions (- strand) ", {
    coord <- "chr6:143828374,143825050-143825211+143825222+143825389,143823259"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr6")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   143828374)
    expect_equal(parsed$`A1 start`, c(143825211, 143825222, 143825389))
    expect_equal(parsed$`A1 end`,   143825050)
    expect_equal(parsed$`C2 start`, 143823259)
})

test_that("parseVastToolsJunctions parses MIC junctions (+ strand)", {
    coord <- "chr1:23385660,23392553-23392564,23395032"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr1")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   23385660)
    expect_equal(parsed$`A1 start`, 23392553)
    expect_equal(parsed$`A1 end`,   23392564)
    expect_equal(parsed$`C2 start`, 23395032)
})

test_that("parseVastToolsJunctions parses MIC junctions (- strand)", {
    coord <- "chr12:6658922,6658642-6658647+6658650,6657991"
    parsed <- parseVastToolsJunctions(coord, "SE")
    expect_equal(parsed$Chromosome, "chr12")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   6658922)
    expect_equal(parsed$`A1 start`, c(6658647, 6658650))
    expect_equal(parsed$`A1 end`,   6658642)
    expect_equal(parsed$`C2 start`, 6657991)
})

test_that("parseVastToolsJunctions parses IR-S junctions (+ strand)", {
    coord <- "chr12:125549925-125550263=125558422-125558525:+"
    parsed <- parseVastToolsJunctions(coord, "RI")
    expect_equal(parsed$Chromosome, "chr12")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 start`, 125549925)
    expect_equal(parsed$`C1 end`,   125550263)
    expect_equal(parsed$`C2 start`, 125558422)
    expect_equal(parsed$`C2 end`,   125558525)
})

test_that("parseVastToolsJunctions parses IR-S junctions (- strand)", {
    coord <- "chr19:58864658-58864693=58864294-58864563:-"
    parsed <- parseVastToolsJunctions(coord, "RI")
    expect_equal(parsed$Chromosome, "chr19")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 start`, 58864693)
    expect_equal(parsed$`C1 end`,   58864658)
    expect_equal(parsed$`C2 start`, 58864563)
    expect_equal(parsed$`C2 end`,   58864294)
})

test_that("parseVastToolsJunctions parses IR-C junctions (+ strand)", {
    coord <- "chr12:9016390-9016604=9020438-9020653:+"
    parsed <- parseVastToolsJunctions(coord, "RI")
    expect_equal(parsed$Chromosome, "chr12")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 start`, 9016390)
    expect_equal(parsed$`C1 end`,   9016604)
    expect_equal(parsed$`C2 start`, 9020438)
    expect_equal(parsed$`C2 end`,   9020653)
})

test_that("parseVastToolsJunctions parses IR-C junctions (- strand)", {
    coord <- "chr19:58861736-58862017=58858719-58859006:-"
    parsed <- parseVastToolsJunctions(coord, "RI")
    expect_equal(parsed$Chromosome, "chr19")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 start`, 58862017)
    expect_equal(parsed$`C1 end`,   58861736)
    expect_equal(parsed$`C2 start`, 58859006)
    expect_equal(parsed$`C2 end`,   58858719)
})

test_that("parseVastToolsJunctions parses Alt3 junctions (+ strand)", {
    coord <- "chr19:36276385,36277798+36277315-36277974"
    parsed <- parseVastToolsJunctions(coord, "A3SS")
    expect_equal(parsed$Chromosome, "chr19")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   36276385)
    expect_equal(parsed$`C2 start`, c(36277798, 36277315))
    expect_equal(parsed$`C2 end`,   36277974)
})

test_that("parseVastToolsJunctions parses Alt3 junctions (- strand)", {
    coord <- "chr17:7133604,7133377-7133474+7133456"
    parsed <- parseVastToolsJunctions(coord, "A3SS")
    expect_equal(parsed$Chromosome, "chr17")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 end`,   7133604)
    expect_equal(parsed$`C2 start`, c(7133474, 7133456))
    expect_equal(parsed$`C2 end`,   7133377)
})

test_that("parseVastToolsJunctions parses Alt3 junctions without downstream
          constitutive exon start", {
    coord <- "chr2:74652618,74652767+74652697-"
    parsed <- parseVastToolsJunctions(coord, "A3SS")
    expect_equal(parsed$Chromosome, "chr2")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 end`,   74652618)
    expect_equal(parsed$`C2 start`, c(74652767, 74652697))
    expect_true(is.na(parsed$`C2 end`))
})

test_that("parseVastToolsJunctions parses Alt5 junctions (+ strand)", {
    coord <- "chr2:74650610-74650654+74650658,74650982"
    parsed <- parseVastToolsJunctions(coord, "A5SS")
    expect_equal(parsed$Chromosome, "chr2")
    expect_equal(parsed$Strand, "+")
    expect_equal(parsed$`C1 start`, 74650610)
    expect_equal(parsed$`C1 end`,   c(74650654, 74650658))
    expect_equal(parsed$`C2 end`,   74650982)
})

test_that("parseVastToolsJunctions parses Alt5 junctions (- strand)", {
    coord <- "chr20:49557666+49557642-49557746,49557470"
    parsed <- parseVastToolsJunctions(coord, "A5SS")
    expect_equal(parsed$Chromosome, "chr20")
    expect_equal(parsed$Strand, "-")
    expect_equal(parsed$`C1 start`, 49557746)
    expect_equal(parsed$`C1 end`,   c(49557666, 49557642))
    expect_equal(parsed$`C2 end`,   49557470)
})