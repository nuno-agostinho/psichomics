context("Parse VAST-TOOLS splicing events")

test_that("parseVastToolsAnnotation parses annotation from VAST-tools", {
    folder <- "extdata/eventsAnnotSample/VASTDB/Hsa/TEMPLATES"
    vastToolsOutput <- system.file(folder, package="psichomics")
    
    vast <- parseVastToolsAnnotation(vastToolsOutput)
    expect_is(vast, "ASevents")
    expect_equal(length(vast), 14)
    expect_equal(unique(vast$Program), "VAST-TOOLS")
    expect_equal(unique(vast$Strand), c("-", "+"))
})

test_that("parseVastToolsEvent parses exon skipping event", {
    events <- read.table(
        text = "
        NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2 0 N,N,N,Bn,S@0,0 0 N,N,N,Bn,S@0,0
        SCYL3 HsaEX0056691 chr1:169839396-169839498 103 chr1:169842837,169839396-169839498,169838180+169838269  S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
        LAP3 HsaEX0035325 chr4:17587574-17587861  288     chr4:17586759,17587574-17587861,17590442        S       0.00    N,N,N,Bn,S@0,0  0.00    N,N,N,Bn,S@0,0
        ")
    parsed <- parseVastToolsEvent(events)
    expect_is(parsed, "data.frame")
    expect_equal(nrow(parsed), 3)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_equal(parsed$Gene.symbol[[1]], "NFYA")
    expect_equal(parsed$Inclusion.level.A[[1]], 0)
    expect_equal(parsed$Inclusion.level.B[[1]], 0)
    expect_equal(parsed$C1.end[[1]], 41040823)
    expect_equal(parsed$Event.type[[1]], "SE")
    expect_equal(parsed$Strand[[1]], "+")

    expect_equal(parsed$Program[[2]], "VAST-TOOLS")
    expect_equal(parsed$Gene.symbol[[2]], "SCYL3")
    expect_equal(parsed$Inclusion.level.A[[2]], 0)
    expect_equal(parsed$Inclusion.level.B[[2]], 0)
    expect_equal(parsed$C1.end[[2]], 169842837)
    expect_equal(parsed$Event.type[[2]], "SE")
    expect_equal(parsed$Strand[[2]], "-")

    expect_equal(parsed$Program[[3]], "VAST-TOOLS")
    expect_equal(parsed$Gene.symbol[[3]], "LAP3")
    expect_equal(parsed$Inclusion.level.A[[3]], 0)
    expect_equal(parsed$Inclusion.level.B[[3]], 0)
    expect_equal(parsed$C1.end[[3]], 17586759)
    expect_equal(parsed$Event.type[[3]], "SE")
    expect_equal(parsed$Strand[[3]], "+")
})

test_that("parseVastToolsEvent parses exon skipping annotation", {
    event <- read.table(
        text = "
        NFYA HsaEX0042823 chr6:41046768-41046903 136 chr6:41040823,41046768-41046903,41051785 C2
        ")
    parsed <- parseVastToolsEvent(event)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_equal(parsed$Gene.symbol[[1]], "NFYA")
    expect_equal(parsed$Event.ID[[1]], "HsaEX0042823")
    expect_equal(parsed$Event.type[[1]], "SE")
    expect_null(parsed$Inclusion.level.A)
    expect_null(parsed$Inclusion.level.B)
    expect_equal(parsed$Chromosome[[1]], "chr6")
    expect_equal(parsed$Strand[[1]], "+")
    expect_equal(parsed$C1.end[[1]],   41040823)
    expect_equal(parsed$A1.start[[1]], 41046768)
    expect_equal(parsed$A1.end[[1]],   41046903)
    expect_equal(parsed$C2.start[[1]], 41051785)
})

test_that("parseVastToolsEvent parses alt. 3' splice site annotation", {
    events <- read.table(
        text = "
        DPM1   ENSG00000000419-8-11,8-10-1/2  chr20:49557402-49557470       0    chr20:49558568,49557402-49557492+49557470    Alt3
        DPM1   ENSG00000000419-8-11,8-10-2/2  chr20:49557402-49557492      22    chr20:49558568,49557402-49557492+49557470    Alt3
        C1orf112 ENSG00000000460-14-19,14-18-1/2 chr1:169773253-169773381       0 chr1:169772450,169773253+169773216-169773381    Alt3
        ")
    parsed <- parseVastToolsEvent(events)
    expect_is(parsed, "data.frame")
    expect_equal(nrow(parsed), 3)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_equal(parsed$C1.end[[1]],   49558568)
    expect_equal(parsed$A1.start[[1]], 49557492)
    expect_equal(parsed$C2.start[[1]], 49557470)
    expect_equal(parsed$C2.end[[1]],   49557402)
    expect_equal(parsed$Event.type[[1]], "A3SS")
    expect_equal(parsed$Strand[[1]], "-")

    expect_equal(parsed$C1.end[[2]],   49558568)
    expect_equal(parsed$A1.start[[2]], 49557492)
    expect_equal(parsed$C2.start[[2]], 49557470)
    expect_equal(parsed$C2.end[[2]],   49557402)
    expect_equal(parsed$Event.type[[2]], "A3SS")
    expect_equal(parsed$Strand[[2]], "-")

    expect_equal(parsed$C1.end[[3]],   169772450)
    expect_equal(parsed$A1.start[[3]], 169773253)
    expect_equal(parsed$C2.start[[3]], 169773216)
    expect_equal(parsed$C2.end[[3]],   169773381)
    expect_equal(parsed$Event.type[[3]], "A3SS")
    expect_equal(parsed$Strand[[3]], "+")
})

test_that("parseVastToolsEvent parses alt. 3' splice site annotation with information missing", {
              events <- read.table(
                  text = "
                  DPM1 ENSG00000000419-8-11,8-10-1/2 chr20:49557402-49557470 0 chr20:49558568,-49557492+49557470 Alt3
                  DPM1 ENSG00000000419-8-11,8-10-2/2 chr20:49557402-49557492 22 chr20:49558568,-49557492+49557470 Alt3
                  C1orf112 ENSG00000000460-14-19,14-18-1/2 chr1:169773253-169773381 0 chr1:169772450,169773253+169773216- Alt3
                  ")
              parsed <- parseVastToolsEvent(events)
              expect_is(parsed, "data.frame")
              expect_equal(nrow(parsed), 3)
              expect_equal(parsed$Program[[1]], "VAST-TOOLS")
              expect_equal(parsed$C1.end[[1]],   49558568)
              expect_equal(parsed$A1.start[[1]], 49557492)
              expect_equal(parsed$C2.start[[1]], 49557470)
              expect_true(is.na(parsed$C2.end[[1]]))
              expect_equal(parsed$Event.type[[1]], "A3SS")
              expect_equal(parsed$Strand[[1]], "-")

              expect_equal(parsed$C1.end[[2]],   49558568)
              expect_equal(parsed$A1.start[[2]], 49557492)
              expect_equal(parsed$C2.start[[2]], 49557470)
              expect_true(is.na(parsed$C2.end[[2]]))
              expect_equal(parsed$Event.type[[2]], "A3SS")
              expect_equal(parsed$Strand[[2]], "-")

              expect_equal(parsed$C1.end[[3]],   169772450)
              expect_equal(parsed$A1.start[[3]], 169773253)
              expect_equal(parsed$C2.start[[3]], 169773216)
              expect_true(is.na(parsed$C2.end[[3]]))
              expect_equal(parsed$Event.type[[3]], "A3SS")
              expect_equal(parsed$Strand[[3]], "+")
              })

test_that("parseVastToolsEvent parses alt. 5' splice site annotation", {
    events <- read.table(
        text = "
        TSPAN6 ENSG00000000003-2-4,3-4-1/2 chrX:99891605-99891686 0 chrX:99891605+99891188-99891686,99890743 Alt5
        DPM1 ENSG00000000419-9-11,10-11-1/2 chr20:49557666-49557746 0 chr20:49557666+49557642-49557746,49557470 Alt5
        NFYA ENSG00000001167-3-4,4-4-1/2 chr6:41048553-41048633 0 chr6:41048553-41048633+41048636,41051785 Alt5
        ")
    parsed <- parseVastToolsEvent(events)
    expect_is(parsed, "data.frame")
    expect_equal(nrow(parsed), 3)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_equal(parsed$C1.start[[1]], 99891686)
    expect_equal(parsed$C1.end[[1]],   99891188)
    expect_equal(parsed$A1.end[[1]],   99891605)
    expect_equal(parsed$C2.start[[1]], 99890743)
    expect_equal(parsed$Event.type[[1]], "A5SS")
    expect_equal(parsed$Strand[[1]], "-")

    expect_equal(parsed$C1.start[[2]], 49557746)
    expect_equal(parsed$C1.end[[2]],   49557642)
    expect_equal(parsed$A1.end[[2]],   49557666)
    expect_equal(parsed$C2.start[[2]], 49557470)
    expect_equal(parsed$Event.type[[2]], "A5SS")
    expect_equal(parsed$Strand[[2]], "-")

    expect_equal(parsed$C1.start[[3]], 41048553)
    expect_equal(parsed$C1.end[[3]],   41048636)
    expect_equal(parsed$A1.end[[3]],   41048633)
    expect_equal(parsed$C2.start[[3]], 41051785)
    expect_equal(parsed$Event.type[[3]], "A5SS")
    expect_equal(parsed$Strand[[3]], "+")
})

test_that("parseVastToolsEvent parses alt. 5' splice site annotation with information missing", {
    events <- read.table(
        text = "
        TSPAN6 ENSG00000000003-2-4,3-4-1/2 chrX:99891605-99891686 0 chrX:99891605+99891188-,99890743 Alt5
        DPM1 ENSG00000000419-9-11,10-11-1/2 chr20:49557666-49557746 0 chr20:49557666+49557642-,49557470 Alt5
        NFYA ENSG00000001167-3-4,4-4-1/2 chr6:41048553-41048633 0 chr6:-41048633+41048636,41051785 Alt5
        ")
    parsed <- parseVastToolsEvent(events)
    expect_is(parsed, "data.frame")
    expect_equal(nrow(parsed), 3)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_true(is.na(parsed$C1.start[[1]]))
    expect_equal(parsed$C1.end[[1]],   99891188)
    expect_equal(parsed$A1.end[[1]],   99891605)
    expect_equal(parsed$C2.start[[1]], 99890743)
    expect_equal(parsed$Event.type[[1]], "A5SS")
    expect_equal(parsed$Strand[[1]], "-")

    expect_true(is.na(parsed$C1.start[[2]]))
    expect_equal(parsed$C1.end[[2]],   49557642)
    expect_equal(parsed$A1.end[[2]],   49557666)
    expect_equal(parsed$C2.start[[2]], 49557470)
    expect_equal(parsed$Event.type[[2]], "A5SS")
    expect_equal(parsed$Strand[[2]], "-")

    expect_true(is.na(parsed$C1.start[[3]]))
    expect_equal(parsed$C1.end[[3]],   41048636)
    expect_equal(parsed$A1.end[[3]],   41048633)
    expect_equal(parsed$C2.start[[3]], 41051785)
    expect_equal(parsed$Event.type[[3]], "A5SS")
    expect_equal(parsed$Strand[[3]], "+")
})

test_that("parseVastToolsEvent parses intron retention annotation", {
    events <- read.table(
        text = "
        A2M    ENSG00000175899-A2M:NM_000014:9 chr12:9258942-9259086    145 chr12:9259087-9259201=9258832-9258941:-    IR-S    A2M:NM_000014:9
        A2ML1  ENSG00000166535-A2ML1:NM_144670:1 chr12:8975310-8975777    468 chr12:8975150-8975309=8975778-8975961:+    IR-S  A2ML1:NM_144670:1
        A2ML1 ENSG00000166535-A2ML1:NM_144670:10 chr12:8991819-8993964   2146 chr12:8991709-8991818=8993965-8994132:+    IR-S A2ML1:NM_144670:10
        ")
    parsed <- parseVastToolsEvent(events)
    expect_is(parsed, "data.frame")
    expect_equal(nrow(parsed), 3)
    expect_equal(parsed$Program[[1]], "VAST-TOOLS")
    expect_equal(parsed$C1.start[[1]], 9259201)
    expect_equal(parsed$C1.end[[1]],   9259087)
    expect_equal(parsed$C2.start[[1]], 9258941)
    expect_equal(parsed$C2.end[[1]],   9258832)
    expect_equal(parsed$Event.type[[1]], "RI")
    expect_equal(parsed$Strand[[1]], "-")

    expect_equal(parsed$C1.start[[2]], 8975150)
    expect_equal(parsed$C1.end[[2]],   8975309)
    expect_equal(parsed$C2.start[[2]], 8975778)
    expect_equal(parsed$C2.end[[2]],   8975961)
    expect_equal(parsed$Event.type[[2]], "RI")
    expect_equal(parsed$Strand[[2]], "+")

    expect_equal(parsed$C1.start[[3]], 8991709)
    expect_equal(parsed$C1.end[[3]],   8991818)
    expect_equal(parsed$C2.start[[3]], 8993965)
    expect_equal(parsed$C2.end[[3]],   8994132)
    expect_equal(parsed$Event.type[[3]], "RI")
    expect_equal(parsed$Strand[[3]], "+")
})

test_that("parseVastToolsSE parses exon skipping junctions", {
    junctions <- read.table(
        text = "
        169768099 169770024 169770112 169771762
        99887482 99885756 99885863 99884983
        ")
    parsed <- parseVastToolsSE(junctions)
    expect_equal(parsed$Strand[[1]], "+")
    expect_equal(parsed$C1.end[[1]],   169768099)
    expect_equal(parsed$A1.start[[1]], 169770024)
    expect_equal(parsed$A1.end[[1]],   169770112)
    expect_equal(parsed$C2.start[[1]], 169771762)
    expect_equal(parsed$Strand[[2]], "-")
    expect_equal(parsed$C1.end[[2]],   99887482)
    expect_equal(parsed$A1.start[[2]], 99885863)
    expect_equal(parsed$A1.end[[2]],   99885756)
    expect_equal(parsed$C2.start[[2]], 99884983)
})

test_that("parseVastToolsRI parses intron retention junctions", {
    junctions <- read.table(
        text = "
        125549925 125550263 125558422 125558525
        58864658 58864693 58864294 58864563
        ")
    parsed <- parseVastToolsRI(junctions, c("+", "-"))
    expect_equal(parsed$Strand[[1]], "+")
    expect_equal(parsed$C1.start[[1]], 125549925)
    expect_equal(parsed$C1.end[[1]],   125550263)
    expect_equal(parsed$C2.start[[1]], 125558422)
    expect_equal(parsed$C2.end[[1]],   125558525)

    expect_equal(parsed$Strand[[2]], "-")
    expect_equal(parsed$C1.start[[2]], 58864693)
    expect_equal(parsed$C1.end[[2]],   58864658)
    expect_equal(parsed$C2.start[[2]], 58864563)
    expect_equal(parsed$C2.end[[2]],   58864294)
})

test_that("parseVastToolsA3SS parses alt. 3' splice site junctions", {
    junctions <- rbind(
        c(36276385, list(c(36277798, 36277315)), 36277974),
        c(7133604, 7133377, list(c(7133474, 7133456)))
    )
    parsed <- parseVastToolsA3SS(junctions)
    expect_equal(parsed$Strand[[1]], "+")
    expect_equal(parsed$C1.end[[1]],   36276385)
    expect_equal(parsed$A1.start[[1]], 36277798)
    expect_equal(parsed$C2.start[[1]], 36277315)
    expect_equal(parsed$C2.end[[1]],   36277974)

    expect_equal(parsed$Strand[[2]], "-")
    expect_equal(parsed$C1.end[[2]],   7133604)
    expect_equal(parsed$A1.start[[2]], 7133474)
    expect_equal(parsed$C2.start[[2]], 7133456)
    expect_equal(parsed$C2.end[[2]],   7133377)
})

test_that("parseVastToolsA5SS parses alt. 5' splice site junctions", {
    junctions <- rbind(
        c(74650610, list(c(74650654, 74650658)), 74650982),
        c(list(c(49557666, 49557642), 49557746, 49557470))
    )
    parsed <- parseVastToolsA5SS(junctions)
    expect_equal(parsed$Strand, c("+", "-"))
    # Plus strand
    expect_equal(parsed$C1.start[[1]], 74650610)
    expect_equal(parsed$C1.end[[1]],   74650658)
    expect_equal(parsed$A1.end[[1]],   74650654)
    expect_equal(parsed$C2.start[[1]], 74650982)
    # Minus strand
    expect_equal(parsed$C1.start[[2]], 49557746)
    expect_equal(parsed$C1.end[[2]],   49557642)
    expect_equal(parsed$A1.end[[2]],   49557666)
    expect_equal(parsed$C2.start[[2]], 49557470)
})