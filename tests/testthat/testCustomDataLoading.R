context("Functions to load and process user-provided data")

test_that("Load SRA metadata", {
    # SRA file with metadata
    sraMetadata <- "data/sraMetadata.txt"
    info        <- loadFile(sraMetadata)
    info_old    <- prepareSRAmetadata(sraMetadata, output=NULL)
    expect_identical(info, info_old)

    # Check attributes
    expect_is(info, "data.frame")
    expect_identical(attr(info, "description"), "SRA sample metadata")
    expect_identical(attr(info, "dataType"), "Sample metadata")
    expect_identical(attr(info, "tablename"), "Sample metadata")
    expect_identical(attr(info, "rows"), "samples")
    expect_identical(attr(info, "columns"), "attributes")
    expect_identical(colnames(info)[1], "Run")
    expect_identical(info[[1]], rownames(info))

    # Check content
    info <- lapply(info, as.character)
    expect_identical(info$Run, paste0("SRR636861", 2:5))
    expect_identical((unique(info$`Assay Type`)), "RNA-Seq")
    expect_identical(unique(info$Cell_Line), "HT29")
    expect_identical(unique(info$`Center Name`), "GEO")
    expect_identical(unique(info$Consent), "public")
    expect_identical(unique(info$LibraryLayout), "PAIRED")
    expect_identical(unique(info$LibrarySelection), "cDNA")
    expect_identical(unique(info$LibrarySource), "TRANSCRIPTOMIC")
    expect_identical(unique(info$Platform), "ILLUMINA")
    expect_identical(unique(info$passage), "16-19")
    expect_identical(unique(info$Organism), "Homo sapiens")
})

test_that("Load VAST-TOOLS' inclusion levels", {
    file <- "data/vasttools_incLevels.tab"
    data <- loadFile(file)
    expect_is(data, "data.frame")
    expect_is(data, "sticky")
    expect_equal(ncol(data), 8)

    # Event data
    rowData <- attr(data, "rowData")
    expect_is(rowData, "data.frame")
    expect_is(rowData, "eventData")
    cols <- c("id", "source", "gene", "coordinates", "length",
              "full coordinates", "type", "subtype", "chrom", "start", "end",
              "strand",
              "constitutive1", "alternative1", "alternative2", "constitutive2")
    expect_true(all(cols == colnames(rowData)[seq(cols)]))
    expect_true(all(paste0(colnames(data), "-Q") %in% colnames(rowData)))
    expect_true(unique(rowData$source) == "vast-tools")
    expect_identical(unique(rowData$gene), c("Cdc45", "Scml2", "Apoh",
                                             "Narf", "Cav2", "Scmh1"))
    expect_identical(unique(rowData$type), "SE")
    expect_identical(unique(rowData$subtype), "S")
    expect_identical(unique(rowData$strand), c("-", "+"))
    expect_true(all(rowData$length >= 66 & rowData$length <= 184))

    # Attributes
    expect_identical(attr(data, "dataType"), "Inclusion levels")
    expect_identical(attr(data, "tablename"), "Inclusion levels")
    expect_identical(attr(data, "description"),
                     "VAST-TOOLS' PSI values per alternative splicing event")
    expect_identical(attr(data, "rows"), "alternative splicing events")
    expect_identical(attr(data, "columns"), "samples")
})

test_that("Load VAST-TOOLS' gene expression", {
    # cRPKMS only
    file   <- "data/vasttools_cRPKM.tab"
    cRPKMs <- loadFile(file, multiple=TRUE)
    expect_is(cRPKMs, "data.frame")

    expect_equal(ncol(cRPKMs), 8)
    expect_identical(attr(cRPKMs, "dataType"), "Gene expression")
    expect_identical(attr(cRPKMs, "tablename"), "Gene expression (cRPKM)")
    expect_identical(attr(cRPKMs, "description"),
                     "VAST-TOOLS' gene expression (cRPKM)")
    expect_identical(attr(cRPKMs, "rows"), "genes")
    expect_identical(attr(cRPKMs, "columns"), "samples")

    # cRPKMs and gene read counts
    file <- "data/vasttools_cRPKM_AND_COUNTS.tab"
    data <- loadFile(file, multiple=TRUE)
    expect_is(data, "list")

    counts <- data[[1]]
    expect_is(counts, "data.frame")
    expect_equal(ncol(counts), 8)
    expect_identical(attr(counts, "dataType"), "Gene expression")
    expect_identical(attr(counts, "tablename"), "Gene expression (read counts)")
    expect_identical(attr(counts, "description"),
                     "VAST-TOOLS' gene expression (read counts)")
    expect_identical(attr(counts, "rows"), "genes")
    expect_identical(attr(counts, "columns"), "samples")

    expect_true(attr(cRPKMs, "filename") != attr(data[[2]], "filename"))
    attr(cRPKMs, "filename") <- attr(data[[2]], "filename") <- NULL
    expect_identical(cRPKMs, data[[2]])
    expect_identical(rownames(counts), rownames(cRPKMs))
    expect_identical(colnames(counts), colnames(cRPKMs))
})
