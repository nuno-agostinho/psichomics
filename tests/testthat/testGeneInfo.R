context("Gene, transcript and protein annotation")

test_that("Query Ensembl API by event", {
    event <- "SE_12_-_7985318_7984360_7984200_7982602_SLC2A14"
    parsed <- parseEvent(event)
    expect_is(parsed, "data.frame")
    expect_equal(parsed$type, "SE")
    expect_equal(parsed$chrom, "12")
    expect_equal(parsed$strand, "-")
    expect_equal(parsed$gene[[1]], "SLC2A14")
    expect_equal(parsed$pos[[1]], c(7982602, 7985318))
    
    info <- queryEnsemblByEvent(event, species="human", assembly="hg19")
    # Check response
    if (!is.null(info)) {
        expect_is(info, "list")
        # Gene information
        expect_equal(info$seq_region_name, parsed$chrom)
        expect_equal(info$display_name, parsed$gene[[1]])
        expect_equal(info$strand, -1)
        expect_equal(info$source, "ensembl_havana")
        expect_equal(info$object_type, "Gene")
        expect_equal(info$logic_name, "ensembl_havana_gene")
        expect_equal(info$version, 7)
        expect_equal(info$species, "human")
        expect_equal(info$start, 7965108)
        expect_equal(info$end, 8043744)
        expect_equal(info$assembly_name, "GRCh37")
        expect_equal(info$id, "ENSG00000173262")
        expect_equal(info$db_type, "core")
        expect_equal(info$biotype, "protein_coding")
        
        expect_is(info$Transcript, "data.frame")
        expect_is(info$Transcript$Exon, "list")
        expect_is(info$Transcript$Translation, "data.frame")
    }
})

test_that("Convert Ensembl protein ID to UniProt ID", {
    protein <- ensemblToUniprot("ENSP00000445929")
    if (!is.null(protein)) {
        expect_is(protein, "character")
        expect_equivalent(protein, "B7ZAC3")
    }
})

test_that("Plot UniProt protein", {
    plot <- tryCatch(plotProtein("B7ZAC3"), error=return)
    if ("error" %in% class(plot))
        skip("Couldn't resolve host name")
    
    expect_is(plot, "highchart")
    expect_equal(plot$x$type, "chart")
    expect_equivalent(plot$x$hc_opts$xAxis[c("min", "max")], c(0, 520))
    expect_equal(plot$x$hc_opts$chart$type, "area")
    expect_equal(plot$x$hc_opts$chart$zoomType, "x")
    expect_length(plot$x$hc_opts$series, 9)
})