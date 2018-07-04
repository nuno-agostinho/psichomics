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

test_that("Plot transcripts", {
    events <- c(
        "A3SS_15_+_63353138_63353912_63353397_TPM1",
        "A3SS_11_-_61118463_61117115_61117894_CYB561A3",
        "A5SS_21_+_48055675_48056459_48056808_PRMT2", 
        "A5SS_1_-_1274742_1274667_1274033_DVL1",
        "AFE_9_+_131902430_131901928_131904724_PPP2R4",
        "AFE_5_-_134686513_134688636_134681747_H2AFY",
        "ALE_12_+_56554104_56554410_56555171_MYL6",
        "ALE_8_-_38314874_38287466_38285953_FGFR1",
        "SE_9_+_6486925_6492303_6492401_6493826_UHRF2",
        "SE_19_-_5218431_5216778_5216731_5215606_PTPRS",
        "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1",
        "MXE_17_-_74090495_74087316_74087224_74086478_74086410_74085401_EXOC7")
    
    for (event in events) {
        info <- queryEnsemblByEvent(event, species="human", assembly="hg19")
        hc   <- plotTranscripts(info, event=event)
        print(hc)
        expect_is(hc, "shiny.tag.list")
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