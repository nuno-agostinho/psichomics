context("Diagram alternative splicing event")

test_that("Plot alternative splicing events", {
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
        "SE_19_-_5218431_5218431_5219478_5215606_PTPRS",
        "MXE_15_+_63335142_63335905_63336030_63336226_63336351_63349184_TPM1",
        "MXE_17_-_74090495_74087316_74087224_74086478_74086410_74085401_EXOC7")

    diagrams <- plotSplicingEvent(events)
    expect_is(diagrams, "list")
    expect_is(diagrams, "splicingEventPlotList")
    for (each in diagrams) {
        expect_is(each, "splicingEventPlot")
        expect_is(each, "character")
    }

    # Event identifiers based on a different exon reference
    events <- c(
        "A3SS_15_+_63353138_63353397_63353912_TPM1",
        "A3SS_11_-_61118463_61117894_61117115_CYB561A3",
        "A5SS_21_+_48056459_48055675_48056808_PRMT2",
        "A5SS_1_-_1274667_1274742_1274033_DVL1",
        "AFE_9_+_131901928_131902430_131904724_PPP2R4",
        "AFE_5_-_134688636_134686513_134681747_H2AFY",
        "ALE_12_+_56554104_56555171_56554410_MYL6",
        "ALE_8_-_38314874_38285953_38287466_FGFR1",
        "MXE_15_+_63335142_63336226_63336351_63335905_63336030_63349184_TPM1",
        "MXE_17_-_74090495_74086478_74086410_74087316_74087224_74085401_EXOC7")

    diagrams <- plotSplicingEvent(events)
    expect_is(diagrams, "list")
    expect_is(diagrams, "splicingEventPlotList")
    for (each in diagrams) {
        expect_is(each, "splicingEventPlot")
        expect_is(each, "character")
    }
})
