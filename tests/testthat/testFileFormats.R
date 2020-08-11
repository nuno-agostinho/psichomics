context("Correctly load files based on their formats")

formats <- loadFileFormats()
test_that("Load generic junction reads", {
    # Parse chromosome numbers, X, Y, M, Z and W
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chr11:102899:115246:+	345	221	921
            chr21:102899:115246:+	643	254	652
            chr4:102899:115246:+	343	524	168
            chrW:102899:115246:-	602	968	653
            chrW:2899:5246:-	342	456	63
            chrZ:102899:115246:-	654	342	625
            chrM:102899:115246:+	362	8	158
            chrY:102899:115246:-	662	242	214"
    table <- loadFile(txt, formats)
    expect_is(table, "data.frame")
    expect_true(all(startsWith(
        rownames(table), 
        paste0("chr", c("11", "21", "4", "W", "W", "Z", "M", "Y")))))
    expect_equal(colnames(table), paste0("SMPL00", 1:3))
    expect_equivalent(rowSums(table), 
                      c(1487, 1549, 1035, 2223, 861, 1621, 528, 1118))
    
    # Parse different format of junction identifiers
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chr11_102899_115246_+	62	431	72
            chr21:102899:115246:+	43	24	68
            chr4 102899 115246 +	43	24	68
            W 102899-115246 -	62	42	63
            chrW:2899:5246 -	62	42	63
            Z:102899:115246:-	62	42	63
            chromosome M, from 102899 to 115246 +	62	42	63
            chrY, 102899-115246, -	62	42	63"
    table <- loadFile(txt, formats)
    expect_identical(rownames(table), c("chr11:102899:115246:+",
                                        "chr21:102899:115246:+",
                                        "chr4:102899:115246:+",
                                        "chrW:102899:115246:-",
                                        "chrW:2899:5246:-",
                                        "chrZ:102899:115246:-",
                                        "chrM:102899:115246:+",
                                        "chrY:102899:115246:-" ))
    
    # Parse with strand
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chr11:102899:115246:+	62	431	72
            chr21:102899:115246:-	74	53	68
            chr4:102899:115246:-	35	24	38
            chrW:102899:115246:+	25	42	43
            chrW:2899:5246:+	23	42	83
            chrZ:102899:115246:-	57	62	43
            chrM:102899:115246:-	63	22	23
            chrY:102899:115246:+	42	42	13"
    table <- loadFile(txt, formats)
    expect_true(all(endsWith(rownames(table),
                             c("+", "-", "-", "+", "+", "-", "-", "+"))))
    
    # Discard random, alt and Un chromosomes
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chr11:102899:115246:-	62	431	72
            chrUn:102899:115246:-	62	431	72
            chr2_JH159137v1_random:216564:119706:-	82	36	44
            chr3_JH159137v1_alt:616564:619706:-	11	32	42
            chr6_JH159137v1_random:2116564:2119706:-	83	70	95
            chrX_JH159137v1_alt:1564:1594:-	234	32	65"
    expect_warning(table <- loadFile(txt, formats))
    expect_equal(nrow(table), 1)
    
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chrUn:102899:115246:-	62	431	72
            chr2_JH159137v1_random:216564:119706:-	82	36	44
            chr3_JH159137v1_alt:616564:619706:-	11	32	42
            chr6_JH159137v1_random:2116564:2119706:-	83	70	95
            chrX_JH159137v1_alt:1564:1594:-	234	32	65"
    expect_warning(table <- loadFile(txt, formats))
    expect_null(table)
    
    # Duplicated junctions are discarded with a warning
    txt <- "Junction ID	SMPL001	SMPL002	SMPL003
            chr11:102899:115246:-	62	431	72
            chr2:102899:115246:-	43	24	68
            chr11:102899:115246:-	43	73	45"
    expect_warning(loadFile(txt, formats))
})
