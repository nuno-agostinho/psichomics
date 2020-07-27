context("Test sticky objects")

compareAttrs <- function(object1, object2, attrs=NULL) {
    attrs1 <- attributes(object1)
    attrs2 <- attributes(object2)
    if (is.null(attrs)) attrs <- names(attrs1) # compare all attributes
    
    expect_true(all(attrs %in% names(attrs1)))
    expect_true(all(attrs %in% names(attrs2)))
    expect_equal(attrs1[attrs], attrs2[attrs])
}

test_that("Sticky vectors preserve attributes", {
    vec1 <- 1:5
    attr(vec1, "colour") <- "orange"
    attr(vec1, "animal") <- "fox"
    
    vec2 <- vec1
    vec2 <- vec2[1:3]
    expect_null(attributes(vec2))
    
    vec3 <- preserveAttributes(vec1)[1:3]
    expect_is(vec3, "sticky")
    expect_equivalent(as.vector(vec3), 1:3)
    compareAttrs(vec1, vec3)
    
    vec4 <- vec3[2]
    expect_is(vec4, "sticky")
    expect_equivalent(as.vector(vec4), 2)
    compareAttrs(vec4, vec3)
    
    vec5 <- preserveAttributes(vec1)
    vec5 <- vec5 + 10
    compareAttrs(vec5, vec4)
})

test_that("Sticky data frames preserve attributes", {
    df1 <- data.frame(1:5, 5:9, 45:49, 54:58, 62:66)
    attr(df1, "colour") <- "orange"
    attr(df1, "animal") <- "fox"
    
    df2 <- df1
    df2 <- df2[1:3]
    attrsdf2 <- attributes(df2)
    expect_null(attr(df2, "colour"))
    expect_null(attr(df2, "animal"))
    expect_true(all(c("names", "row.names", "class") %in% names(attrsdf2)))
    
    df3 <- preserveAttributes(df1)[1:3]
    expect_is(df3, "sticky")
    compareAttrs(df3, df2, c("names", "row.names"))
    compareAttrs(df3, df1, c("colour", "animal"))
    
    df4 <- preserveAttributes(df1)[ , 1:3]
    expect_is(df4, "sticky")
    compareAttrs(df4, df3)
    expect_equal(df4, df3)
    
    df5 <- preserveAttributes(df1)[1:3, ]
    expect_is(df5, "sticky")
    expect_identical(data.frame(df5), data.frame(df1[1:3, ]))
    compareAttrs(df5, df1[1:3, ], c("names", "row.names", "colour", "animal"))
    
    df6 <- preserveAttributes(df1)[1, ]
    df7 <- preserveAttributes(df1)[1, , drop=FALSE]
    df8 <- preserveAttributes(df1)[1, , drop=TRUE]
    expect_identical(df6, df7)
    compareAttrs(df6, df7)
    compareAttrs(df6, df8, c("names", "animal", "colour"))
    
    df9  <- preserveAttributes(df1)[ , 1]
    df10 <- preserveAttributes(df1)[ , 1, drop=TRUE]
    df11 <- preserveAttributes(df1)[ , 1, drop=FALSE]
    expect_identical(df9, df10)
    compareAttrs(df9, df10)
    compareAttrs(df9, df11, c("colour", "animal"))
    
    # df12 <- preserveAttributes(df1)
    # df13 <- df12 + 10
    # compareAttrs(df12, df13)
})

test_that("Sticky matrices preserve attributes", {
    mat1 <- as.matrix(data.frame(1:5, 5:9, 45:49, 54:58, 62:66))
    attr(mat1, "colour") <- "orange"
    attr(mat1, "animal") <- "fox"
    
    mat2 <- mat1
    mat2 <- mat2[ , 1:3]
    attrsmat2 <- attributes(mat2)
    expect_null(attr(mat2, "colour"))
    expect_null(attr(mat2, "animal"))
    expect_true(all(c("dim", "dimnames") %in% names(attrsmat2)))
    
    mat3 <- preserveAttributes(mat1)[ , 1:3]
    expect_is(mat3, "sticky")
    compareAttrs(mat3, mat2, c("dim", "dimnames"))
    compareAttrs(mat3, mat1, c("colour", "animal"))
    
    mat5 <- preserveAttributes(mat1)[1:3, ]
    expect_is(mat5, "sticky")
    expect_identical(data.frame(mat5), data.frame(mat1[1:3, ]))
    compareAttrs(mat5, mat1[1:3, ], c("dim", "dimnames"))
    compareAttrs(mat5, mat1, c("colour", "animal"))
    
    mat6 <- preserveAttributes(mat1)[1, ]
    mat7 <- preserveAttributes(mat1)[1, , drop=TRUE]
    mat8 <- preserveAttributes(mat1)[1, , drop=FALSE]
    expect_identical(mat6, mat7)
    compareAttrs(mat6, mat7)
    compareAttrs(mat6, mat8, c("animal", "colour"))
    
    mat9  <- preserveAttributes(mat1)[ , 1]
    mat10 <- preserveAttributes(mat1)[ , 1, drop=TRUE]
    mat11 <- preserveAttributes(mat1)[ , 1, drop=FALSE]
    expect_identical(mat9, mat10)
    compareAttrs(mat9, mat10)
    compareAttrs(mat9, mat11, c("colour", "animal"))
})

test_that("Subset rowData and colData based on extracted sticky data", {
    df1 <- data.frame(1:5, 5:9, 45:49, 54:58, 62:66)
    colnames(df1) <- paste0("sample", 1:5)
    rownames(df1) <- paste0("event", 6:10)
    
    colData <- data.frame(
        cols=colnames(df1), 
        type=c("cancer", "normal", "cancer", "cancer", "normal"),
        sex=c("male", "female", "male", "female", "male"))
    rownames(colData) <- colData[[1]]
    attr(df1, "colData") <- colData
    
    rowData <- data.frame(
        rows=rownames(df1),
        type=c("SE", "RI", "RI", "SE", "SE"),
        gene=c("AAA", "BBB", "CCC", "DDD", "EEE")
    )
    rownames(rowData) <- rowData[[1]]
    attr(df1, "rowData") <- rowData
    
    attr(df1, "colour") <- "orange"
    attr(df1, "animal") <- "fox"
    df1 <- preserveAttributes(df1)
    
    rows <- 2:4
    cols <- 3:5
    df2  <- df1[rows, cols]
    expect_identical(data.frame(df2), data.frame(df1[rows, cols]))
    expect_identical(rownames(df2), rownames(df1)[rows])
    expect_identical(rownames(df2), rownames(attr(df2, "rowData")))
    expect_identical(colnames(df2), colnames(df1)[cols])
    expect_identical(colnames(df2), rownames(attr(df2, "colData")))
    
    df3 <- df1[1, 5]
})

test_that("Transposing a sticky object preserves their attributes", {
    df <- data.frame(1:5, 5:9, 45:49, 54:58, 62:66)
    attr(df, "colour") <- "orange"
    attr(df, "animal") <- "fox"
    
    df1 <- t(df)
    attrsdf1 <- attributes(df1)
    expect_null(attr(df1, "colour"))
    expect_null(attr(df1, "animal"))
    expect_equal(attr(df1, "dimnames")[[1]], colnames(df))
    
    df2 <- t(preserveAttributes(df))
    expect_is(df2, "sticky")
    compareAttrs(df2, df, c("colour", "animal"))
    expect_equal(attr(df2, "dimnames")[[1]], colnames(df))
    
    # Transpose rowData and colData
    colnames(df) <- paste0("sample", 1:5)
    rownames(df) <- paste0("event", 6:10)
    
    colData <- data.frame(
        cols=colnames(df), 
        type=c("cancer", "normal", "cancer", "cancer", "normal"),
        sex=c("male", "female", "male", "female", "male"))
    rownames(colData) <- colData[[1]]
    attr(df, "colData") <- colData
    
    rowData <- data.frame(
        rows=rownames(df),
        type=c("SE", "RI", "RI", "SE", "SE"),
        gene=c("AAA", "BBB", "CCC", "DDD", "EEE")
    )
    rownames(rowData) <- rowData[[1]]
    attr(df, "rowData") <- rowData
    
    attr(df, "colour") <- "orange"
    attr(df, "animal") <- "fox"
    
    df3 <- t(df)
    expect_null(attr(df3, "rowData"))
    expect_null(attr(df3, "colData"))
    
    df4 <- t(preserveAttributes(df))
    expect_is(df4, "sticky")
    expect_equal(attr(df4, "rowData"), attr(df, "colData"))
    expect_equal(attr(df4, "colData"), attr(df, "rowData"))
})
