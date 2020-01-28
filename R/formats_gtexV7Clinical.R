gtexV7ClinicalInfoFormat <- function() {
    list(
        tablename   = "Clinical data",    # Name of the created table
        filename    = "GTEx_v7_Annotations_SubjectPhenotypesDS.txt",
        description = "Clinical data of GTEx subjects",
        dataType    = "Clinical data",    # General category for the data
        
        # Transpose the data? This is the first step before parsing the information!
        # After transposition, a row of the current data equals a column of the original
        skip        = 1,     # Rows to skip when parsing file
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of the row or column used to check the format
        
        # File string to check
        check = c("SUBJID", "SEX", "AGE", "DTHHRDY"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1,    # Columns to ignore
        ignoreRows  = 1,    # Rows to ignore
        commentChar = NULL, # String to identify comments (these lines will be ignored)
        
        # Other options
        unique = FALSE, # Remove duplicated rows
        
        # Identity of rows and columns
        rows    = "subjects",
        columns = "attributes",
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Replace gender values with their meaning
            gender <- c("1"="Male", "2"="Female")
            value <- as.character(data[ , "SEX"])
            data[ , "SEX"] <- as.factor(gender[value])
            
            # Replace death circumstance (4-point hardy scale) values with their
            # meaning
            dthhrdy <- c("0"="Ventilator Case", "1"="Violent and fast death",
                         "2"="Fast death of natural causes", 
                         "3"="Intermediate death", "4"="Slow death")
            value <- as.character(data[ , "DTHHRDY"])
            data[ , "DTHHRDY"] <- as.factor(dthhrdy[value])
            
            # Correctly name columns
            match <- c("SUBJID"="Subject ID", "SEX"="Sex", "AGE"="Age",
                       "DTHHRDY"="Death Circumstances")
            colnames(data) <- match[colnames(data)]
            
            return(data)
        }
    )
}

attr(gtexV7ClinicalInfoFormat, "loader") <- "formats"