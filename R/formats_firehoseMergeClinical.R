firehoseMergeClinicalFormat <- function() {
    list(
        tablename   = "Clinical data",      # Name of the created table
        filename    = "clin.merged.txt",    # Name of the file
        description = "Clinical data of the patients",
        dataType    = "Clinical data",      # General category for the data
        
        # Transpose the data? This is the first step before parsing the information!
        # After transposition, a row of the current data equals a column of the original
        transpose   = TRUE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of the row or column used to check the format
        
        # File string to check
        check = c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
                   "admin.disease_code", "admin.file_uuid",
                   "admin.month_of_dcc_upload"),
        
        # Parsing information
        delim       = "\t",  # Delimiter used to separate fields
        colNames    = 1,     # Row to use for column names
        rowNames    = "patient.bcr_patient_barcode",  # Column to use for row names
        ignoreCols  = NULL,  # Columns to ignore
        ignoreRows  = 1:2,   # Rows to ignore
        commentChar = NULL,  # String to identify comments (these lines will be ignored)
        
        # Other options
        unique = FALSE,   # Remove duplicated rows
        
        # Default columns to show (NULL to show all)
        show = c("patient.stage_event.pathologic_stage_tumor_stage",
                 "patient.vital_status",
                 "patient.days_to_death", "patient.days_to_last_followup", 
                 "patient.radiation_therapy", "patient.gender", 
                 "patient.clinical_cqcf.histological_type", "patient.race", 
                 "patient.ethnicity", "patient.race_list.race"),
        
        process = function(data) {
            # Modify column name to be more suggestive
            col <- grep("stage.*pathologic_stage", colnames(data))
            colnames(data)[col] <- paste0(colnames(data)[col], "_tumor_stage")
            
            # Transform patient identifiers to upper case
            rownames(data) <- toupper(rownames(data))
            
            # Remove columns only containing missing values
            onlyNA <- colSums(is.na(data)) == nrow(data)
            data <- data[ , !onlyNA]
            return(data)
        }
    )
}

attr(firehoseMergeClinicalFormat, "loader") <- "formats"