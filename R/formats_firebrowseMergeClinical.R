firebrowseMergeClinicalFormat <- function() {
    list(
        tablename   = "Clinical data",      # Name of the created table
        filename    = "clin.merged.txt",    # Name of the file
        description = "Clinical data of the patients",
        dataType    = "Clinical data",      # General category for the data
        
        # Transpose data? Occurs before parsing: after transposing data, be
        # careful with what is a row/column in the following options
        transpose   = TRUE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using row (TRUE) or column (FALSE)
        checkIndex  = 1,     # Index of row/column to check the format
        
        # File string to check
        check = c("admin.batch_number", "admin.bcr", "admin.day_of_dcc_upload",
                  "admin.disease_code", "admin.file_uuid",
                  "admin.month_of_dcc_upload"),
        
        # Parsing information
        delim       = "\t",   # Delimiter used to separate fields
        colNames    = 1,      # Row to use for column names
        rowNames    = "patient.bcr_patient_barcode",  # Column for row names
        ignoreCols  = NULL,   # Columns to ignore
        ignoreRows  = seq(2), # Rows to ignore
        commentChar = NULL,   # String to identify comments (which are ignored)
        
        # Remove duplicated rows
        unique = FALSE,
        
        # Identity of rows and columns
        rows    = "patients",
        columns = "attributes",
        
        # Default columns to show (NULL to show all)
        show = c(
            # breast cancer-specific information
            "patient.breast_carcinoma_estrogen_receptor_status",
            "patient.breast_carcinoma_progesterone_receptor_status",
            "patient.lab_procedure_her2_neu_in_situ_hybrid_outcome_type",
            "patient.lab_proc_her2_neu_immunohistochemistry_receptor_status",
            "patient.breast_carcinoma_surgical_procedure_name",
            "patient.breast_carcinoma_primary_surgical_procedure_name",
            
            # patient information
            "patient.vital_status",
            "patient.age_at_initial_pathologic_diagnosis",
            "patient.gender", 
            "patient.height",
            "patient.race", 
            "patient.ethnicity", 
            "patient.race_list.race",
            "patient.tumor_samples.tumor_sample.country",
            "patient.menopause_status",
            "patient.barretts_esophagus",
            "patient.h_pylori_infection",
            
            # maybe include all columns with history in the name?
            "patient.tobacco_smoking_history",
            "patient.seizure_history",
            "patient.headache_history",
            "patient.asthma_history",
            "patient.eczema_history",
            "patient.reflux_history",
            "patient.family_history_of_cancer",
            "patient.family_history_of_primary_brain_tumor",
            "patient.family_history_of_stomach_cancer",
            
            # tumour information
            "patient.stage_event.pathologic_stage_tumor_stage",
            "admin.disease_code",
            "patient.tumor_tissue_site",
            "patient.tumor_location",
            "patient.histological_type",
            "patient.primary_pathology.histological_type",
            "patient.clinical_cqcf.histological_type", 
            "patient.neoplasm_histologic_grade",
            paste0("patient.stage_event.tnm_categories.",
                   "pathologic_categories.pathologic_", c("m", "n", "t")),
            "patient.new_tumor_events.new_tumor_event_after_initial_treatment",
            "patient.karnofsky_performance_score",
            "patient.radiation_therapy", 
            "patient.residual_tumor",
            "patient.days_to_death", 
            "patient.days_to_last_followup",
            
            "patient.lymph_node_examined_count",
            "patient.number_of_lymphnodes_positive_by_he",
            "patient.primary_lymph_node_presentation_assessment",
            
            "patient.person_neoplasm_cancer_status",
            "patient.primary_therapy_outcome_success",
            "patient.targeted_molecular_therapy",
            "patient.tissue_source_site"),
        
        process = function(data) {
            # Modify column name to be more suggestive
            col <- grep("stage.*pathologic_stage", colnames(data))
            colnames(data)[col] <- paste0(colnames(data)[col], "_tumor_stage")
            
            # Transform subject identifiers to upper case
            rownames(data) <- toupper(rownames(data))
            
            # Remove columns only containing missing values
            onlyNA <- colSums(is.na(data)) == nrow(data)
            data <- data[ , !onlyNA]
            
            # Replace smoking history numeric values with respective description
            smkCol <- grep("patient.tobacco_smoking_history", colnames(data))
            if (length(smkCol) == 1) {
                smk <- data[ , smkCol]
                smk[smk == 1] <- "Lifelong Non-Smoker"
                smk[smk == 2] <- "Current Smoker"
                smk[smk == 3] <- "Reformed Smoker for > 15 years"
                smk[smk == 4] <- "Reformed Smoker for <= 15 years"
                smk[smk == 5] <- "Reformed Smoker, Unspecified Duration"
                smk[smk == 6] <- "Smoker at Diagnosis"
                smk[smk == 7] <- "Smoking history not documented"
                data[ , smkCol] <- smk
            }
            return(data)
        }
    )
}

attr(firebrowseMergeClinicalFormat, "loader") <- "formats"