gtexSampleInfoFormat <- function() {
    list(
        tablename   = "Sample metadata",      # Name of the created table
        filename    = "GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
        description = "Metadata for GTEx samples",
        dataType    = "Sample metadata",      # General category for the data
        
        # Transpose the data? This is the first step before parsing the information!
        # After transposition, a row of the current data equals a column of the original
        skip        = 1,     # Rows to skip when parsing file
        transpose   = FALSE,
        
        # Format checker information
        rowCheck    = TRUE,  # Check format using a row (TRUE) or a column (FALSE)
        checkIndex  = 1,     # Index of the row or column used to check the format
        
        # File string to check
        check = c("SAMPID", "SMATSSCR", "SMCENTER", "SMPTHNTS", "SMRIN",
                  "SMTS"),
        
        # Parsing information
        delim       = "\t", # Delimiter used to separate fields
        colNames    = 1,    # Row to use for column names
        rowNames    = 1,    # Column to use for row names
        ignoreCols  = 1,    # Columns to ignore
        ignoreRows  = 1,    # Rows to ignore
        commentChar = NULL, # String to identify comments (these lines will be ignored)
        
        # Other options
        unique = FALSE,   # Remove duplicated rows
        
        # Default columns to show (NULL to show all)
        show = NULL,
        
        process = function(data) {
            # Replace autolysis values with their meaning
            autolysis <- c("0"="None", "1"="Mild", "2"="Moderate", "3"="Severe")
            value <- as.character(data[ , "SMATSSCR"])
            data[ , "SMATSSCR"] <- as.factor(autolysis[value])

            # Correctly name columns
            match <- c("SAMPID"="Sample ID", "SMATSSCR"="Autolysis Score",
                       "SMNABTCH"="Nucleic Acid Isolation Batch ID",
                       "SMNABTCHT"="Type of nucleic acid isolation batch",
                       "SMNABTCHD"="Date of nucleic acid isolation batch",
                       "SMGEBTCH"="Genotype or Expression Batch ID",
                       "SMGEBTCHD"="Date of genotype or expression batch",
                       "SMGEBTCHT"="Type of genotype or expression batch",
                       "SMCENTER"="BSS collection site",
                       "SMPTHNTS"="Pathology Notes", "SMRIN"="RIN Number",
                       "SMTS"="Tissue Type (area of retrieval)",
                       "SMTSD"="Tissue Type (detail)", "SMUBRID"="Uberon ID",
                       "SMTSISCH"="Total Ischemic time",
                       "SMTSPAX"="PAXgene fixative time",
                       "SMTSTPTREF"="Period of sample procurement",
                       "SMAFRZE"="Samples in GTEx Analysis Freeze",
                       "SMGTC"="Genotype GTC file",
                       "SME2MPRT"="End 2 Mapping Rate",
                       "SMCHMPRS"="Chimeric Pairs",
                       "SMNTRART"="Intragenic Rate",
                       "SMNUMGPS"="Number of Gaps", "SMMAPRT"="Mapping Rate",
                       "SMEXNCRT"="Exonic Rate",
                       "SM550NRM"="5' 50-based normalization",
                       "SMGNSDTC"="Genes Detected",
                       "SMUNMPRT"="Unique Rate of Mapped",
                       "SM350NRM"="3' 50-base normalization",
                       "SMRDLGTH"="Read Length",
                       "SMMNCPB"="Mean Coverage Per Base",
                       "SME1MMRT"="End 1 Mismatch Rate",
                       "SMSFLGTH"="Fragment Length StdDev",
                       "SMESTLBS"="Estimated library size",
                       "SMMPPD"="Total mapped reads",
                       "SMNTERRT"="Intergenic Rate",
                       "SMRRNANM"="rRNA reads", "SMRDTTL"="Total reads",
                       "SMVQCFL"="Failed Vendor QC Check",
                       "SMMNCV"="Mean coefficient of variation",
                       "SMTRSCPT"="Transcripts Detected",
                       "SMMPPDPR"="Mapped Pairs",
                       "SMCGLGTH"="Cumulative Gap Length",
                       "SMGAPPCT"="Gap Percentage",
                       "SMUNPDRD"="Unpaired Reads", "SMNTRNRT"="Intronic Rate",
                       "SMMPUNRT"="Mapped Unique Rate of Total",
                       "SMEXPEFF"="Expression Profiling Efficiency",
                       "SMMPPDUN"="Mapped Unique",
                       "SME2MMRT"="End 2 Mismatch Rate",
                       "SME2ANTI"="End 2 Antisense",
                       "SMALTALG"="Alternative Aligments",
                       "SME2SNSE"="End 2 Sense",
                       "SMMFLGTH"="Fragment Length Mean",
                       "SMSPLTRD"="Split Reads", "SME1ANTI"="End 1 Antisense",
                       "SMBSMMRT"="Base Mismatch Rate",
                       "SME1SNSE"="End 1 Sense",
                       "SME1PCTS"="End 1 % Sense", "SMRRNART"="rRNA Rate",
                       "SME1MPRT"="End 1 Mapping Rate",
                       "SMNUM5CD"="Number Covered 5'",
                       "SMDPMPRT"="Duplication Rate of Mapped",
                       "SME2PCTS"="End 2 % Sense")

            colnames(data) <- match[colnames(data)]

            return(data)
        }
    )
}

attr(gtexSampleInfoFormat, "loader") <- "formats"