#' @rdname appUI
#'
#' @importFrom shiny textInput
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsPopover
localDataUI <- function(id, panel) {
    ns <- NS(id)

    sampleInfoBrowser <- fileBrowserInput(
        ns("sampleInfo"), "Sample information",
        placeholder="No file selected", clearable=TRUE,
        info=TRUE, infoFUN=bsPopover, infoTitle="Sample information",
        infoContent=paste(
            tags$ul(
                class="popover-list",
                tags$li("Tab-separated values (TSV)"),
                tags$li(
                    "Sample identifiers (rows) and their attributes (columns)"),
                tags$li("The first column must contain sample identifiers",
                        "and be named", tags$kbd("Sample ID")),
                tags$li("Optionally, indicate the subject associated to",
                        "each sample in a column named",
                        tags$kbd("Subject ID"))),
            tags$hr(), helpText("Example:"), tags$table(
                class="table table-condensed",
                tags$thead(
                    tableRow("Sample ID", "Type", "Tissue", "Subject ID",
                             th=TRUE)),
                tags$tbody(
                    tableRow("SMP-01", "Tumour", "Lung", "SUBJ-03"),
                    tableRow("SMP-02", "Normal", "Blood", "SUBJ-12"),
                    tableRow("SMP-03", "Normal", "Blood", "SUBJ-25")))))
    subjectInfoBrowser <- fileBrowserInput(
        ns("subjectInfo"), "Subject information",
        placeholder="No file selected", clearable=TRUE,
        info=TRUE, infoFUN=bsPopover, infoTitle="Subject information",
        infoContent=paste(
            tags$ul(
                class="popover-list",
                tags$li("Tab-separated values (TSV)"),
                tags$li("Subject identifiers (rows) and their attributes",
                        "(columns)"),
                tags$li("The first column must contain subject identifiers and",
                        "be named", tags$kbd("Subject ID"))),
            tags$hr(),
            helpText("Example:"), tags$table(
                class="table table-condensed",
                tags$thead(
                    tableRow("Subject ID", "Age", "Gender", "Race",
                             th=TRUE)),
                tags$tbody(
                    tableRow("SUBJ-01", "34", "Female", "Black"),
                    tableRow("SUBJ-02", "22", "Male", "Black"),
                    tableRow("SUBJ-03", "58", "Female", "Asian")))))
    junctionQuantBrowser <- fileBrowserInput(
        ns("junctionQuant"), "Exon-exon junction read counts",
        placeholder="No file selected", clearable=TRUE,
        info=TRUE, infoFUN=bsPopover,
        infoTitle="Exon-exon junction read counts",
        infoContent=paste(
            tags$ul(
                class="popover-list",
                tags$li("Tab-separated values (TSV)"),
                tags$li("Read counts of exon-exon junctions (rows) across",
                        "samples (columns)"),
                tags$li(
                    "The first column must contain junction identifiers",
                    "and be named", tags$kbd("Junction ID")),
                tags$li(
                    "Only chromosome number and capital letters X, Y, Z, W",
                    "and M, followed by the genomic regions are supported;",
                    "acceptable junction identifiers include:",
                    tags$kbd("10_18748_21822"), ",",
                    tags$kbd("chromosome 10 (18748 to 21822)"), "and",
                    tags$kbd("chr10:18748-21822")),
                tags$li(
                    "Optionally, indicate the strand with", tags$kbd("+"),
                    "or", tags$kbd("-"),
                    "at the end of the junction identifier; e.g.",
                    tags$kbd("10:3213:9402:+"), "and",
                    tags$kbd("chr10:3213-9402 -")),
                tags$li(
                    "Rows whose junction identifiers contain",
                    tags$kbd("alt"), ",", tags$kbd("random"), "or",
                    tags$kbd("Un"), "in chromosome names are discarded")),
            tags$hr(), helpText("Example:"), tags$table(
                class="table table-condensed",
                tags$thead(
                    tableRow("Junction ID", "SMP-18", "SMP-03", th=TRUE)),
                tags$tbody(
                    tableRow("10:6752-7393", "4", "0"),
                    tableRow("10:18748-21822", "8", "46"),
                    tableRow("10:24257-25325", "83", "65")))))
    addMultipleFiles <- tagList(
        helpText("All fields below are optional."),
        h3(icon("vial"), "Metadata"),
        sampleInfoBrowser,
        subjectInfoBrowser,
        tags$hr(),
        h3(icon("dna"), "Molecular data"),
        geneExprFileInput(ns("geneExpr"), clearable=TRUE),
        junctionQuantBrowser,
        ASquantFileInput(ns("ASquant"), clearable=TRUE),
        tags$hr(),
        textInput(ns("userFilesCategory"), label="Dataset name", width = "100%",
                  value="User dataset", placeholder="Name to identify dataset"),
        processButton(ns("loadMultipleFiles"), "Load files"))

    addFolder <- tagList(
        helpText("For your convenience, move all files to a single folder.",
                 "Only files supported by psichomics will be loaded."),
        fileBrowserInput(ns("localFolder"), "Folder where data is stored",
                         placeholder="No folder selected",
                         value=getDownloadsFolder()),
        textInput(ns("localCategory"), label="Dataset name",
                  placeholder="Name to identify dataset", width = "100%"),
        selectizeInput(ns("localIgnore"), "Files/directories to ignore",
                       choices=getFirebrowseDataTypes(), width = "100%",
                       multiple=TRUE, options=list(
                           # Allow to add new items
                           create=TRUE, createOnBlur=TRUE,
                           placeholder="These files will not be loaded")),
        processButton(ns("acceptFile"), "Load files"))

    panel(style="info", title=list(icon("plus-circle"),
                                   "User-provided data loading"),
          value="Load local files",
          uiOutput(ns("localDataModal")),
          tabsetPanel(
              tabPanel("File by file", addMultipleFiles),
              tabPanel("By folder", addFolder)))
}

#' Prepare user-provided files to be loaded into psichomics
#'
#' @param file Character: path to file
#' @param output Character: path of output file (if \code{NULL}, only returns
#' the data without saving it to a file)
#'
#' @importFrom data.table fread fwrite
#'
#' @return Prepared file (if \code{output != NULL}) and object
#' @export
prepareSRAmetadata <- function(file, output="psichomics_metadata.txt") {
    data <- loadFile(file, SraRunTableSampleInfoFormat())
    if (!is.null(output)) {
        fwrite(data, output, sep="\t")
        return(invisible(data))
    } else {
        return(data)
    }
}

#' Process SRA quantification data
#'
#' @param files Character: path to SRA quantification files
#' @param data Data frame: processed quantification data
#' @param IDcolname Character: name of the column containing the identifiers
#'
#' @importFrom utils askYesNo
#'
#' @return Process file
#' @keywords internal
processSRAdata <- function(files, data, IDcolname) {
    # Add sample names
    samples <- names(files)
    if (is.null(samples)) {
        # Remove STAR filename end
        if (IDcolname == "Gene ID")
            filenameSuffix <- "ReadsPerGene"
        else if (IDcolname == "Junction ID")
            filenameSuffix <- "SJ"
        filenameSuffix <- paste0(filenameSuffix, "\\.out\\.tab$")

        samples <- gsub(filenameSuffix, "", unlist(files))
        samples <- basename(samples)
    }
    colnames(data) <- as.character(samples)

    quant <- cbind(rownames(data), data)
    setnames(quant, "V1", IDcolname)
    return(quant)
}

#' Save processed SRA data in file
#'
#' @param data Object to save
#' @param output Character: output filename (if \code{NULL}, no file is saved)
#'
#' @return If \code{output = NULL}, save input to a file and return it as
#'   invisible; otherwise, just return the input
#' @keywords internal
saveProcessedSRAdata <- function(data, output=NULL) {
    # Save data to given path
    if (!is.null(output)) {
        if (file.exists(output)) {
            msg <- sprintf(
                "File %s already exists. Do you want to overwrite it?", output)
            allowOverwrite <- askYesNo(msg, default=FALSE,
                                       prompts=c("Overwrite", "No", "Cancel"))
            if (!allowOverwrite || is.na(allowOverwrite))
                return(invisible(NULL))
        } else {
            allowOverwrite <- FALSE
        }
        fwrite(data, output, sep="\t", na=0, quote=FALSE)
        message(sprintf("File %s was %s", output,
                        ifelse(allowOverwrite, "overwritten", "created")))
        return(invisible(data))
    } else {
        return(data)
    }
}

#' @rdname prepareSRAmetadata
#'
#' @param ... Character: path of (optionally named) input files (see Examples)
#' @param startOffset Numeric: value to offset start position
#' @param endOffset Numeric: value to offset end position
#'
#' @importFrom data.table fwrite
#' @export
#'
#' @examples
#' \dontrun{
#' prepareJunctionQuant("Control rep1"=junctionFile1,
#'                      "Control rep2"=junctionFile2,
#'                      "KD rep1"=junctionFile3,
#'                      "KD rep2"=junctionFile4)
#' }
prepareJunctionQuant <- function(..., output="psichomics_junctions.txt",
                                 startOffset=NULL, endOffset=NULL) {
    # Detect splice-aware aligner used
    # TODO(NunoA): support TopHat
    files <- list(...)

    # Prepare junction quantification accordingly
    data  <- prepareJunctionQuantSTAR(..., startOffset=startOffset,
                                      endOffset=endOffset)
    quant <- processSRAdata(files, data, "Junction ID")
    quant <- saveProcessedSRAdata(quant, output)
    return(invisible(quant))
}

#' @inherit prepareSRAmetadata
#' @importFrom data.table fread setnames setkeyv setorderv
#' @keywords internal
prepareJunctionQuantSTAR <- function(..., startOffset=-1, endOffset=+1) {
    if (is.null(startOffset)) startOffset <- -1
    if (is.null(endOffset))   endOffset   <- +1

    files <- list(...)
    if (length(files) == 1) files <- unlist(files)
    joint <- NULL
    for (file in files) {
        cat(sprintf("Processing %s...", file), fill=TRUE)
        table    <- fread(file)[, c(1, 2, 3, 4, 7)]
        table$V2 <- table$V2 + startOffset
        table$V3 <- table$V3 + endOffset
        joint    <- c(joint, list(table))
    }

    index <- 0
    lapply(joint, function(table) {
        index <<- index + 1
        setnames(table, "V7", paste0("col", index))
        setkeyv(table, c("V1", "V2", "V3", "V4"))
    })

    # Merge together junction quantification from different samples
    cat("Merging junction quantification files...", fill=TRUE)
    junctionQuant <- Reduce(function(...) merge(..., all=TRUE), joint)
    # setorderv(junctionQuant, cols=c("V1", "V2", "V3"))

    # Use splice junction location as row names
    # TODO (NunoA): what to do in case the strand is 0 (i.e. undefined)? Maybe
    # duplicate entry and append both a positive and a negative sign
    strand <- ifelse(junctionQuant$V4 == "1", "+", "-")
    cat("Preparing event identifiers...", fill=TRUE)
    ns     <- with(junctionQuant, paste(V1, V2, V3, strand, sep=":"))
    junctionQuant <- junctionQuant[ , -c(1, 2, 3, 4)]
    rownames(junctionQuant) <- ns
    return(junctionQuant)
}

#' @rdname prepareSRAmetadata
#'
#' @param strandedness Character: strandedness of RNA-seq protocol; may be one
#' of the following: \code{unstraded}, \code{stranded} or
#' \code{stranded (reverse)}
#'
#' @importFrom data.table fwrite
#' @export
#'
#' @examples
#' \dontrun{
#' prepareGeneQuant("Control rep1"=geneCountFile1,
#'                  "Control rep2"=geneCountFile2,
#'                  "KD rep1"=geneCountFile3,
#'                  "KD rep2"=geneCountFile4)
#' }
prepareGeneQuant <- function(..., output="psichomics_gene_counts.txt",
                             strandedness=c("unstranded", "stranded",
                                            "stranded (reverse)")) {
    strandedness <- match.arg(strandedness)

    # Detect splice-aware aligner used
    # TODO(NunoA): support TopHat
    files <- list(...)

    # Prepare file accordingly
    data  <- prepareGeneQuantSTAR(..., strandedness=strandedness)
    quant <- processSRAdata(files, data, "Gene ID")
    quant <- saveProcessedSRAdata(quant, output)
    return(invisible(quant))
}

#' @rdname prepareJunctionQuantSTAR
#' @importFrom data.table fread setnames setkeyv setorderv
prepareGeneQuantSTAR <- function(..., strandedness=c("unstranded", "stranded",
                                                     "stranded (reverse)")) {
    strandedness <- match.arg(strandedness)
    strandedness <- switch(strandedness,
                           "unstranded"=2, "stranded"=3, "stranded (reverse)"=4)

    files <- list(...)
    if (length(files) == 1) files <- unlist(files)
    joint <- NULL
    for (file in files) {
        cat(sprintf("Processing %s...", file), fill=TRUE)
        table  <- fread(file, skip=4)
        table  <- table[ , c(1, strandedness), with=FALSE]
        joint  <- c(joint, list(table))
    }

    index <- 0
    lapply(joint, function(table) {
        index <<- index + 1
        setnames(table, colnames(table)[[2]], paste0("col", index))
        setkeyv(table, "V1")
    })

    # Merge together files from different samples
    cat("Merging gene read count files...", fill=TRUE)
    geneQuant <- Reduce(function(...) merge(..., all=TRUE), joint)
    # setorderv(geneQuant, cols=c("V1", "V2", "V3"))

    # Use gene names as row names
    cat("Preparing event identifiers...", fill=TRUE)
    ns <- geneQuant$V1
    geneQuant <- geneQuant[ , -1]
    rownames(geneQuant) <- ns
    return(geneQuant)
}

# Remove datasets containing the same information
removeRedundantDatasets <- function(data) {
    toRemove <- NULL
    if (length(data) > 1) {
        for (i in seq(length(data) - 1)) {
            for (j in seq(i + 1, length(data))) {
                isRedundant <- isTRUE(all.equal(data[[i]], data[[j]],
                                                check.attributes=FALSE))
                if (isRedundant) toRemove <- c(toRemove, j)
            }
        }
        if (!is.null(toRemove)) data <- data[-c(toRemove)]
    }
    return(data)
}

#' Load local files
#'
#' @param folder Character: path to folder containing files of interest
#' @param name Character: name of the category containing all loaded datasets
#' @param ignore Character: skip folders and filenames that match the expression
#' @param verbose Boolean: detail steps while parsing?
#'
#' @importFrom stats setNames
#'
#' @family functions to load local files
#' @family functions to load data
#' @return List of data frames from valid files
#' @export
#'
#' @examples
#' \dontrun{
#' folder <- "~/Downloads/ACC 2016"
#' data <- loadLocalFiles(folder)
#'
#' ignore <- c(".aux.", ".mage-tab.", "junction quantification")
#' loadLocalFiles(folder, ignore)
#' }
loadLocalFiles <- function(folder, ignore=c(".aux.", ".mage-tab."),
                           name="Data", verbose=FALSE) {
    if (!dir.exists(folder)) stop("Folder does not exist.")

    time <- Sys.time()
    # Get all files in the specified directory and subdirectories
    files <- list.files(folder, recursive=TRUE, full.names=TRUE)

    # Exclude undesired subdirectories or files
    files <- files[!dir.exists(files)]
    ignore <- paste(ignore, collapse = "|")
    if (ignore != "") files <- files[!grepl(ignore, files)]

    updateProgress(sprintf("Browsing files in %s...", folder),
                   divisions=length(files))

    loaded <- list()
    formats <- loadFileFormats()
    for (each in seq_along(files)) {
        updateProgress("Processing file", detail = basename(files[each]))
        loadedFile <- suppressWarnings(tryCatch(
            loadFile(files[each], formats, verbose=verbose, multiple=TRUE),
            error=return))
        if (is(loadedFile, "error")) {
            warning(sprintf("Error while reading %s:\n    %s",
                            files[each], loadedFile$message))
        } else if (is.data.frame(loadedFile)) {
            loaded[[length(loaded) + 1]] <- loadedFile
        } else {
            loaded <- append(loaded, loadedFile)
        }
    }
    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- Filter(length, loaded)
    closeProgress()

    if (length(loaded) == 0) {
        compressed <- grep("tar.gz$|tar$|zip$", files, value=TRUE,
                           ignore.case=TRUE)
        anyCompressed <- length(compressed) > 0
        msg <- "No supported files were found in the given folder."
        if (anyCompressed) {
            msg <- paste(msg, "\n\nIf applicable, extract compressed archive",
                         "files (e.g. ZIP and TAR).")
        }
        warning(msg)
        data <- NULL
    } else {
        loaded <- removeRedundantDatasets(loaded)
        loaded <- list(loaded)
        loaded <- loadTCGAsampleMetadata(loaded)

        data <- setNames(loaded, name)
        data <- processDatasetNames(data)
        message("Files loaded in ", format(round(Sys.time() - time, 2)))
    }
    return(data)
}

#' Load local files
#' @inheritParams appServer
#' @param replace Boolean: replace loaded data?
#'
#' @importFrom shinyjs disable enable
#'
#' @inherit psichomics return
#' @keywords internal
setLocalData <- function(input, output, session, replace=TRUE) {
    time <- startProcess("acceptFile")

    folder <- input$localFolder
    category <- input$localCategory
    if (identical(category, "")) category <- "User dataset"
    ignore <- c(".aux.", ".mage-tab.", input$localIgnore)

    # Load valid local files
    data <- tryCatch(loadLocalFiles(folder, name=category, ignore),
                     warning=return, error=return)
    if (is(data, "warning") || is(data, "error")) {
        warningModal(session, "No files available to load", data$message,
                     modalId="localDataModal", caller="Load local data")
    } else if (!is.null(data)) {
        if (!replace) data <- processDatasetNames(c(getData(), data))
        setData(data)
    } else {
        errorModal(session, "No data loaded", "Something went wrong...",
                   modalId="localDataModal", caller="Load local data")
    }
    endProcess("acceptFile", time)
}

#' @rdname setLocalData
setMultipleFilesData <- function(input, output, session, replace=TRUE) {
    category <- input$userFilesCategory
    if (identical(category, "")) category <- "User dataset"

    # Load files
    files <- c("Sample metadata"        =input$sampleInfo,
               "Clinical data"          =input$subjectInfo,
               "Gene expression"        =input$geneExpr,
               "Junction quantification"=input$junctionQuant,
               "Inclusion levels"       =input$ASquant)
    files <- files[files != ""]
    ASquant <- input$ASquant

    # Check if at least one input file was given
    if (length(files) == 0) {
        errorModal(session, "No files selected",
                   "Please provide at least one file.",
                   modalId="localDataModal", caller="Load local data")
        return(NULL)
    } else if (any(!isFile(files))) {
        formatFileInfo <- function(item, files) {
            filepath <- prepareWordBreak(files[[item]])
            tagList(tags$b(names(files[item])), tags$br(),
                    tags$kbd(filepath), tags$br(), tags$br())
        }

        nonExisting   <- files[!isFile(files)]
        filesNotFound <- do.call(
            tagList, lapply(names(nonExisting), formatFileInfo, nonExisting))

        errorModal(session, "Files not found",
                   "The following files were not found:", tags$br(), tags$br(),
                   filesNotFound, modalId="localDataModal",
                   caller="Load local data")
        return(NULL)
    }

    time <- startProcess("loadMultipleFiles")

    loaded <- list()
    allFormats <- loadFileFormats()
    updateProgress("Loading files...", divisions=length(files))
    for (each in seq(files)) {
        file <- files[[each]]
        if (file == "") next

        # Get appropriate format for each type
        dataType <- names(files[each])
        formats <- allFormats[sapply(allFormats, "[[", "dataType") == dataType]

        updateProgress("Processing file", detail=basename(file))
        loadedFile <- tryCatch(loadFile(file, formats, multiple=TRUE),
                               warning=return, error=return)
        if (!is(loadedFile, "warning") && !is(loadedFile, "error")) {
            if (is.data.frame(loadedFile)) {
                loaded[[length(loaded) + 1]] <- loadedFile
            } else {
                loaded <- append(loaded, loadedFile)
            }
        }
    }

    names(loaded) <- sapply(loaded, attr, "tablename")
    loaded <- list(Filter(length, loaded))
    loaded <- loadTCGAsampleMetadata(loaded)

    data <- setNames(loaded, category)
    data <- processDatasetNames(data)

    if (!is.null(data)) {
        if(replace) {
            setData(data)
        } else {
            data <- processDatasetNames(c(getData(), data))
            setData(data)
        }
        if (all(ASquant != "")) setCategory(category)
    }
    endProcess("loadMultipleFiles", time)
}

#' @rdname appServer
#' @importFrom shiny updateTextInput
localDataServer <- function(input, output, session) {
    ns <- session$ns

    ## Load data based on individual files
    prepareFileBrowser(session, input, "sampleInfo")
    prepareFileBrowser(session, input, "subjectInfo")
    prepareFileBrowser(session, input, "geneExpr")
    prepareFileBrowser(session, input, "junctionQuant")
    prepareFileBrowser(session, input, "ASquant")

    ## Load data based on folder
    prepareFileBrowser(session, input, "localFolder", directory=TRUE)

    # # The button is only enabled if it meets the conditions that follow
    # observe(toggleState("acceptFile", input$species != ""))

    # Update category name input based on given folder
    observe({
        folder <- input$localFolder
        if (!is.null(folder))
            updateTextInput(session, "localCategory", value=basename(folder))
    })

    # If data is loaded, let user replace or append to loaded data
    observeEvent(input$loadMultipleFiles, {
        if (!is.null(getData())) {
            loadedDataModal(session,
                            "localDataModal",
                            "multipleFilesReplace",
                            "multipleFilesLocalAppend")
        } else {
            setMultipleFilesData(input, output, session)
        }
    })

    # Load data when the user presses to replace data
    observeEvent(input$multipleFilesReplace,
                 setMultipleFilesData(input, output, session, replace=TRUE))

    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$multipleFilesLocalAppend,
                 setMultipleFilesData(input, output, session, replace=FALSE))

    # If data is loaded, let user replace or append to loaded data
    observeEvent(input$acceptFile, {
        folder <- input$localFolder
        if (!dir.exists(folder)) {
            # Folder not found
            errorModal(session, "Folder not found",
                       "Check if the folder path is correct.",
                       modalId="localDataModal", caller="Load local data")
            enable("acceptFile")
            return(NULL)
        } else if (!is.null(getData())) {
            loadedDataModal(session,
                            "localDataModal",
                            "localReplace",
                            "localAppend")
        } else {
            setLocalData(input, output, session)
        }
    })

    # Load data when the user presses to replace data
    observeEvent(input$localReplace,
                 setLocalData(input, output, session, replace=TRUE))

    # Load data when the user presses to load new data (keep previously loaded)
    observeEvent(input$localAppend,
                 setLocalData(input, output, session, replace=FALSE))
}

attr(localDataUI, "loader") <- "data"
attr(localDataServer, "loader") <- "data"
