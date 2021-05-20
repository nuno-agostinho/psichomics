## Auxiliary functions used throughout the program

#' Print startup message
#'
#' @param pkgname Character: package name
#' @param libname Character: library name
#'
#' @importFrom utils packageVersion
#'
#' @return Startup message
.onAttach <- function(libname, pkgname) {
    version <- utils::packageVersion(pkgname, libname)
    msg <- sprintf(
        "psichomics %s: start the visual interface by running psichomics()",
        version)
    site <- "https://nuno-agostinho.github.io/psichomics"
    msg2 <- paste("Full documentation and tutorials at", site)
    packageStartupMessage(paste("", msg, msg2, sep = "\n"))
}

#' Preserve attributes when extracting values
#'
#' Add object to class \code{sticky}
#'
#' @param x Object
#'
#' @return Object with class \code{sticky}
#'
#' @keywords internal
preserveAttributes <- function(x) {
    class(x) <- c("sticky", class(x))
    return(x)
}

#' Preserve attributes of \code{sticky} objects when extracting or transposing
#' object
#'
#' Most attributes - with the exception of \code{names}, \code{dim},
#' \code{dimnames}, \code{class} and \code{row.names} - are preserved in simple
#' transformations of objects from class \code{sticky}
#'
#' @return Transformed object with most attributes preserved
#'
#' @export
t.sticky <- function(x) {
    attrs <- attributes(x)
    attrs <- attrs[!names(attrs) %in% c("names", "dim", "dimnames", "class",
                                        "row.names")]
    res <- NextMethod()
    res <- addObjectAttrs(res, attrs)

    tmp <- attr(res, "rowData")
    attr(res, "rowData") <- attr(res, "colData")
    attr(res, "colData") <- tmp

    if (!is(res, "sticky")) res <- preserveAttributes(res)
    return(res)
}

#' @rdname t.sticky
#'
#' @param x Object
#' @param i,j,... Numeric or character: indices of elements to extract
#'
#' @export
`[.sticky` <- function(x, i, j, ...) {
    attrs <- attributes(x)
    attrs <- attrs[!names(attrs) %in% c("names", "dim", "dimnames", "class",
                                        "row.names")]
    res <- NextMethod()
    res <- addObjectAttrs(res, attrs)

    # Subset rowData and colData based on row and column names, respectively
    subsetInfo <- function(data, names) {
        areValidNames <- !is.null(data) && !is.null(names) &&
            all(names %in% rownames(data))
        if (areValidNames) data <- data[names, ]
        return(data)
    }

    attr(res, "rowData") <- subsetInfo(attr(res, "rowData"), rownames(res))
    attr(res, "colData") <- subsetInfo(attr(res, "colData"), colnames(res))

    if (!is(res, "sticky")) res <- preserveAttributes(res)
    return(res)
}

#' Round down/up the minimum/maximum value
#' @param x Numeric: values
#' @param digits Numeric: number of maximum digits
#'
#' @return Rounded numeric value
#' @keywords internal
roundMinDown <- function(x, digits=0) floor  (min(x) * 10^digits) / 10^digits

#' @rdname roundMinDown
roundMaxUp   <- function(x, digits=0) ceiling(max(x) * 10^digits) / 10^digits

#' Get psichomics file inside a given directory
#' @inheritParams base::system.file
#' @return Loaded file
#' @keywords internal
insideFile <- function(...) {
    return(system.file(..., package="psichomics"))
}

#' Check if files exist
#'
#' @param files Character: vector of filepaths to check
#'
#' @return Boolean vector stating whether each file exists or not
#' @keywords internal
isFile <- function(files) {
    fileExists <- file.exists(files) & !dir.exists(files)
    names(fileExists) <- files
    return(fileExists)
}

#' Load psichomics-specific file
#' @param file Character: path to the file
#'
#' @return Loaded file
#' @export
#'
#' @examples
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
readFile <- function(file) {
    readRDS(insideFile("extdata", file))
}

#' Trims whitespace from a word
#'
#' @param word Character to trim
#'
#' @return Character without whitespace
#' @keywords internal
#'
#' @examples
#' psichomics:::trimWhitespace("    hey   there     ")
#' psichomics:::trimWhitespace(c("pineapple    ", "one two three",
#'                               " sunken    ship   "))
trimWhitespace <- function(word) {
    # Remove leading and trailing whitespace
    word <- gsub("^\\s+|\\s+$", "", word)
    # Replace multiple spaces between words with one single space
    word <- gsub("\\s+", " ", word)
    return(word)
}

#' Filter \code{NULL} elements from a vector or a list
#'
#' @param v Vector or list
#'
#' @return Filtered vector or list with no \code{NULL} elements; if \code{v} is
#' a vector composed of \code{NULL} elements, returns a \code{NULL}; if \code{v}
#' is a list of \code{NULL} elements, returns an empty list
#' @keywords internal
rm.null <- function(v) Filter(Negate(is.null), v)

#' Escape symbols for use in regular expressions
#'
#' @param ... Characters to be pasted with no space
#'
#' @return Escaped string
#' @keywords internal
escape <- function(...) {
    # return(gsub("/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]", "\\$&", string))
    return(gsub("(\\W)", "\\\\\\1", paste0(...)))
}

#' Check if a number is whole
#'
#' @param x Object to be tested
#' @param tol Numeric: tolerance used for comparison
#'
#' @return TRUE if number is whole; otherwise, FALSE
#' @keywords internal
is.whole <- function(x, tol=.Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}

# Template for using faster version of row- or column-wise statistics from Rfast
fasterDataStats <- function(mat, fastFUN, classicFUN=NULL, na.rm=FALSE,
                            byRow=TRUE, ...) {
    if (is.vector(mat)) {
        mat <- matrix(mat, nrow=1)
    } else if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
    }

    if (!is.null(classicFUN)) {
        # Use classicFUN to calculate missing values
        res <- fastFUN(mat, ...)
        if (na.rm) {
            nas      <- is.na(res)
            res[nas] <- classicFUN(mat[nas, , drop=FALSE], na.rm=na.rm,
                                   fast=FALSE)
        }
    } else {
        res <- fastFUN(mat, na.rm=na.rm, ...)
    }

    if (byRow) {
        ns <- rownames
    } else {
        ns <- colnames
    }
    if (!is.null(ns(mat))) names(res) <- ns(mat)
    return(res)
}

#' Calculate statistics for each row or column of a matrix
#'
#' @param mat Matrix
#' @param na.rm Boolean: remove missing values (\code{NA})?
#' @param fast Boolean: use \code{Rfast} functions? They may return different
#' results from R built-in functions
#'
#' @importFrom Rfast rowmeans
#'
#' @return Vector of selected statistic
#' @keywords internal
#'
#' @examples
#' df <- rbind("Gene 1"=c(3, 5, 7), "Gene 2"=c(8, 2, 4), "Gene 3"=c(9:11))
#' psichomics:::customRowMeans(df)
#' psichomics:::customRowVars(df, fast=TRUE)
customRowMeans <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowMeans, customRowMeans, na.rm=na.rm)
    } else if ( !is.null(dim(mat)) ) {
        nas <- 0
        if (na.rm) nas <- rowSums(is.na(mat))
        res <- rowSums(mat, na.rm=na.rm) / (ncol(mat) - nas)
    } else {
        res <- mean(mat, na.rm=na.rm)
    }
    return(res)
}

#' @rdname customRowMeans
#' @importFrom Rfast rowMedians
customRowMedians <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowMedians, classicFUN=NULL, na.rm=na.rm)
    } else {
        res <- apply(mat, 1, median, na.rm=na.rm)
    }
    return(res)
}

#' @rdname customRowMeans
#' @importFrom Rfast rowVars
customRowVars <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowVars, classicFUN=NULL, na.rm=na.rm)
        res[res < 0] <- 0 # Safeguard for when log-transforming values
    } else if ( !is.null(dim(mat)) ) {
        means      <- customRowMeans(mat, na.rm=na.rm, fast=FALSE)
        meansSqDev <- (mat - means) ** 2
        squaresSum <- rowSums(meansSqDev, na.rm=na.rm)

        nas <- 0
        if (na.rm) nas <- rowSums(is.na(mat))
        dem <- ncol(mat) - nas - 1
        res <- squaresSum / dem
    } else {
        res <- var(mat, na.rm=na.rm)
    }
    return(res)
}

#' @rdname customRowMeans
#' @importFrom Rfast rowMins
customRowMins <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowMins, customRowMins, na.rm=na.rm,
                               value=TRUE)
    } else {
        suppressWarnings(apply(mat, 1, min, na.rm=na.rm))
    }
}

#' @rdname customRowMeans
#' @importFrom Rfast rowMaxs
customRowMaxs <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowMaxs, customRowMaxs, na.rm=na.rm,
                               value=TRUE)
    } else {
        res <- suppressWarnings(apply(mat, 1, max, na.rm=na.rm))
    }
    return(res)
}

#' @rdname customRowMeans
#' @importFrom Rfast rowrange
customRowRanges <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, rowrange, customRowRanges, na.rm=na.rm)
    } else {
        maxs <- customRowMaxs(mat, na.rm=na.rm, fast=FALSE)
        mins <- customRowMins(mat, na.rm=na.rm, fast=FALSE)
        res  <- maxs - mins
    }
    return(res)
}

#' @rdname customRowMeans
#' @importFrom Rfast colMedians
customColMedians <- function(mat, na.rm=FALSE, fast=FALSE) {
    if (fast) {
        res <- fasterDataStats(mat, colMedians, classicFUN=NULL, na.rm=na.rm,
                               byRow=FALSE)
    } else {
        res <- apply(mat, 2, median, na.rm=na.rm)
    }
    return(res)
}

#' Rename vector to avoid duplicated values with another vector
#'
#' Renames values by adding an index to the end of duplicates. This allows to
#' prepare unique values in two vectors before a merge, for instance.
#'
#' @param check Character: values to rename if duplicated
#' @param comp Character: values to compare with
#'
#' @importFrom utils head
#'
#' @return Character vector with renamed values if duplicated; else, it
#' returns the usual values. It does not return the comparator values.
#' @keywords internal
#'
#' @examples
#' psichomics:::renameDuplicated(check = c("blue", "red"), comp = c("green",
#'                                                                  "blue"))
renameDuplicated <- function(check, comp) {
    # If there's nothing to compare with, return the values
    # if (length(comp) == 0) return(check)

    repeated <- check %in% comp | duplicated(check)
    uniq <- check[!repeated]

    for (dup in which(repeated)) {
        # Locate matches (don't forget the counter)
        all <- c(comp, head(check[seq(dup)], -1))
        expr <- paste0(escape(check[dup]), " \\([0-9]+\\)|", escape(check[dup]))
        locate <- grep(expr, all, value = TRUE)

        # Get the maximum counter and add one
        counter <- sub(".* \\(([0-9]+)\\)", "\\1", locate)

        # Replace strings with 0
        counter[grep("^[0-9]*$", counter, invert =TRUE)] <- 0
        check[dup] <- sprintf("%s (%i)", check[dup],
                              max(as.numeric(counter)) + 1)
    }
    return(check)
}

#' Get subjects from given samples
#'
#' @param sampleId Character: sample identifiers
#' @param patientId Character: subject identifiers to filter by (optional; if a
#' matrix or data frame is given, its rownames will be used to infer the subject
#' identifiers)
#' @param na Boolean: return \code{NA} for samples with no matching subjects
#' @param sampleInfo Data frame or matrix: sample information containing the
#' sample identifiers as rownames and a column named "Subject ID" with the
#' respective subject identifiers
#'
#' @aliases getPatientFromSample
#' @family functions for data grouping
#' @return Character: subject identifiers corresponding to the given samples
#'
#' @export
#' @examples
#' samples <- paste0("GTEX-", c("ABC", "DEF", "GHI", "JKL", "MNO"), "-sample")
#' getSubjectFromSample(samples)
#'
#' # Filter returned samples based on available subjects
#' subjects <- paste0("GTEX-", c("DEF", "MNO"))
#' getSubjectFromSample(samples, subjects)
getSubjectFromSample <- function(sampleId, patientId=NULL, na=FALSE,
                                 sampleInfo=NULL) {
    if (!is.null(patientId) &&
        (is.matrix(patientId) || is.data.frame(patientId))) {
        patientId <- rownames(patientId)
    }

    # Extract subject identifiers from sample ID and then retrieve their index
    extractSubjectIndex <- function(pattern, samples, allSubjects) {
        subject <- gsub(pattern, "\\1", samples)
        names(subject) <- samples

        if (!is.null(allSubjects)) {
            # Filter by subjects of interest
            if (na) {
                # Return NA as the corresponding subject
                subject[!subject %in% allSubjects] <- NA
            } else {
                subject <- subject[subject %in% allSubjects]
            }
        }
        return(subject)
    }

    if ( any(grepl("^TCGA", sampleId)) ) {
        # Retrieve TCGA subject index
        extractSubjectIndex("(TCGA-.*?-.*?)-.*", sampleId, patientId)
    } else if ( any(grepl("^GTEX", sampleId)) ) {
        # Retrieve GTEx subject index
        extractSubjectIndex("(GTEX-.*?)-.*", sampleId, patientId)
    } else if ( "Subject ID" %in% colnames(sampleInfo) ) {
        # Based on user-provided files
        subjects <- as.character(sampleInfo[ , "Subject ID"])
        names(subjects) <- rownames(sampleInfo)
        return(subjects)
    } else {
        return(NULL)
    }
}

#' @export
getPatientFromSample <- getSubjectFromSample

#' Get samples matching the given subjects
#'
#' @param patients Character or list of characters: subject identifiers
#' @param samples Character: sample identifiers
#' @param clinical Data frame or matrix: clinical dataset
#' @param rm.NA Boolean: remove missing values?
#' @param match Integer: vector of subject index with the sample identifiers as
#' name to save time (optional)
#' @param showMatch Boolean: show matching subject index?
#'
#' @aliases getSampleFromPatient getMatchingSamples
#' @family functions for data grouping
#' @return Names of the matching samples (if \code{showMatch = TRUE},
#' a character with the subjects as values and their respective samples as names
#' is returned)
#' @export
#'
#' @examples
#' subjects <- c("GTEX-ABC", "GTEX-DEF", "GTEX-GHI", "GTEX-JKL", "GTEX-MNO")
#' samples <- paste0(subjects, "-sample")
#' clinical <- data.frame(samples=samples)
#' rownames(clinical) <- subjects
#' getSampleFromSubject(subjects[c(1, 4)], samples, clinical)
getSampleFromSubject <- function(patients, samples, clinical=NULL, rm.NA=TRUE,
                               match=NULL, showMatch=FALSE) {
    if (is.null(match))
        match <- getSubjectFromSample(samples, clinical)

    if (is.list(patients)) {
        samples <- lapply(patients, function(i) {
            res <- match[match %in% i]
            if (!showMatch) res <- unique(names(res))
            if (rm.NA) res <- res[!is.na(res)]
            return(res)
        })
    } else {
        samples <- match[match %in% patients]
        if (!showMatch) samples <- unique(names(samples))
        if (rm.NA) samples <- samples[!is.na(samples)]
    }
    return(samples)
}

#' @export
getSampleFromPatient <- getSampleFromSubject

#' @export
getMatchingSamples <- getSampleFromSubject

#' Assign one group to each element
#'
#' @param groups List of integers: groups of elements
#' @param elem Character: all elements available
#' @param outerGroupName Character: name to give to outer group (if \code{NULL},
#' only show elements matched to their respective groups)
#'
#' @family functions for data grouping
#' @return Character vector where each element corresponds to the group of the
#' respective element
#' @export
#'
#' @examples
#' groups <- list(1:3, 4:7, 8:10)
#' names(groups) <- paste("Stage", 1:3)
#' groupPerElem(groups)
groupPerElem <- function(groups, elem=NULL, outerGroupName=NA) {
    if (length(groups) == 0) {
        singleGroup <- "Single group"
        if (length(elem) > 0)
            return(rep(singleGroup, length(elem)))
        else
            return(singleGroup)
    }

    all <- unlist(groups)
    names(all) <- rep(names(groups), sapply(groups, length))

    finalGroups <- NULL
    if (!is.null(elem) && !is.null(outerGroupName)) {
        finalGroups <- rep(outerGroupName, length(elem))
        names(finalGroups) <- elem
    }

    colour <- NULL
    assignedColours <- attr(groups, "Colour")
    for (each in unique(all)) {
        each <- as.character(each) # Force to use numeric identifiers as names
        groupName <- names(all[all == each])
        groupNameStr <- paste(groupName, collapse=", ")
        if (!is.null(assignedColours) && !groupNameStr %in% names(colour)) {
            cols   <- assignedColours[groupName]
            colour <- c(colour, setNames(Reduce(blendColours, cols),
                                         groupNameStr))
        }
        finalGroups[each] <- groupNameStr
    }
    attr(finalGroups, "Colour") <- colour
    return(finalGroups)
}

#' @rdname missingDataModal
#'
#' @param modal Character: modal identifier
#'
#' @keywords internal
loadRequiredData <- function( modal=NULL ) {
    modal <- ifelse(is.null(modal), "null", modal)
    return(sprintf("showDataPanel('#%s');", modal))
}

#' Get the path to the Downloads folder
#'
#' @family functions associated with TCGA data retrieval
#' @family functions associated with GTEx data retrieval
#' @family functions associated with SRA data retrieval
#' @return Path to Downloads folder
#' @export
#'
#' @examples
#' getDownloadsFolder()
getDownloadsFolder <- function() {
    if (Sys.info()['sysname'] == "Windows")
        folder <- dirname("~")
    else
        folder <- path.expand("~")
    folder <- file.path(folder, "Downloads")
    folder <- paste0(folder, .Platform$file.sep)
    return(folder)
}

#' Get number of significant digits
#'
#' @param n Numeric: number to round
#'
#' @importFrom shiny isolate
#'
#' @return Formatted number with a given number of significant digits
#' @keywords internal
signifDigits <- function(n) {
    return(isolate(formatC(n, getSignificant(), format="g")))
}

#' Round by the given number of digits
#'
#' @param n Numeric: number to round
#'
#' @return Formatted number with a given numeric precision
#' @keywords internal
roundDigits <- function(n) {
    return(isolate(formatC(n, getPrecision(), format="f")))
}

#' Blend two HEX colours
#'
#' @param colour1 Character: HEX colour
#' @param colour2 Character: HEX colour
#' @param colour1Percentage Character: percentage of colour 1 mixed in blended
#' colour
#'
#' @source Code modified from \url{https://stackoverflow.com/questions/5560248}
#'
#' @return Character representing an HEX colour
#' @keywords internal
#'
#' @examples
#' psichomics:::blendColours("#3f83a3", "#f48000")
blendColours <- function (colour1, colour2, colour1Percentage=0.5) {
    colour1 <- gsub("#", "", colour1)
    colour1 <- as.hexmode(colour1)
    R1 <- bitwShiftR(colour1, 16)
    G1 <- bitwAnd(bitwShiftR(colour1, 8), 0x00FF)
    B1 <- bitwAnd(colour1, 0x0000FF)

    colour2 <- gsub("#", "", colour2)
    colour2 <- as.hexmode(colour2)
    R2 <- bitwShiftR(colour2, 16)
    G2 <- bitwAnd(bitwShiftR(colour2, 8), 0x00FF)
    B2 <- bitwAnd(colour2, 0x0000FF)

    # Round to biggest integer if ending in .5
    mround <- function(x) trunc(x + 0.5)

    red   <- 0x1000000 + (mround((R2 - R1) * colour1Percentage) + R1) * 0x10000;
    green <- (mround((G2 - G1) * colour1Percentage) + G1) * 0x100;
    blue  <- mround((B2 - B1) * colour1Percentage) + B1;
    blended <- substr(as.hexmode(red + green + blue), 2, 16)
    return(paste0("#", blended));
}

#' Plot survival curves
#'
#' @param object \code{survfit} object as returned from
#' \code{\link{survfit.survTerms}()} function
#' @inheritDotParams highcharter::hc_add_series -hc -data
#' @param fun Name of function or function used to transform the survival curve:
#' \code{log} will put y axis on log scale, \code{event} plots cumulative events
#' (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard function (f(y) =
#' -log(y)), and \code{cloglog} creates a complimentary log-log survival plot
#' (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param markTimes Label curves marked at each censoring time?
#' @param symbol Symbol to use as marker
#' @param markerColor Colour of the marker; if \code{NULL}, the respective
#' colour of each series are used
#' @param ranges Plot interval ranges?
#' @param rangesOpacity Opacity of the interval ranges
#'
#' @importFrom highcharter %>% hc_add_series highchart hc_tooltip hc_yAxis
#' hc_plotOptions fa_icon_mark JS
#' @importFrom stats setNames
#'
#' @return \code{highchart} object to plot survival curves
#' @keywords internal
#'
#' @examples
#'
#' # Plot Kaplan-Meier curves
#' require("survival")
#' require("highcharter")
#' leukemia.surv <- survfit(Surv(time, status) ~ x, data = aml)
#' hchart(leukemia.surv)
#'
#' # Plot the cumulative hazard function
#' lsurv2 <- survfit(Surv(time, status) ~ x, aml, type='fleming')
#' hchart(lsurv2, fun="cumhaz")
#'
#' # Plot the fit of a Cox proportional hazards regression model
#' fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian)
#' ovarian.surv <- survfit(fit, newdata=data.frame(age=60))
#' hchart(ovarian.surv, ranges = TRUE)
hchart.survfit <- function(object, ..., fun = NULL, markTimes = TRUE,
                           symbol = "plus", markerColor = "black",
                           ranges = FALSE, rangesOpacity = 0.3) {
    groups <- NULL
    # Check if there are groups
    if (is.null(object$strata))
        strata <- c("Series 1" = length(object$time))
    else
        strata <- object$strata

    # Modify data according to functions (adapted from survival:::plot.survfit)
    if (is.character(fun)) {
        tfun <- switch(fun,
                       log = function(x) x,
                       event = function(x) 1 - x,
                       cumhaz = function(x) -log(x),
                       cloglog = function(x) log(-log(x)),
                       pct = function(x) x * 100,
                       logpct = function(x) 100 * x,
                       identity = function(x) x,
                       function(x) x)
    } else if (is.function(fun)) {
        tfun <- fun
    } else {
        tfun <- function(x) x
    }

    firsty <- tfun(1)
    object$surv <- tfun(object$surv)
    if (ranges && !is.null(object$upper)) {
        object$upper <- tfun(object$upper)
        object$lower <- tfun(object$lower)
    }

    # Prepare data
    data <- data.frame(x=object$time, y=object$surv,
                       up=object$upper, low=object$lower,
                       group=rep(names(strata), strata),
                       stringsAsFactors = FALSE)
    # Data markers
    marker <- list(list(fillColor=markerColor, symbol=symbol, enabled=TRUE))
    if(markTimes)
        mark <- object$n.censor == 1
    else
        mark <- FALSE

    # Adjust Y axis range
    yValues <- object$surv
    ymin <- ifelse(min(yValues) >= 0, 0, min(yValues))
    ymax <- ifelse(max(yValues) <= 1, 1, max(yValues))

    hc <- highchart() %>%
        hc_tooltip(pointFormat="{point.y}") %>%
        hc_yAxis(min=ymin, max=ymax) %>%
        hc_plotOptions(line = list(marker = list(enabled = FALSE)))

    count <- 0

    # Process groups by columns (CoxPH-like) or in a single column
    if(!is.null(ncol(object$surv))) {
        groups <- seq(ncol(object$surv))
    } else {
        groups <- names(strata)
    }

    summ <- summary(object)$table
    for (name in groups) {
        if (!is.null(ncol(object$surv))) {
            df <- df[c("x", paste(c("y", "low", "up"), col, sep="."))]
            names(df) <- c("x", "y", "low", "up")
            submark <- mark
        } else {
            this <- data$group == name
            df <- data[this, ]
            submark <- mark[this]
        }

        # Add first value if there is no value for time at 0 in the data
        first <- NULL
        if (!0 %in% df$x) first <- list(list(x=0, y=firsty))

        # Mark events
        ls <- lapply(seq(nrow(df)), function(i) as.list(df[i, , drop=FALSE]))
        if (markTimes) ls[submark] <- lapply(ls[submark], c, marker=marker)

        if (is.matrix(summ))
            curveSumm <- summ[name, ]
        else
            curveSumm <- summ

        # Prepare group colours
        colour <- attr(object, "Colour")
        if (!is.null(colour))
            colour <- unname(colour[name])
        else
            colour <- JS("Highcharts.getOptions().colors[", count, "]")

        hc <- do.call(hc_add_series, c(list(
            hc, data=c(first, ls), step="left", name=name, zIndex=1,
            color=colour, ...), curveSumm))

        if (ranges && !is.null(object$upper)) {
            # Add interval range
            range <- lapply(ls, function(i)
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=rangesOpacity,
                lineWidth=0, color=colour, ...)
        }
        count <- count + 1
    }

    return(hc)
}

#' Check unique rows of a data frame based on a set of its columns
#'
#' @param data Data frame or matrix
#' @param ... Name of columns
#'
#' @return Data frame with unique values based on set of columns
#' @keywords internal
uniqueBy <- function(data, ...) {
    sub <- subset(data, select=c(...))
    uniq <- !duplicated(sub)
    return(data[uniq, ])
}

#' Add an exporting feature to a \code{highcharts} object
#'
#' @param hc A \code{highcharts} object
#' @param fill Character: colour fill
#' @param text Character: button text
#'
#' @importFrom highcharter hc_exporting JS
#'
#' @return A \code{highcharts} object with an export button
#' @keywords internal
export_highcharts <- function(hc, fill="transparent", text="Export") {
    createJSExport <- function(type) {
        JS(paste0("function () { this.exportChart({ type: '", type, "' }); }"))
    }

    export <- list(
        list(text="PNG image",  onclick=createJSExport("image/png")),
        list(text="JPEG image", onclick=createJSExport("image/jpeg")),
        list(text="SVG vector image", onclick=createJSExport("image/svg+xml")),
        list(text="PDF document", onclick=createJSExport("application/pdf")),
        list(separator=TRUE),
        list(text="CSV document",
             onclick=JS("function () { this.downloadCSV(); }")),
        list(text="XLS document",
             onclick=JS("function () { this.downloadXLS(); }")))

    hc_exporting(hc, enabled=TRUE, formAttributes=list(target="_blank"),
                 buttons=list(contextButton=list(text=text, menuItems=export,
                                                 theme=list(fill=fill))))
}

#' Create scatter plot
#'
#' Create a scatter plot using \code{highcharter}
#'
#' @param hc \code{Highchart} object
#' @param x Numeric: X axis
#' @param y Numeric: Y axis
#' @param z Numeric: Z axis to set the bubble size (optional)
#' @param label Character: data label for each point (optional)
#' @param showInLegend Boolean: show the data in the legend box?
#' @param color Character: series colour
#' @inheritDotParams highcharter::hc_add_series -hc -data
#'
#' @importFrom highcharter hc_add_series list_parse
#'
#' @return \code{highcharter} object containing information for a scatter plot
#' @keywords internal
hc_scatter <- function (hc, x, y, z=NULL, label=NULL, showInLegend=FALSE,
                        color=NULL, ...) {
    df <- data.frame(x, y)
    if (!is.null(z)) df <- cbind(df, z=z)
    if (!is.null(label)) df <- cbind(df, label=label)

    args <- list(...)
    for (i in seq_along(args)) {
        if (!is.list(args[[i]]) && length(x) == length(args[[i]]) &&
            names(args[i]) != "name") {
            attr <- list(args[i])
            names(attr) <- names(args)[i]
            df <- cbind(df, attr)
            args[[i]] <- character(0)
        }
    }

    args <- Filter(length, args)
    ds <- list_parse(df)
    type <- ifelse(!is.null(z), "bubble", "scatter")

    if (!is.null(label))
        dlopts <- list(enabled=TRUE, format="{point.label}")
    else
        dlopts <- list(enabled=FALSE)

    do.call("hc_add_series", c(list(hc, data=ds, type=type, color=color,
                                    showInLegend=showInLegend,
                                    dataLabels=dlopts), args))
}

#' Check if running in RStudio Server
#'
#' @return Boolean stating whether running in RStudio Server
#' @keywords internal
isRStudioServer <- function() {
    tryCatch(
        rstudioapi::isAvailable() && rstudioapi::versionInfo()$mode == "server",
        error=function(i) FALSE)
}
