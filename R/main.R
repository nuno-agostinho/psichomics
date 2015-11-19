# Load PSI files
# source("R/readPSIfiles.R")

# What is important to keep from each event?
# Attributes of interest:
#   chromosome
#   strand
#   gene / mRNA
#   event type
#   psi
#   exon 1 to 4 (start, end)

parseSuppaEventA <- function(event) {
  # create vector in the end resorting to many applies
  
  # split event ID by semicolon and colon symbols
  tmp <- strsplit(event, ";|:")
  gene <- sapply(tmp, "[[", 1)
  event_type <- sapply(tmp, "[[", 2)
  chromosome <- sapply(tmp, "[[", 3)
  strand <- sapply(tmp, function(x) x[[length(x)]])
  
  # split string by minus symbol (which indicates the span of an exon)
  # tmp <- strsplit(tmp, "-")
  return(c(gene = gene,
           event_type = event_type,
           chromosome = chromosome,
           strand = strand))
}

parseSuppaEventB <- function(event) {
  # create vector in the end resorting to a for loop
  
  # split event ID by semicolon and colon symbols
  tmp <- strsplit(event, ";|:")
  gene <- character()
  event_type <- character()
  chromosome <- character()
  strand <- character()
  for (each in tmp) {
    gene <- c(gene, each[[1]])
    event_type <- c(event_type, each[[2]])
    chromosome <- c(chromosome, each[[3]])
    
    len <- length(each)
    strand <- c(strand, each[[len]])
  }
  return(gene)
  
  # split string by minus symbol (which indicates the span of an exon)
  # tmp <- strsplit(tmp, "-")
}

parseSuppaEventC <- function(event) {
  # create vector in the end resorting to a single apply
  
  # split event ID by semicolon and colon symbols
  tmp <- strsplit(event, ";|:")
  gene <- character()
  event_type <- character()
  chromosome <- character()
  strand <- character()
  f <- lapply (tmp, function(each) {
    gene <- c(gene, each[[1]])
    event_type <- c(event_type, each[[2]])
    chromosome <- c(chromosome, each[[3]])
    
    len <- length(each)
    strand <- c(strand, each[[len]])
  }) # lapply
  return(f)
  
  # split string by minus symbol (which indicates the span of an exon)
  # tmp <- strsplit(tmp, "-")
}

parseSuppaEventD <- function(event) {
  # split event ID by semicolon and colon symbols
  tmp <- strsplit(event, ";|:")
  a <- t(data.frame(tmp, row.names = c("gene", "event_type", "chromosome",
                                       "exon 1", "exon 2", "strand")))
  row.names(a) <- NULL
  
  # split string by minus symbol (which indicates the span of an exon)
  # tmp <- strsplit(tmp, "-")
}

#' Parses an event ID from SUPPA
#'
#' @param event Event ID
#'
#' @return List with the event attributes (chromosome, strand, event type and
#' the position of the exon boundaries)
#' @export
#'
#' @examples
#' event <- "ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-"
#' parseSuppaEventE(event)
#' 
#' events <- c("ENSG00000000419;A3:20:49557492-49557642:49557470-49557642:-",
#'             "ENSG00000000419;A3:20:49557492-49558568:49557470-49558568:-",
#'             "ENSG00000000457;A3:1:69838269-69839396:69838180-69839396:-")
#' t(sapply(events, parseSuppaEventE))
parseSuppaEventE <- function(event) {
  # Split event ID by semicolons and colons
  tmp <- strsplit(event, ";|:")[[1]]
  len <- length(tmp)
  
  # Get the attributes of the event
  event_attrs <- as.list(tmp[c(1:3, len)])
  names(event_attrs) <- c("Gene", "Event type", "Chromosome", "Strand")
  
  # Get the junction positions for each exon and parse them
  junctions <- tmp[4:(len-1)]
  parsed_junctions <- parseSuppaJunctions(event_attrs[["Event type"]],
                                          event_attrs[["Strand"]],
                                          junctions)
  return(c(event_attrs, parsed_junctions))
}

#' Parses junctions
#'
#' @param event_type 
#' @param junctions 
#'
#' @return
#' @export
#'
#' @examples
parseSuppaJunctions <- function(event_type, strand, junctions) {
  # Split junctions by the hyphen
  junctions <- strsplit(junctions, "-")
  junctions <- unlist(junctions)       
  
  # If minus strand, reverse junctions
  if(strand == "-") junctions <- rev(junctions)
  
  # Fill list of parsed junctions with NAs
  parsed_junctions <- as.list(rep(NA, 8))
  names(parsed_junctions) <- c("C1 start",
                               "C1 end",
                               "A1 start",
                               "A1 end",
                               "A2 start",
                               "A2 end",
                               "C2 start",
                               "C2 end"
  )
  
  # Parse junction positions according to event type
  parseJunctionsByType <- switch(
    event_type,
    "A3" = parseSuppaJunctionsA3,
    "A5" = parseSuppaJunctionsA5,
    "SE" = parseSuppaJunctionsSE,
    "MX" = parseSuppaJunctionsMX,
    "RI" = parseSuppaJunctionsRI,
    "AF" = parseSuppaJunctionsAF,
    "AL" = parseSuppaJunctionsAL
  )
  parsed_junctions <- parseJunctionsByType(junctions, parsed_junctions)
  return(parsed_junctions)
}

parseSuppaJunctionsA3 <- function(junctions, parsed_junctions) {
  parsed_junctions[["C1 end"]] <- junctions[1]
  parsed_junctions[["C2 start"]] <- c(junctions[c(2,4)])
  return(parsed_junctions)
}

#### BENCHMARKING
benchmarking <- function() {
  library(rbenchmark)
  events <- rownames(suppa_A3.tab)
  event <- events[1]
  
  ## Different methods os parsing events
  benchmark(replications = 10000,
            parseSuppaEventA(event), # 2.589s for 10000 reps of 1 event
            parseSuppaEventB(event), # 0.482s for 10000 reps of 1 event
            parseSuppaEventD(event), # 9.424s for 10000 reps of 1 event
            parseSuppaEventE(event)) # 0.541s for 10000 reps of 1 event
  
  ## Different loop techniques to process many events
  c <- character()
  benchmark(replications = 1,
            for (e in events)
              t <-c(t,parseSuppaEventB(e)),  # 62.209s for 16861 events
            sapply(events, parseSuppaEventB, 
                   USE.NAMES = FALSE),       #  0.969s for 16861 events
            parseSuppaEventB(events) )       # 10.572s for 16861 events
  
  benchmark(replications = 1,
            sapply(events, parseSuppaEventB, # 0.833s for 16861 events
                   USE.NAMES = FALSE),
            sapply(events, parseSuppaEventE, # 0.908s for 16861 events
                   USE.NAMES = FALSE))
  
  ## Timing to get the tranpose of interest
  benchmark(replications = 1,
            a <- sapply(events, parseSuppaEventE),
            t(a))
  
  ## Timing sorting methods
  tmp <- strsplit(event, ";|:")[[1]]
  junctions <- tmp[4:(len-1)]
  junctions <- strsplit(junctions, "-")
  
  revsort() <- function() {
    lapply(junctions, rev)
    rev(junctions) }
  shellsort() <- function() {
    lapply(junctions, sort, method = "shell", decreasing = TRUE)
    rev(junctions) }
  quicksort() <- function() {
    lapply(junctions, sort, method = "quick", decreasing = TRUE)
    rev(junctions) }
  
  benchmark(replications = 100000, revsort(), shellsort(), quicksort())
}