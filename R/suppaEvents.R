# Load PSI files
# source("R/readPSIfiles.R")

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

#' Parses a splicing event ID from SUPPA (more information available at
#' https://bitbucket.org/regulatorygenomicsupf/suppa)
#'
#' @param event Chracter: Splicing event ID
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
#'             "ENSG00000000003;A5:X:99890743-99891188:99890743-99891605:-")
#' lapply(events, parseSuppaEventE)
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

#' Parses splicing junctions
#'
#' NOTE: In case the -b V (Variable) option is selected (see below), some
#' variability is allowed in some of the boundaries. This is not accounted at
#' the moment.
#' 
#' @param event_type Chatacter: Type of the splicing event
#' @param junctions Character vector: Splicing junctions 
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
  parsed <- as.list(rep(NA, 8))
  names(parsed) <- c("C1 start", "C1 end",
                               "A1 start", "A1 end",
                               "A2 start", "A2 end",
                               "C2 start", "C2 end"
  )
  
  # Parse junction positions according to event type
  switch(event_type,
         "A3" = {
           parsed[["C1 end"]]   <- junctions[1]
           parsed[["C2 start"]] <- junctions[c(2,4)]
         },
         "A5" = {
           parsed[["C1 end"]]   <- junctions[c(3,1)]
           parsed[["C2 start"]] <- junctions[2]
         },
         "SE" =
           parsed[c("C1 end", "A1 start", "A1 end",
                              "C2 start")] <- junctions,
         "MX" =
           parsed[c("C1 end", "A1 start", "A1 end", "A2 start",
                              "A2 end", "C2 start")] <- junctions[-c(4,5)],
         "RI" =
           parsed[c("C1 start", "C1 end",
                              "C2 start", "C2 end")] <- junctions,
         "AF" =
           parsed[c("C1 start", "C1 end", "C2 end",
                              "A1 start", "A1 end")] <- junctions[1:5],
         "AL" = 
           parsed[c("C2 start", "C2 end", "C1 start",
                              "A1 start", "A1 end")] <- junctions[2:6]
  )
  return(parsed)
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