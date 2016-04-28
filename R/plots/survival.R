## TODO(NunoA): Select different time ranges (allow to change between days,
## weeks, months and years) or to even select custom

## TODO(NunoA): How to select groups? If there are intersections, merge them
## (maybe warn the user? maybe not needed as he'll see it?)

## TODO(NunoA): How to correctly do interval censoring?

# The name used for the plot must be unique
plot <- "Survival plots"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarPanel(
        selectizeInput(id("timeStart"), choices = NULL,
                       "Clinical attribute for starting time"),
        selectizeInput(id("timeStop"), choices = NULL,
                       "Clinical attribute for ending time (optional)"),
        helpText("In case there's no record for a sample, the days to last",
                 "follow up will be used instead."),
        selectizeInput(id("event"), choices = NULL, "Event of interest"),
        conditionalPanel(paste0("input.", id("timeStop"), "==''"),
                         radioButtons(id("censoring"), "Data censoring",
                                      selected="right",
                                      choices = c(Left="left", Right="right"))),
        actionButton(id("coxModel"), "Plot Cox Model"),
        actionButton(id("survivalCurves"), class="btn-primary",
                     "Plot survival curves")
    ),
    mainPanel(
        plotOutput(id(plot))
    )
)

#' Get the max number of days recorded of each sample for a given column
#'
#' @param col Character: column of interest
#' @param clinical Data.frame: clinical data
#'
#' @return Numeric vector with days recorded for given column (with NA when data
#' is not given)
sampleDays <- function(col, clinical) {
    cols <- grep(col, names(clinical))
    row <- apply(clinical[cols], 1, function(i)
        if(!all(is.na(i))) max(as.numeric(i), na.rm = TRUE) else NA)
    return(row)
}

server <- function(input, output, session) {
    # Update every time clinical data changes
    observe({
        clinical <- getClinicalData()
        if (!is.null(clinical)) {
            # Allow the user to select any "days_to" attribute available
            daysTo <- grep("days_to_", names(clinical), value=TRUE, fixed=TRUE)
            subDaysTo <- gsub(".*(days_to_.*)", "\\1", daysTo)
            choices <- unique(subDaysTo)
            names(choices) <- gsub("_", " ", choices, fixed=TRUE)
            names(choices) <- R.utils::capitalize(names(choices))
            updateSelectizeInput(session, id("timeStart"), choices=choices)
            updateSelectizeInput(
                session, id("timeStop"), choices = choices, options=list(
                    placeholder = "Select a column to use interval data",
                    onInitialize = I('function() { this.setValue(""); }')))
            names(choices) <- gsub("Days to ", "", names(choices), fixed=TRUE)
            names(choices) <- R.utils::capitalize(names(choices))
            updateSelectizeInput(session, id("event"), choices=choices)
        }
    })
    
    # Plot survival curve
    observeEvent(input[[id("survivalCurves")]], {
        output[[id(plot)]] <- renderPlot({
            isolate({
                # Get clinical data and column of interest
                clinical  <- getClinicalData()
                timeStart <- input[[id("timeStart")]]
                timeStop  <- input[[id("timeStop")]]
                dataEvent <- input[[id("event")]]
                censoring <- input[[id("censoring")]]
            })
            if (is.null(clinical)) {
                errorModal(session, "Clinical data missing",
                           "Insert clinical data first.")
            } else {
                rownames(clinical) <- toupper(
                    clinical$patient.bcr_patient_barcode)
                
                # Save the days from columns of interest in a data frame
                if (timeStop == "") timeStop <- NULL
                cols <- c(followup = "days_to_last_followup", start = timeStart,
                          stop = timeStop, event = dataEvent)
                colsDays <- lapply(cols, sampleDays, clinical)
                colsDays <- as.data.frame(colsDays)
                
                # Create new time using the starting time replacing the NAs with
                # days to last follow up
                nas <- is.na(colsDays$start)
                colsDays$time <- colsDays$start
                colsDays$time[nas] <- colsDays$followup[nas]
                
                # Indicate event of interest and groups
                colsDays$event <- ifelse(!is.na(colsDays$event), 1, 0)
                colsDays$groups <- clinical$patient.stage_event.pathologic_stage
                
                # Estimate and plot survival curves by groups
                if (is.null(timeStop))
                    form <- Surv(time, event, type = censoring) ~ groups
                else {
                    # Create new time using the ending time replacing the NAs
                    # with days to last follow up
                    nas <- is.na(colsDays$stop)
                    colsDays$time2 <- colsDays$stop
                    colsDays$time2[nas] <- colsDays$followup[nas]
                    
                    View(colsDays)
                    form <- Surv(time, time2, event, type = "interval") ~ groups
                }
                surv <- survfit(form, data = colsDays)
                plot(surv, lty=2:5, ylab="Proportion of individuals",
                     xlab="Time in days", col=1:4)
                legend("topright", sort(unique(colsDays$groups)), col=1:4,
                       lty=2:5)
            }
        })
    })
}