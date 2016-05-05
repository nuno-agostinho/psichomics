## TODO(NunoA): Select different time ranges (allow to change between days,
## weeks, months and years) or to even select custom

## TODO(NunoA): Should groups be merged if there's an intersection?

## TODO(NunoA): How to correctly do interval censoring?

## TODO(NunoA): Icon symbol used for cross are moved when hiding/showing series
## and besides they are never really centered...

# The name used for the plot must be unique
plot <- "Survival plots Highcharts"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarPanel(
        radioButtons(id("censoring"), "Data censoring", selected="right",
                     inline=TRUE, choices=c(Left="left",
                                            Right="right",
                                            Interval="interval",
                                            "Interval 2" = "interval2")),
        selectizeInput(id("timeStart"), choices = NULL, "Follow up time"),
        # If the chosen censoring contains the word 'interval', show this input
        conditionalPanel(
            paste0("input.", id("censoring"), ".indexOf('interval') > -1"),
            selectizeInput(id("timeStop"), choices = NULL, "Ending time")),
        helpText("In case there's no record for a patient, the days to last",
                 "follow up will be used instead."),
        selectizeInput(id("event"), choices = NULL, "Event of interest"),
        fluidRow(
            column(9, selectizeInput(id("dataGroups"), "Clinical groups to use",
                                     choices = NULL, multiple = TRUE)),
            column(2, actionButton(id("dataGroups_selectAll"), "Select all",
                                   class="inline_selectize"))),
        textAreaInput(id("formula"), "Insert formula"),
        uiOutput(id("formulaAutocomplete")),
        checkboxInput(id("showOutGroup"), "Show data outside chosen groups",
                      value = FALSE),
        checkboxInput(id("coxModel"), "Use Cox proportional hazards model"),
        actionButton(id("survivalCurves"), class="btn-primary", 
                     "Plot survival curves")
    ),
    mainPanel(
        highchartOutput(id(plot))
    )
)

#' Process survival data
#' 
#' @param time Integer: starting time of the interval or follow up time
#' @param time2 Integer: ending time of the interval
#' @param timeEvent Integer: time of the event of interest
#' @param clinical Data.frame: clinical data
#' 
#' @details The event time will only be used to determine whether the event has
#' happened (1) or not in case of NAs (0)
#' 
#' @return Data frame with 
processSurvData <- function(timeStart, timeStop, event, groups, clinical) {
    cols <- c(followup = "days_to_last_followup", start = timeStart,
              stop = timeStop, event = event)
    survTime <- lapply(cols, timePerPatient, clinical)
    survTime <- as.data.frame(survTime)
    
    # Create new time using the starting time replacing the NAs with
    # days to last follow up
    nas <- is.na(survTime$start)
    survTime$time <- survTime$start
    survTime$time[nas] <- survTime$followup[nas]
    
    # Indicate event of interest and groups
    survTime$event <- ifelse(!is.na(survTime$event), 1, 0)
    survTime$groups <- groups
    
    if (!is.null(timeStop)) {
        # Create new time using the ending time replacing the NAs
        # with days to last follow up
        nas <- is.na(survTime$stop)
        survTime$time2 <- survTime$stop
        survTime$time2[nas] <- survTime$followup[nas]
    }
    return(survTime)
}

#' Get all columns matching a given string and return a single vector with the
#' max time for each patient if available
#'
#' @param col Character: column of interest
#' @param clinical Data.frame: clinical data
#'
#' @return Numeric vector with days recorded for columns of interest
timePerPatient <- function(col, clinical) {
    cols <- grep(col, names(clinical))
    row <- apply(clinical[cols], 1, function(i)
        if(!all(is.na(i))) max(as.numeric(i), na.rm = TRUE) else NA)
    return(row)
}

server <- function(input, output, session) {
    output[[id("formulaAutocomplete")]] <- renderUI({
        words <- names(getClinicalData())
        textComplete(id("formula"), words)
    })
    
    # Update available group choices to select
    observe({
        groups <- getGroupsFrom("Clinical data")
        updateSelectizeInput(
            session, id("dataGroups"), choices = groups[, "Names"],
            options = list(placeholder =
                               ifelse(length(groups) > 0,
                                      "Click 'Select all' to select all groups",
                                      "No groups created")))
    })
    
    # Select all data groups when pressing the respective "Select all" button
    observeEvent(input[[id("dataGroups_selectAll")]], {
        updateSelectizeInput(
            session, id("dataGroups"), 
            selected = getGroupsFrom("Clinical data")[, "Names"])
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        label <- "Follow up time"
        if (grepl("interval", input[[id("censoring")]], fixed=TRUE))
            label <- "Starting time"
        updateSelectizeInput(session, id("timeStart"), label=label)
    })
    
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
            updateSelectizeInput(session, id("timeStart"), choices=choices,
                                 selected="days_to_death")
            updateSelectizeInput(
                session, id("timeStop"), choices = choices, options=list(
                    onInitialize = I('function() { this.setValue(""); }')))
            names(choices) <- gsub("Days to ", "", names(choices), fixed=TRUE)
            names(choices) <- R.utils::capitalize(names(choices))
            updateSelectizeInput(session, id("event"), 
                                 choices=list(
                                     "Suggested events"=choices,
                                     "All clinical data columns"=names(clinical)),
                                 selected="days_to_death")
        }
    })
    
    # Plot survival curve
    observeEvent(input[[id("survivalCurves")]], {
        output[[id(plot)]] <- renderHighchart({
            isolate({
                # Get user input
                clinical  <- getClinicalData()
                timeStart <- input[[id("timeStart")]]
                timeStop  <- input[[id("timeStop")]]
                dataEvent <- input[[id("event")]]
                censoring <- input[[id("censoring")]]
                showOther <- input[[id("showOutGroup")]]
                coxModel  <- input[[id("coxModel")]]
                
                # Get chosen groups
                chosen <- input[[id("dataGroups")]]
                dataGroups <- getGroupsFrom("Clinical data")[chosen, , drop=F]
            })
            
            if (is.null(clinical)) {
                errorModal(session, "Clinical data missing",
                           "Insert clinical data first.")
            } else if (nrow(dataGroups) > 0 && 
                       anyDuplicated(unlist(dataGroups[, "Rows"])) > 0) {
                # If the chosen groups have any intersections
                errorModal(session, "Clinical groups intercept",
                           "There is an interception between clinical groups.")
            } else {
                # Save the days from columns of interest in a data frame
                fillGroups <- groupPerPatient(dataGroups, nrow(clinical), 
                                              showOther)
                
                # Ignore timeStop if interval-censoring is not selected
                if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "") 
                    timeStop <- NULL
                
                survTime <- processSurvData(timeStart, timeStop, dataEvent,
                                            fillGroups, clinical)
                
                # Estimate survival curves by groups
                if (is.null(timeStop))
                    form <- Surv(time, event, type=censoring) ~ groups
                else
                    form <- Surv(time, time2, event, type=censoring) ~ groups
                
                if (coxModel) {
                    fit <- coxph(form, data = survTime)
                    surv <- survfit(fit)
                } else {
                    surv <- survfit(form, data = survTime)
                }
                
                # Plot survival curves
                groups <- sort(unique(survTime$groups))
                nos <- seq_along(groups)
                plotSurvCurves(surv)
            }
        })
    })
}

#' Plot survival curves using Highcharts
#' 
#' @param surv survfit object: survival curves returned from the \code{survfit}
#' function
#' @param fun Character (name of function) or function: used to transform the
#' survival curve. \code{log} will put y axis on log scale, \code{event} plots 
#' cumulative events (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard 
#' function (f(y) = -log(y)), and \code{cloglog} creates a complimentary log-log
#' survival plot (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param ymin Integer: minimum Y to plot if all Y values are higher than the 
#' indicated (0 by default)
#' @param ymax Integer: maximum Y to plot if all Y values are lower than the 
#' indicated (1 by default)
#' @param markTimes Boolean: should times be marked? TRUE by default
#' @param markerSymbol Character: symbol to use as marker (plus sign by default)
#' @param markerColor Character: color of the marker ("black" by default); use
#' NULL to color using the series color
#' @param ranges Boolean: plot interval ranges? FALSE by default
#' 
#' @return Highchart object to plot survival curves
plotSurvCurves <- function(surv, fun = NULL, ymin=0, ymax=1, markTimes=TRUE,
                           markerSymbol=fa_icon_mark("plus"),
                           markerColor="black", ranges=FALSE) {
    # Check if there are groups
    if (is.null(surv$strata))
        group <- c("Series 1" = length(surv$time))
    else
        group <- surv$strata
    
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
    } else 
        tfun <- function(x) x
    
    firsty <- tfun(1)
    surv$surv <- tfun(surv$surv)
    if (ranges && !is.null(surv$upper)) {
        surv$upper <- tfun(surv$upper)
        surv$lower <- tfun(surv$lower)
    }
    
    # Data markers
    noMarker <- list(list(enabled=FALSE))
    marker <- ifelse(
        markTimes,
        list(list(fillColor=markerColor, symbol=markerSymbol, enabled=TRUE)),
        noMarker)
    
    # Prepare data
    mark <- ifelse(surv$n.censor == 1, 1, 0)
    data <- data.frame(x=surv$time, y=surv$surv, mark,
                       up=surv$upper, low=surv$lower,
                       group=rep(names(group), group), 
                       stringsAsFactors = FALSE)
    
    # Adjust Y axis range
    yValues <- data$y
    ymin <- ifelse(min(yValues) >= ymin, ymin, min(yValues))
    ymax <- ifelse(max(yValues) <= ymax, ymax, max(yValues))
    
    hc <- highchart() %>%
        hc_chart(zoomType="xy") %>%
        hc_tooltip(shared = TRUE) %>%
        hc_yAxis(min=ymin, max=ymax, 
                 title=list(text="Proportion of individuals")) %>%
        hc_xAxis(title=list(text="Time in days"))
    
    count <- 0
    for (name in names(group)) {
        df <- subset(data, group == name)
        
        # Add first value if there is no x=0 in the data
        first <- NULL
        if (!0 %in% df$x)
            first <- list(list(x=0, y=firsty, marker=noMarker))
        
        # Mark events
        ls <- list.parse3(df)
        if (markTimes) {
            ls <- lapply(ls, function(i)
                c(i, marker=ifelse(i$mark, marker, noMarker)))
        }
        
        hc <- hc %>% hc_add_series(
            data=c(first, ls), step="left", name=name, zIndex=1,
            color=JS("Highcharts.getOptions().colors[", count, "]"))
        
        if (ranges) {
            # Add interval range
            range <- lapply(ls, function(i) 
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=0.3, lineWidth=0,
                color=JS("Highcharts.getOptions().colors[", count, "]"))
        }
        count <- count + 1
    }
    
    return(hc)
}