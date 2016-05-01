## TODO(NunoA): Select different time ranges (allow to change between days,
## weeks, months and years) or to even select custom

## TODO(NunoA): Should groups be merged if there's an intersection?

## TODO(NunoA): How to correctly do interval censoring?

## TODO(NunoA): Icon symbol used for cross are moved when hiding/showing series
## and besides they are never really centered...

# The name used for the plot must be unique
plot <- "Survival plots high"
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
            updateSelectizeInput(session, id("event"), choices=choices,
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
                if (timeStop == "") timeStop <- NULL
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

plotSurvCurves <- function(surv) {
    # The X axis will be the time and the Y axis the survival probability
    # These variables don't separate by groups, so we'll have to split them
    x <- surv$time
    y <- surv$surv
    
    # Check which points should be marked
    mark <- ifelse(surv$n.censor, 1, 0)
    
    # Check if there are groups
    if (is.null(surv$strata)) {
        group <- c("Test" = length(surv$time))
    } else {
        group <- surv$strata
        names(group) <- gsub(".*=", "", names(group))
    }
    
    marker <- list(list(fillColor="black", symbol=fa_icon_mark("plus"),
                        enabled=TRUE))
    dont <- list(list(enabled=FALSE))
    data <- list.parse3(data.frame(x, y, mark, group=rep(names(group), group), 
                                   stringsAsFactors = FALSE))
    data <- lapply(data, function(i) c(i, marker=ifelse(i$mark, marker, dont)))
    
    hc <- highchart() %>%
        hc_chart(zoomType="xy") %>%
        hc_yAxis(min=0, max=1, title=list(text="Proportion of individuals")) %>%
        hc_xAxis(title=list(text="Time in days"))
    
    for (name in names(group)) {
        ls <- lapply(data, function(i) if (i$group == name) i)
        ls <- Filter(Negate(is.null), ls)
        
        first <- NULL
        if (!0 %in% vapply(ls, "[[", "x", FUN.VALUE = numeric(1))) {
            # Add first value if there is no x = 0 in the data
            first <- list(list(x=0, y=1, marker=dont))
        }
        
        hc <- hc %>% hc_add_series(data=c(first, ls), step="left", name=name)
    }
    return(hc)
}