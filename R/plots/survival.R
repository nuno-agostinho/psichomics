## TODO(NunoA): Select different time ranges (allow to change between days,
## weeks, months and years) or to even select custom

## TODO(NunoA): Should groups be merged if there's an intersection?

## TODO(NunoA): How to correctly do interval censoring?

## TODO(NunoA): Icon symbol used for cross in survival curves are moved when
## hiding/showing series and they're not centered...

# The name used for the plot must be unique
plot <- "Survival curves"
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
        radioButtons(id("modelTerms"), selected="groups", inline=TRUE,
                     "Select model terms of the right-hand using",
                     choices=c("Clinical groups"="groups",
                               "Formula"="formula")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", id("modelTerms"), "groups"),
            fluidRow(
                column(9, selectizeInput(id("dataGroups"), 
                                         "Clinical groups to use",
                                         choices = NULL, multiple = TRUE)),
                column(2, actionButton(id("dataGroups_selectAll"), "Select all",
                                       class="inline_selectize"))),
            checkboxInput(id("showOutGroup"), "Show data outside chosen groups",
                          value = FALSE)),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", id("modelTerms"), "formula"),
            textAreaInput(id("formula"), "Formula for right-hand side"),
            uiOutput(id("formulaAutocomplete"))),
        checkboxInput(id("ranges"), "Show interval ranges", value = FALSE),
        actionButton(id("coxModel"), "Fit Cox PH model"),
        actionButton(id("survivalCurves"), class="btn-primary", 
                     "Plot survival curves")
    ),
    mainPanel(
        highchartOutput(id(plot)),
        uiOutput(id("coxphUI"))
    )
)

#' Process survival data to calculate survival curves
#' 
#' @param timeStart Integer: starting time of the interval or follow up time
#' @param timeStop Integer: ending time of the interval
#' @param timeEvent Integer: time of the event of interest
#' @param clinical Data.frame: clinical data
#' 
#' @details The event time will only be used to determine whether the event has
#' happened (1) or not in case of NAs (0)
#' 
#' @return Data frame with terms needed to calculate survival curves
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

#' Process survival curves terms to calculate survival curves
#'
#' @param session Session object from Shiny function
#' @param group Character: group of each individual 
#' @param clinical Data frame: clinical data
#' @param outGroup Boolean: show group with outsiders (FALSE by default)
#' @param censoring Character: censor using "left", "right", "interval" or
#' "interval2"
#' @param timeStart Integer: staring time
#' @param timeStop Integer: ending time (needed only for interval-censored data)
#' @param dataEvent Character: event
#' @param modelTerms Character: use "groups" or "formula" for the survival 
#' curves?
#' @param formulaStr Character: formula to use
#' @param coxph Boolean: fit a Cox proportional hazards regression model? FALSE 
#' by default
#'
#' @return A list with a \code{formula} object and a data frame with terms
#' needed to calculate survival curves
processSurvTerms <- function(session, group, clinical, outGroup, censoring, 
                             timeStart, timeStop, dataEvent, modelTerms, 
                             formulaStr, coxph=FALSE) {
    # Save the days from columns of interest in a data frame
    fillGroups <- groupPerPatient(group, nrow(clinical), outGroup)
    
    # Ignore timeStop if interval-censoring is not selected
    if (!grepl("interval", censoring, fixed=TRUE) || timeStop == "") 
        timeStop <- NULL
    
    # Check if using or not interval-censored data
    formulaSurv <- ifelse(is.null(timeStop),
                          "Surv(time, event, type=censoring) ~", 
                          "Surv(time, time2, event, type=censoring) ~")
    
    survTime <- processSurvData(timeStart, timeStop, dataEvent, fillGroups, 
                                clinical)
    
    # Estimate survival curves by groups or using formula
    if (modelTerms == "groups") {
        formulaTerms <- "groups"
    } else if (modelTerms == "formula") {
        formulaTerms <- formulaStr
        if (formulaTerms == "") {
            errorModal(session, "Error in formula",
                       "The formula field can't be empty.")
            return(NULL)
        }
        survTime <- cbind(survTime, clinical)
    }
    
    form <- tryCatch(formula(paste(formulaSurv, formulaTerms)), error = return)
    if ("simpleError" %in% class(form)) {
        errorModal(session, "Formula error",
                   "Maybe you misplaced a ", tags$kbd("+"), ", ", tags$kbd(":"), 
                   " or ", tags$kbd("*"), "?", br(), br(),  
                   "The following error was raised:", br(), 
                   tags$code(form$message))
        return(NULL)
    }
    
    if (coxph)
        res <- tryCatch(coxph(form, data=survTime), error=return)
    else
        res <- list(form=form, survTime=survTime)
    return(res)
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
    # Update available clinical data attributes to use in a formula
    output[[id("formulaAutocomplete")]] <- renderUI({
        attributes <- names(getClinicalData())
        textComplete(id("formula"), attributes)
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
                clinical   <- getClinicalData()
                timeStart  <- input[[id("timeStart")]]
                timeStop   <- input[[id("timeStop")]]
                dataEvent  <- input[[id("event")]]
                censoring  <- input[[id("censoring")]]
                outGroup   <- input[[id("showOutGroup")]]
                modelTerms <- input[[id("modelTerms")]]
                formulaStr <- input[[id("formula")]]
                intRanges  <- input[[id("ranges")]]
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
                # Calculate survival curves
                survTerms <- processSurvTerms(session, dataGroups, clinical, 
                                              outGroup, censoring, timeStart, 
                                              timeStop, dataEvent, modelTerms, 
                                              formulaStr)
                form <- survTerms$form
                data <- survTerms$survTime
                surv <- tryCatch(survfit(form, data = data),
                                 error = return)
                
                if ("simpleError" %in% class(surv)) {
                    errorModal(session, "Formula error",
                               "The following error was raised:", br(),
                               tags$code(surv$message))
                    return(NULL)
                }
                
                # If there's an error with survdiff, show p-value of 0
                pvalue <- tryCatch({
                    # Test the difference between survival curves
                    diff <- survdiff(form, data = data)
                    
                    # Calculate p-value with 5 significant numbers
                    pvalue <- 1 - stats::pchisq(diff$chisq, length(diff$n) - 1)
                    signif(pvalue, 5)
                }, error = function(e) 0)
                
                # Plot survival curves
                hc <- hchart(surv, ranges = intRanges) %>%
                    hc_chart(zoomType="xy") %>%
                    hc_yAxis(title=list(text="Proportion of individuals")) %>%
                    hc_xAxis(title=list(text="Time in days")) %>% 
                    hc_tooltip(headerFormat='Time: {point.x}<br>') %>%
                    hc_credits(href=NA, enabled=TRUE,
                               text=paste("p-value:", pvalue)) %>%
                    hc_tooltip(crosshairs=TRUE)
            }
        })
    })
    
    # Plot cox model
    observeEvent(input[[id("coxModel")]], {
        isolate({
            # Get user input
            clinical   <- getClinicalData()
            timeStart  <- input[[id("timeStart")]]
            timeStop   <- input[[id("timeStop")]]
            dataEvent  <- input[[id("event")]]
            censoring  <- input[[id("censoring")]]
            outGroup   <- input[[id("showOutGroup")]]
            modelTerms <- input[[id("modelTerms")]]
            formulaStr <- input[[id("formula")]]
            intRanges  <- input[[id("ranges")]]
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
            # Calculate survival curves
            survTerms <- processSurvTerms(session, dataGroups, clinical, 
                                          outGroup, censoring, timeStart, 
                                          timeStop, dataEvent, modelTerms, 
                                          formulaStr, cox=TRUE)
            if ("simpleError" %in% class(survTerms)) {
                errorModal(session, "Formula error",
                           "The following error was raised:", br(),
                           tags$code(survTerms$message))
                return(NULL)
            }
            
            surv <- survfit(survTerms)
            summary <- summary(survTerms)
            print(summary)
            
            output[[id("coxphUI")]] <- renderUI({
                # highchartOutput(id("coxPlot"))
                list(
                    dataTableOutput(id("coxGroups")),
                    dataTableOutput(id("coxTests"))
                )
            })
            
            output[[id("coxGroups")]] <- renderDataTable({
                groups <- cbind(rownames(summary$coefficients),
                                signif(summary$coefficients, 5),
                                signif(summary$conf.int[ , 2:4], 5))
                return(groups)
            }, options = list(scrollX = TRUE))
            
            output[[id("coxTests")]] <- renderDataTable({
                tests <- rbind("Wald test"=summary$waldtest, 
                               "Log test"=summary$logtest,
                               "Score (logrank) test"=summary$sctest)
                tests <- cbind(rownames(tests), tests)
                colnames(tests) <- c("Statistical test", "Value", 
                                     "Degrees of freedom", "p-value")
                return(tests)
            }, options = list(scrollX = TRUE))
            
            # output[[id("coxPlot")]] <- renderHighchart({
            #     # Plot survival curves
            #     hchart(surv, ranges = intRanges) %>%
            #         hc_chart(zoomType="xy") %>%
            #         hc_yAxis(title=list(text="Proportion of individuals")) %>%
            #         hc_xAxis(title=list(text="Time in days"))
            # })
        }
    })
}