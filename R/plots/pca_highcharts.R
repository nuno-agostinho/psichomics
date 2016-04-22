## TODO(NunoA): project either individuals (as default) or events
## TODO(NunoA): add histogram in/above percentage of NAs per row to remove
## TODO(NunoA): add brushing capabilities (brush function from Shiny? only
## rectangle selection available?)
## TODO(NunoA): create clusters and use those clusters as groups of data

# The name used for the plot must be unique
plot <- "PCA highcharts"
id <- function(value) objectId(name, plot, value)

getMatchingRowNames <- function(selected, clinicalGroups, clinicalMatches) {
    # Get selected groups from clinical data
    rows <- clinicalGroups[selected, "Rows"]
    
    # Get names of the matching rows with the clinical data
    ns <- names(clinicalMatches[clinicalMatches %in% unlist(rows)])
    ns <- toupper(unique(ns))
    return(ns)
}

psiPCA <- function(psi, center = TRUE, scale. = FALSE, naTolerance = 30) {
    # Get individuals (rows) with less than a given percentage of NAs
    nas <- apply(psi, 1, function(row) sum(is.na(row)))
    psi <- psi[nas/ncol(psi)*100 <= naTolerance, ]
    if (nrow(psi) == 0) {
        print("No rows...")
        return(NULL)
    }
    
    # Replace NAs with the medians for each individual (row)
    medians <- miscTools::rowMedians(psi, na.rm = T)
    nas <- apply(psi, 1, function(row) sum(is.na(row)))
    psi[is.na(psi)] <- rep(medians, nas)
    
    # Perform principal component analysis on resulting data
    return(prcomp(psi, center = center, scale. = scale.))
}

ui <- list(
    sidebarPanel(
        checkboxGroupInput(id("preprocess"), "Preprocessing",
                           c("Center values" = "center",
                             "Scale values" = "scale"),
                           selected = c("center")),
        sliderInput(id("naTolerance"), "Percentage of NAs per row to tolerate",
                    min = 0, max=100, value=30, post="%"),
        fluidRow(
            column(8,
                   selectizeInput(id("dataGroups"), "Groups to perform PCA",
                                  choices = NULL, multiple = TRUE, 
                                  options = list(placeholder = "No groups created"))),
            column(2,
                   actionButton(id("dataGroups_selectAll"), "Select all", 
                                class = "inline_selectize"))),
        actionButton(id("editGroups"), "Edit groups"),
        actionButton(id("calculate"), class = "btn-primary", "Calculate PCA"),
        uiOutput(id("selectPC"))
    ), mainPanel(
        highchartOutput(id("scatterplot"))
    )
)

server <- function(input, output, session) {
    observeEvent(input[[id("editGroups")]], {
        env <- new.env()
        sys.source("R/data/2-groups.R", envir = env)
        env$server(input, output, session)
        
        showModal(session, "Groups", env$ui(),
                  iconName = "object-group", size = "small", style = "info",
                  printMessage = FALSE)
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
    
    # Performs principal component analysis (PCA)
    observeEvent(input[[id("calculate")]], {
        psi <- getInclusionLevels()
        
        # Subset data by the selected clinical groups
        selected <- isolate(input[[id("dataGroups")]])
        if (!is.null(selected)) {
            clinical <- getGroupsFrom("Clinical data")
            match <- getClinicalMatchFrom("Inclusion levels")
            ns <- getMatchingRowNames(selected, clinical, match)
            psi <- psi[ , ns]
        }
        
        # Raise error if data has no rows
        if (nrow(psi) == 0) stop("No data!")
        
        # Transpose the data to have individuals as rows
        psi <- t(psi)
        
        # Perform principal component analysis (PCA) on the subset data
        preprocess <- isolate(input[[id("preprocess")]])
        naTolerance <- isolate(input[[id("naTolerance")]])
        
        pca <- psiPCA(psi,
                      center = "center" %in% preprocess,
                      scale. = "scale" %in% preprocess,
                      naTolerance = naTolerance)
        sharedData$inclusionLevelsPCA <- pca
        
        # Save the original subset data used before performing PCA
        sharedData$inclusionLevelsPCA$original <- psi
    })
    
    # Interface and plots to help to select principal components
    output[[id("selectPC")]] <- renderUI({
        if (is.null(sharedData$inclusionLevelsPCA)) 
            return(NULL)
        
        sdevSq <- sharedData$inclusionLevelsPCA$sdev ^ 2
        perc <- sdevSq / sum(sdevSq) * 100
        names(perc) <- colnames(sharedData$inclusionLevelsPCA$x)
        
        label <- sprintf("%s (%s%% explained variance)", 
                         names(perc), round(perc, 2))
        choices <- setNames(names(perc), label)
        
        groups <- getGroupsFrom("Clinical data")
        
        list(
            hr(),
            fluidRow(
                column(6, selectizeInput(id("pcX"), "Choose X axis",
                                         choices = choices)),
                column(6, selectizeInput(id("pcY"), "Choose Y axis",
                                         choices = choices,
                                         selected = choices[[2]]))),
            fluidRow(
                column(8,
                       selectizeInput(
                           id("colorGroups"), "Groups to color the PCA",
                           multiple = TRUE, 
                           choices = getGroupsFrom("Clinical data")[, "Names"],
                           selected = getGroupsFrom("Clinical data")[, "Names"],
                           options = list(
                               placeholder = paste("Click in 'Select all' to",
                                                   "select all groups")))),
                column(2,
                       actionButton(id("colorGroups_selectAll"), "Select all", 
                                    class = "inline_selectize"))),
            actionButton(id("showVariancePlot"), "Show variance plot"),
            actionButton(id("plot"), class = "btn-primary", "Plot PCA")
        )
    })
    
    # Show variance plot
    observeEvent(input[[id("showVariancePlot")]], {
        showModal(session, "Variance plot", highchartOutput(id("variancePlot")), 
                  size = "large", style = "info", icon = "info-circle",
                  printMessage = FALSE)
    })
    
    # Select all color groups when pressing the respective "Select all" button
    observeEvent(input[[id("colorGroups_selectAll")]], {
        updateSelectizeInput(
            session, id("colorGroups"), 
            selected = getGroupsFrom("Clinical data")[, "Names"])
    })
    
    # Plots the explained variance plot
    output[[id("variancePlot")]] <- renderHighchart({
        if (is.null(sharedData$inclusionLevelsPCA)) 
            return(NULL)
        
        sdevSq <- sharedData$inclusionLevelsPCA$sdev ^ 2
        highchart() %>%
            hc_chart(zoomType = "xy", backgroundColor = NULL) %>%
            hc_title(text = paste("Explained variance by each",
                                  "Principal Component (PC)")) %>%
            hc_add_series(name = "PCs", data = sdevSq, type = "waterfall") %>%
            hc_plotOptions(series = list(dataLabels = list(
                align = "center",
                verticalAlign = "top",
                enabled = TRUE,
                formatter = JS(
                    "function() {",
                    "var total = ", sum(sdevSq), ";",
                    "var perc = (this.y/total) * 100;",
                    "return (Highcharts.numberFormat(this.y) +'<br/>'+",
                    "Highcharts.numberFormat(perc) + '%')}")))) %>%
            hc_xAxis(categories = colnames(sharedData$inclusionLevelsPCA$x),
                     crosshair = TRUE) %>%
            hc_yAxis(title = list(text = "Explained variance")) %>%
            hc_legend(enabled = FALSE) %>%
            hc_tooltip(pointFormat = '{point.name} {point.y:.5f}') %>%
            hc_exporting(enabled = TRUE,
                         buttons = list(contextButton = list(
                             text = "Export", y = -50,
                             verticalAlign = "bottom",
                             theme = list(fill = NULL))))
    })
    
    # Plots the principal component analysis
    observeEvent(input[[id("plot")]], {
        output[[id("scatterplot")]] <- renderHighchart({
            pca <- sharedData$inclusionLevelsPCA
            if (is.null(pca)) 
                return(NULL)
            
            isolate({
                xAxis <- input[[id("pcX")]]
                yAxis <- input[[id("pcY")]]
                selected <- input[[id("colorGroups")]]  
            })
            
            if (!is.null(xAxis) & !is.null(yAxis)) {
                imp <- as.data.frame(summary(pca)$importance)[2, ]
                perc <- as.numeric(imp)
                names(perc) <- names(imp)
                
                label <- sprintf("%s (%s%% explained variance)", 
                                 names(perc[c(xAxis, yAxis)]), 
                                 round(perc[c(xAxis, yAxis)], 2))
                
                hc <- highchart() %>%
                    hc_chart(zoomType = "xy") %>%
                    hc_xAxis(title = list(text = label[1])) %>%
                    hc_yAxis(title = list(text = label[2]))
                
                df <- data.frame(pca[["x"]])
                if (is.null(selected)) {
                    hc <- hc %>%
                        hc_scatter(df[[xAxis]], df[[yAxis]],
                                   sample = rownames(df)) %>%
                        hc_tooltip(pointFormat = "{point.sample}")
                } else {
                    # Subset data by the selected clinical groups
                    clinical <- getGroupsFrom("Clinical data")
                    match <- getClinicalMatchFrom("Inclusion levels")
                    
                    for (groupName in selected) {
                        ns <- getMatchingRowNames(groupName, clinical, match)
                        hc <- hc %>% 
                            hc_scatter(df[ns, xAxis], df[ns, yAxis],
                                       name = groupName,
                                       sample = rownames(df[ns, ]),
                                       showInLegend = TRUE) %>%
                            hc_tooltip(pointFormat = "{point.sample}")
                    }
                }
                return(hc)
            }
        })
    })
}