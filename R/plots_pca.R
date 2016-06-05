## TODO(NunoA): project either individuals (as default) or events
## TODO(NunoA): add histogram in/above percentage of NAs per row to remove
## TODO(NunoA): add brushing capabilities (brush function from Shiny? only
## rectangle selection available?)
## TODO(NunoA): logarithmic values
## TODO(NunoA): BoxCox transformation

## TODO(NunoA): create clusters and use those clusters as groups of data
##
##  km <- kmeans(scores, centers = 7, nstart = 5)
##  ggdata <- data.frame(scores, Cluster=km$cluster, Species=df$Species)
##  ggplot(ggdata) +
##      geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
##      stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
##          geom="polygon", level=0.95, alpha=0.2) +
##      guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster"))
## 
## Source:
## http://stackoverflow.com/questions/20260434/test-significance-of-clusters-on-a-pca-plot

# The name used for the plot must be unique
plot <- "Principal component analysis"
id <- function(value) objectId(name, plot, value)

getMatchingRowNames <- function(selected, clinicalGroups, clinicalMatches) {
    # Get selected groups from clinical data
    rows <- clinicalGroups[selected, "Rows"]
    
    # Get names of the matching rows with the clinical data
    ns <- names(clinicalMatches[clinicalMatches %in% unlist(rows)])
    ns <- toupper(unique(ns))
    return(ns)
}

#' @importFrom stats prcomp
#' @importFrom miscTools rowMedians
psiPCA <- function(psi, center = TRUE, scale. = FALSE, naTolerance = 30) {
    # Get individuals (rows) with less than a given percentage of NAs
    nas <- apply(psi, 1, function(row) sum(is.na(row)))
    # hist(nas/ncol(psi)*100)
    psi <- psi[nas/ncol(psi)*100 <= naTolerance, , drop = FALSE]
    if (nrow(psi) == 0) return(NULL)
    
    # Replace NAs with the medians for each individual (row)
    medians <- rowMedians(psi, na.rm=TRUE)
    nas <- apply(psi, 1, function(row) sum(is.na(row)))
    psi[is.na(psi)] <- rep(medians, nas)
    
    # Perform principal component analysis on resulting data
    pca <- prcomp(psi, center = center, scale. = scale.)
    # PCA is useless if done with only one point
    if (nrow(pca$x) == 1) return(NULL)
    return(pca)
}

#' @importFrom highcharter highchartOutput
ui <- function() {
    list(
        sidebarPanel(
            checkboxGroupInput(id("preprocess"), "Preprocessing",
                               c("Center values" = "center",
                                 "Scale values" = "scale"),
                               selected = c("center")),
            sliderInput(id("naTolerance"), "Percentage of NAs per individual to tolerate",
                        min = 0, max=100, value=30, post="%"),
            fluidRow(
                column(9, selectizeInput(id("dataGroups"),
                                         "Clinical groups to perform PCA",
                                         choices = NULL, multiple = TRUE)),
                column(2, actionButton(id("dataGroups_selectAll"), "Select all",
                                       class="inline_selectize"))),
            actionButton(id("editGroups"), "Edit groups"),
            actionButton(id("calculate"), class = "btn-primary", "Calculate PCA"),
            uiOutput(id("selectPC"))
        ), mainPanel(
            highchartOutput(id("scatterplot"))
        )
    )
}

server <- function(input, output, session) {
    observeEvent(input[[id("editGroups")]], {
        env <- new.env()
        sys.source("R/data/2-groups.R", envir = env)
        env$server(input, output, session)
        
        showModal(session, "Groups", env$ui(),  iconName = "object-group",
                  style = "info")
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
        psi <- isolate(getInclusionLevels())
        
        if (is.null(psi)) {
            errorModal(session, "Inclusion levels missing",
                       "Insert or calculate exon/intron inclusion levels.")
        } else {
            # Subset data by the selected clinical groups
            selected <- isolate(input[[id("dataGroups")]])
            if (!is.null(selected)) {
                clinical <- isolate(getGroupsFrom("Clinical data"))
                match <- getClinicalMatchFrom("Inclusion levels")
                ns <- getMatchingRowNames(selected, clinical, match)
                psi <- psi[ , ns]
            }
            
            # Raise error if data has no rows
            if (nrow(psi) == 0)
                errorModal(session, "No data!", paste(
                    "Calculation returned nothing. Check if everything is as",
                    "expected and try again."))
            
            # Transpose the data to have individuals as rows
            psi <- t(psi)
            
            # Perform principal component analysis (PCA) on the subset data
            isolate({
                preprocess <- input[[id("preprocess")]]
                naTolerance <- input[[id("naTolerance")]]
            })
            
            pca <- psiPCA(psi, naTolerance = naTolerance,
                          center = "center" %in% preprocess,
                          scale. = "scale" %in% preprocess)
            if (is.null(pca)) {
                errorModal(session, "No individuals to plot PCA", 
                           "Try increasing the tolerance of NAs per individual.")
            }
            sharedData$inclusionLevelsPCA <- pca
        }
    })
    
    # Show variance plot
    observeEvent(input[[id("showVariancePlot")]], {
        infoModal(session, "Variance plot", highchartOutput(id("variancePlot")), 
                  size = "large")
    })
    
    # Select all color groups when pressing the respective "Select all" button
    observeEvent(input[[id("colorGroups_selectAll")]], {
        updateSelectizeInput(
            session, id("colorGroups"), 
            selected = getGroupsFrom("Clinical data")[, "Names"])
    })
    
    observe({
        pca <- sharedData$inclusionLevelsPCA
        if (is.null(pca)) return(NULL)
        
        imp <- summary(pca)$importance[2, ]
        perc <- as.numeric(imp)
        names(perc) <- names(imp)
        
        # Interface and plots to help to select principal components
        output[[id("selectPC")]] <- renderUI({
            label <- sprintf("%s (%s%% explained variance)", 
                             names(perc), round(perc * 100, 2))
            choices <- setNames(names(perc), label)
            groups <- getGroupsFrom("Clinical data")
            
            tagList(
                hr(),
                selectizeInput(id("pcX"), "Choose X axis", choices = choices),
                selectizeInput(id("pcY"), "Choose Y axis", choices = choices,
                               selected = choices[[2]]),
                fluidRow(
                    column(9,
                           selectizeInput(
                               id("colorGroups"), 
                               "Clinical groups to color the PCA",
                               choices = groups[, "Names"], multiple = TRUE,
                               options = list(placeholder = ifelse(
                                   length(groups) > 0,
                                   "Click 'Select all'to select all groups",
                                   "No groups created")))),
                    column(2,
                           actionButton(id("colorGroups_selectAll"),
                                        "Select all", 
                                        class = "inline_selectize"))),
                checkboxGroupInput(id("plotShow"), "Show in plot",
                                   c("Individuals", "Loadings"),
                                   selected = c("Individuals")),
                actionButton(id("showVariancePlot"), "Show variance plot"),
                actionButton(id("plot"), class = "btn-primary", "Plot PCA")
            )
        })
        
        # Plots the explained variance plot
        output[[id("variancePlot")]] <- renderHighchart({
            sdevSq <- pca$sdev ^ 2
            
            highchart() %>%
                hc_chart(zoomType = "xy", backgroundColor = NULL) %>%
                hc_title(text = paste("Explained variance by each",
                                      "Principal Component (PC)")) %>%
                hc_add_series(name = "PCs", data = sdevSq,
                              type = "waterfall") %>%
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
                hc_xAxis(categories = colnames(pca[["x"]]), 
                         crosshair = TRUE) %>%
                hc_yAxis(title = list(text = "Explained variance")) %>%
                hc_legend(enabled = FALSE) %>%
                hc_tooltip(pointFormat = 
                               '{point.name} {point.y:.2f} {point.perc}') %>%
                hc_exporting(enabled = TRUE,
                             buttons = list(contextButton = list(
                                 text = "Export", y = -50,
                                 verticalAlign = "bottom",
                                 theme = list(fill = NULL))))
        })
        
        # Plots the principal component analysis
        observeEvent(input[[id("plot")]], {
            output[[id("scatterplot")]] <- renderHighchart({
                isolate({
                    xAxis <- input[[id("pcX")]]
                    yAxis <- input[[id("pcY")]]
                    selected <- input[[id("colorGroups")]]
                    show <- input[[id("plotShow")]]
                })
                
                if (!is.null(xAxis) & !is.null(yAxis)) {
                    label <- sprintf("%s (%s%% explained variance)", 
                                     names(perc[c(xAxis, yAxis)]), 
                                     round(perc[c(xAxis, yAxis)]*100, 2))
                    
                    hc <- highchart() %>%
                        hc_chart(zoomType = "xy") %>%
                        hc_xAxis(title = list(text = label[1])) %>%
                        hc_yAxis(title = list(text = label[2])) %>%
                        hc_tooltip(pointFormat = "{point.sample}")
                    
                    if ("Individuals" %in% show) {
                        df <- data.frame(pca$x)
                        if (is.null(selected)) {
                            hc <- hc %>% hc_scatter(df[[xAxis]], df[[yAxis]],
                                                    sample = rownames(df))
                        } else {
                            # Subset data by the selected clinical groups
                            clinical <- getGroupsFrom("Clinical data")
                            match <- getClinicalMatchFrom("Inclusion levels")
                            
                            for (groupName in selected) {
                                ns <- getMatchingRowNames(groupName, clinical,
                                                          match)
                                ns <- ns[ns %in% rownames(df)]
                                hc <- hc %>%
                                    hc_scatter(df[ns, xAxis],
                                               df[ns, yAxis],
                                               name = groupName,
                                               sample = rownames(df[ns, ]),
                                               showInLegend = TRUE)
                            }
                        }
                    }
                    if ("Loadings" %in% show) {
                        m <- data.frame(pca$rotation)
                        # For loadings, add series (don't add to legend)
                        hc <- hc %>% hc_scatter(m[[xAxis]], m[[yAxis]])
                    }
                    return(hc)
                }
            })
        })
    })
}