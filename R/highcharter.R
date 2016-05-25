#' Add scatter
#' 
#' @importFrom dplyr data_frame
#' @importFrom assertthat assert_that
#' @importFrom stats setNames
#' @importFrom highcharter %>% hc_add_series
#' 
#' @export
hc_scatter <- function (hc, x, y, z = NULL, color = NULL, label = NULL, 
                        showInLegend = FALSE, viridis.option = "D", ...) {
    assert_that(length(x) == length(y), 
                is.numeric(x), is.numeric(y))
    df <- data_frame(x, y)
    if (!is.null(z)) {
        assert_that(length(x) == length(z))
        df <- df %>% mutate(z = z)
    }
    if (!is.null(color)) {
        assert_that(length(x) == length(color))
        assert_that(viridis.option %in% c("A", "B", "C", "D"))
        cols <- colorize_vector(color, option = viridis.option)
        df <- df %>% mutate(valuecolor = color, color = cols)
    }
    if (!is.null(label)) {
        assert_that(length(x) == length(label))
        df <- df %>% mutate(label = label)
    }
    # Add arguments to data points if they match the length of the data
    args <- list(...)
    for (i in seq_along(args)) {
        if (length(x) == length(args[[i]]) && names(args[i]) != "name") {
            df <- cbind(df, setNames(list(args[i]), names(args[i])))
            args[[i]] <- character(0)
        }
    }
    
    ds <- list.parse3(df)
    type <- ifelse(!is.null(z), "bubble", "scatter")
    if (!is.null(label)) {
        dlopts <- list(enabled = TRUE, format = "{point.label}")
    }
    else {
        dlopts <- list(enabled = FALSE)
    }
    args <- Filter(length, args)
    do.call("hc_add_series", c(list(hc, data = ds, type = type, 
                                    showInLegend = showInLegend, 
                                    dataLabels = dlopts), args))
}

#' Plot survival curves using Highcharts
#' 
#' @param object A survfit object as returned from the \code{survfit} function
#' @param ... Extra parameters to pass to \code{hc_add_series} function
#' @param fun Name of function or function used to transform the survival curve:
#' \code{log} will put y axis on log scale, \code{event} plots cumulative events
#' (f(y) = 1-y), \code{cumhaz} plots the cumulative hazard function (f(y) =
#' -log(y)), and \code{cloglog} creates a complimentary log-log survival plot
#' (f(y) = log(-log(y)) along with log scale for the x-axis.
#' @param markTimes Label curves marked at each censoring time? TRUE by default
#' @param symbol Symbol to use as marker (plus sign by default)
#' @param markerColor Color of the marker ("black" by default); use NULL to use
#' the respective color of each series
#' @param ranges Plot interval ranges? FALSE by default
#' @param rangesOpacity Opacity of the interval ranges (0.3 by default)
#' 
#' @importFrom highcharter %>% hc_add_series highchart hc_tooltip hc_yAxis
#' hc_plotOptions
#' @return Highcharts object to plot survival curves
#' 
#' @examples
#' 
#' # Plot Kaplan-Meier curves
#' require("survival")
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
#' 
#' @export
hchart.survfit <- function(object, ..., fun = NULL, markTimes = TRUE,
                           symbol = fa_icon_mark("plus"), markerColor = "black",
                           ranges = FALSE, rangesOpacity = 0.3) {
    
    group <- NULL
    
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
        hc_tooltip(shared = TRUE) %>%
        hc_yAxis(min=ymin, max=ymax) %>%
        hc_plotOptions(line = list(marker = list(enabled = FALSE)))
    
    count <- 0
    
    # Process groups by columns (CoxPH-like) or in a single column
    if(!is.null(ncol(object$surv))) {
        groups <- seq(ncol(object$surv))
    } else {
        groups <- names(strata)
    }
    
    for (name in groups) {
        if (!is.null(ncol(object$surv))) {
            df <- df[c("x", paste(c("y", "low", "up"), col, sep="."))]
            names(df) <- c("x", "y", "low", "up")
            submark <- mark
        } else {
            df <- subset(data, group == name)
            submark <- mark[data$group == name]
        }
        
        # Add first value if there is no value for time at 0 in the data
        if (!0 %in% df$x)
            first <- list(list(x=0, y=firsty))
        else
            first <- NULL
        
        # Mark events
        ls <- list.parse3(df)
        if (markTimes)
            ls[submark] <- lapply(ls[submark], c, marker=marker)
        
        hc <- hc %>% hc_add_series(
            data=c(first, ls), step="left", name=name, zIndex=1,
            color=JS("Highcharts.getOptions().colors[", count, "]"),
            ...)
        
        if (ranges && !is.null(object$upper)) {
            # Add interval range
            range <- lapply(ls, function(i) 
                setNames(i[c("x", "low", "up")], NULL))
            hc <- hc %>% hc_add_series(
                data=range, step="left", name="Ranges", type="arearange",
                zIndex=0, linkedTo=':previous', fillOpacity=rangesOpacity, 
                lineWidth=0,
                color=JS("Highcharts.getOptions().colors[", count, "]"),
                ...)
        }
        count <- count + 1
    }
    
    return(hc)
}


#' Shorcut to create a density plot 
#' 
#' @param hc A \code{highchart} \code{htmlwidget} object. 
#' @param x A numeric vector
#' @param area A boolean value to show or not the area
#' @param ... Aditional shared arguments for the data series
#'   (\url{http://api.highcharts.com/highcharts#series}).
#' 
#' @importFrom stats density
#' @importFrom highcharter %>% hc_add_series list.parse3
#' @examples
#' 
#' require(`highchart`)
#' highchart() %>%
#'   hc_add_series_density(rnorm(1000)) %>%
#'   hc_add_series_density(rexp(1000), area = TRUE)
#'  
#' @export
hc_add_series_density <- function (hc, x, area = FALSE, ...) {
    if(is.numeric(x)) x <- density(x)
    type <- ifelse(area, "areaspline", "spline")
    data <- list.parse3(data.frame(cbind(x = x$x, y = x$y)))
    return(hc %>% hc_add_series(data = data, type = type, ...))
}