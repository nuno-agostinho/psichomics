# Plot p-values of survival difference between groups based on multiple cutoffs

Plot p-values of survival difference between groups based on multiple
cutoffs

## Usage

``` r
plotSurvivalPvaluesByCutoff(
  clinical,
  data,
  censoring,
  event,
  timeStart,
  timeStop = NULL,
  followup = "days_to_last_followup",
  significance = 0.05,
  cutoffs = seq(0, 0.99, 0.01)
)
```

## Arguments

- clinical:

  Data frame: clinical data

- data:

  Numeric: elements of interest to test against the cutoff

- censoring:

  Character: censor using `left`, `right`, `interval` or `interval2`

- event:

  Character: name of column containing time of the event of interest

- timeStart:

  Character: name of column containing starting time of the interval or
  follow up time

- timeStop:

  Character: name of column containing ending time of the interval (only
  relevant for interval censoring)

- followup:

  Character: name of column containing follow up time

- significance:

  Numeric: significance threshold

- cutoffs:

  Numeric: cutoffs to test

## Value

p-value plot

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md),
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
clinical <- read.table(text = "2549   NA ii  female
                                840   NA i   female
                                 NA 1204 iv    male
                                 NA  383 iv  female
                               1293   NA iii   male")
names(clinical) <- c("patient.days_to_last_followup",
                     "patient.days_to_death",
                     "patient.stage_event.pathologic_stage",
                     "patient.gender")
clinical <- do.call(rbind, rep(list(clinical), 5))
rownames(clinical) <- paste("Subject", seq(nrow(clinical)))

# Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
annot <- readFile("ex_splicing_annotation.RDS")
junctionQuant <- readFile("ex_junctionQuant.RDS")

psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 
#> Using 3 of 3 events (100%) whose junctions are present in junction quantification data...
#>   |                                        |   0%   |========                                |  20%   |================                        |  40%   |========================                |  60%   |================================        |  80%   |========================================| 100% 

# Match between subjects and samples
match <- c("Cancer 1"="Subject 3",
           "Cancer 2"="Subject 17",
           "Cancer 3"="Subject 21")

eventData <- assignValuePerSubject(psi[3, ], match)

event      <- "days_to_death"
timeStart  <- "days_to_death"
plotSurvivalPvaluesByCutoff(clinical, eventData, censoring="right",
                            event=event, timeStart=timeStart)

{"x":{"hc_opts":{"chart":{"reflow":true,"height":"100px","zoomType":"x"},"title":{"text":null},"yAxis":{"title":{"text":null},"crosshair":{"color":"gray","width":1,"dashStyle":"shortdash"},"labels":{"enabled":false},"gridLineWidth":0,"plotLines":[{"value":0.05,"color":"Highcharts.getOptions().colors[0]","dashStyle":"shortdash","width":1,"label":{"align":"left","text":"p < 0.05","style":{"color":"Highcharts.getOptions().colors[0]"}}}]},"credits":{"enabled":false},"exporting":{"enabled":false},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0,"cursor":"pointer","point":{"events":{"click":"function () { setPSIcutoffSlider(this.x) }"}},"marker":{"radius":2}},"treemap":{"layoutAlgorithm":"squarified"}},"series":[{"data":[{"x":0,"y":-0,"patients1":null,"patients2":null},{"x":0.01,"y":-0,"patients1":null,"patients2":null},{"x":0.02,"y":-0,"patients1":null,"patients2":null},{"x":0.03,"y":-0,"patients1":null,"patients2":null},{"x":0.04,"y":-0,"patients1":null,"patients2":null},{"x":0.05,"y":-0,"patients1":null,"patients2":null},{"x":0.06,"y":-0,"patients1":null,"patients2":null},{"x":0.07000000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.08,"y":-0,"patients1":null,"patients2":null},{"x":0.09,"y":-0,"patients1":null,"patients2":null},{"x":0.1,"y":-0,"patients1":null,"patients2":null},{"x":0.11,"y":-0,"patients1":null,"patients2":null},{"x":0.12,"y":-0,"patients1":null,"patients2":null},{"x":0.13,"y":-0,"patients1":null,"patients2":null},{"x":0.14,"y":-0,"patients1":null,"patients2":null},{"x":0.15,"y":-0,"patients1":null,"patients2":null},{"x":0.16,"y":-0,"patients1":null,"patients2":null},{"x":0.17,"y":-0,"patients1":null,"patients2":null},{"x":0.18,"y":-0,"patients1":null,"patients2":null},{"x":0.19,"y":-0,"patients1":null,"patients2":null},{"x":0.2,"y":-0,"patients1":null,"patients2":null},{"x":0.21,"y":-0,"patients1":null,"patients2":null},{"x":0.22,"y":-0,"patients1":null,"patients2":null},{"x":0.23,"y":-0,"patients1":null,"patients2":null},{"x":0.24,"y":-0,"patients1":null,"patients2":null},{"x":0.25,"y":-0,"patients1":null,"patients2":null},{"x":0.26,"y":-0,"patients1":null,"patients2":null},{"x":0.27,"y":-0,"patients1":null,"patients2":null},{"x":0.28,"y":-0,"patients1":null,"patients2":null},{"x":0.29,"y":-0,"patients1":null,"patients2":null},{"x":0.3,"y":-0,"patients1":null,"patients2":null},{"x":0.31,"y":-0,"patients1":null,"patients2":null},{"x":0.32,"y":-0,"patients1":null,"patients2":null},{"x":0.33,"y":-0,"patients1":null,"patients2":null},{"x":0.34,"y":-0,"patients1":null,"patients2":null},{"x":0.35,"y":-0,"patients1":null,"patients2":null},{"x":0.36,"y":-0,"patients1":null,"patients2":null},{"x":0.37,"y":-0,"patients1":null,"patients2":null},{"x":0.38,"y":-0,"patients1":null,"patients2":null},{"x":0.39,"y":-0,"patients1":null,"patients2":null},{"x":0.4,"y":-0,"patients1":null,"patients2":null},{"x":0.41,"y":-0,"patients1":null,"patients2":null},{"x":0.42,"y":0.4989407377822485,"patients1":2,"patients2":1},{"x":0.43,"y":0.4989407377822485,"patients1":2,"patients2":1},{"x":0.44,"y":0.4989407377822485,"patients1":2,"patients2":1},{"x":0.45,"y":-0,"patients1":1,"patients2":2},{"x":0.46,"y":-0,"patients1":null,"patients2":null},{"x":0.47,"y":-0,"patients1":null,"patients2":null},{"x":0.48,"y":-0,"patients1":null,"patients2":null},{"x":0.49,"y":-0,"patients1":null,"patients2":null},{"x":0.5,"y":-0,"patients1":null,"patients2":null},{"x":0.51,"y":-0,"patients1":null,"patients2":null},{"x":0.52,"y":-0,"patients1":null,"patients2":null},{"x":0.53,"y":-0,"patients1":null,"patients2":null},{"x":0.54,"y":-0,"patients1":null,"patients2":null},{"x":0.55,"y":-0,"patients1":null,"patients2":null},{"x":0.5600000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.5700000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.58,"y":-0,"patients1":null,"patients2":null},{"x":0.59,"y":-0,"patients1":null,"patients2":null},{"x":0.6,"y":-0,"patients1":null,"patients2":null},{"x":0.61,"y":-0,"patients1":null,"patients2":null},{"x":0.62,"y":-0,"patients1":null,"patients2":null},{"x":0.63,"y":-0,"patients1":null,"patients2":null},{"x":0.64,"y":-0,"patients1":null,"patients2":null},{"x":0.65,"y":-0,"patients1":null,"patients2":null},{"x":0.66,"y":-0,"patients1":null,"patients2":null},{"x":0.67,"y":-0,"patients1":null,"patients2":null},{"x":0.68,"y":-0,"patients1":null,"patients2":null},{"x":0.6900000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.7000000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.71,"y":-0,"patients1":null,"patients2":null},{"x":0.72,"y":-0,"patients1":null,"patients2":null},{"x":0.73,"y":-0,"patients1":null,"patients2":null},{"x":0.74,"y":-0,"patients1":null,"patients2":null},{"x":0.75,"y":-0,"patients1":null,"patients2":null},{"x":0.76,"y":-0,"patients1":null,"patients2":null},{"x":0.77,"y":-0,"patients1":null,"patients2":null},{"x":0.78,"y":-0,"patients1":null,"patients2":null},{"x":0.79,"y":-0,"patients1":null,"patients2":null},{"x":0.8,"y":-0,"patients1":null,"patients2":null},{"x":0.8100000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.8200000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.8300000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.84,"y":-0,"patients1":null,"patients2":null},{"x":0.85,"y":-0,"patients1":null,"patients2":null},{"x":0.86,"y":-0,"patients1":null,"patients2":null},{"x":0.87,"y":-0,"patients1":null,"patients2":null},{"x":0.88,"y":-0,"patients1":null,"patients2":null},{"x":0.89,"y":-0,"patients1":null,"patients2":null},{"x":0.9,"y":-0,"patients1":null,"patients2":null},{"x":0.91,"y":-0,"patients1":null,"patients2":null},{"x":0.92,"y":-0,"patients1":null,"patients2":null},{"x":0.93,"y":-0,"patients1":null,"patients2":null},{"x":0.9400000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.9500000000000001,"y":-0,"patients1":null,"patients2":null},{"x":0.96,"y":-0,"patients1":null,"patients2":null},{"x":0.97,"y":-0,"patients1":null,"patients2":null},{"x":0.98,"y":-0,"patients1":null,"patients2":null},{"x":0.99,"y":-0,"patients1":null,"patients2":null}],"zones":[{"value":1.301029995663981,"color":"lightgray"}]}],"xAxis":{"tickInterval":0.1,"showLastLabel":true,"endOnTick":true,"min":0,"max":1,"minorGridLineWidth":0,"gridLineWidth":0},"legend":[null],"tooltip":{"formatter":"function() { return getPvaluePlotTooltip(this); }"}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":["hc_opts.yAxis.plotLines.0.color","hc_opts.yAxis.plotLines.0.label.style.color","hc_opts.plotOptions.series.point.events.click","hc_opts.tooltip.formatter"],"jsHooks":[]}
```
