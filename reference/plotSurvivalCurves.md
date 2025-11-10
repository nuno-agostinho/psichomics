# Plot survival curves

Plot survival curves

## Usage

``` r
plotSurvivalCurves(
  surv,
  mark = TRUE,
  interval = FALSE,
  pvalue = NULL,
  title = "Survival analysis",
  scale = NULL,
  auto = TRUE
)
```

## Arguments

- surv:

  Survival object

- mark:

  Boolean: mark times?

- interval:

  Boolean: show interval ranges?

- pvalue:

  Numeric: p-value of the survival curves

- title:

  Character: plot title

- scale:

  Character: time scale (default is `days`)

- auto:

  Boolean: return the plot automatically prepared (`TRUE`) or only the
  bare minimum (`FALSE`)?

## Value

Plot of survival curves

## See also

Other functions to analyse survival:
[`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md),
[`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md),
[`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md),
[`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md),
[`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md),
[`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md),
[`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md),
[`survfit.survTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md),
[`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)

## Examples

``` r
require("survival")
fit <- survfit(Surv(time, status) ~ x, data = aml)
plotSurvivalCurves(fit)

{"x":{"hc_opts":{"chart":{"reflow":true,"zoomType":"xy"},"title":{"text":"Survival analysis"},"yAxis":{"title":{"text":"Survival proportion"},"min":0,"max":1,"crosshair":true},"credits":{"enabled":false},"exporting":{"enabled":false},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0,"stickyTracking":false},"treemap":{"layoutAlgorithm":"squarified"},"line":{"marker":{"enabled":false}}},"tooltip":{"pointFormat":"Survival proportion: {point.y:.3f} <br/> Records: {series.options.records} <br/> Events: {series.options.events} <br/> Median: {series.options.median}","headerFormat":"<small>\n  {point.x:.3f}\n  days\n<\/small> <br/> <span style=\"color:{point.color}\">● <\/span> <b>{series.name}<\/b> <br/>"},"series":[{"data":[{"x":0,"y":1},{"x":9,"y":0.9090909090909091,"up":1,"low":0.7541338450815255,"group":"x=Maintained"},{"x":13,"y":0.8181818181818181,"up":1,"low":0.6192489873993635,"group":"x=Maintained","marker":{"fillColor":"black","symbol":"plus","enabled":true}},{"x":18,"y":0.7159090909090908,"up":1,"low":0.4884262874221285,"group":"x=Maintained"},{"x":23,"y":0.6136363636363635,"up":0.9991576002284727,"low":0.3768670595016798,"group":"x=Maintained"},{"x":28,"y":0.6136363636363635,"up":0.9991576002284727,"low":0.3768670595016798,"group":"x=Maintained","marker":{"fillColor":"black","symbol":"plus","enabled":true}},{"x":31,"y":0.4909090909090909,"up":0.9455849552004292,"low":0.2548599511993173,"group":"x=Maintained"},{"x":34,"y":0.3681818181818182,"up":0.8752606781439076,"low":0.1548771178971911,"group":"x=Maintained"},{"x":45,"y":0.3681818181818182,"up":0.8752606781439076,"low":0.1548771178971911,"group":"x=Maintained","marker":{"fillColor":"black","symbol":"plus","enabled":true}},{"x":48,"y":0.1840909090909091,"up":0.943525769474552,"low":0.03591789848918526,"group":"x=Maintained"},{"x":161,"y":0.1840909090909091,"up":0.943525769474552,"low":0.03591789848918526,"group":"x=Maintained","marker":{"fillColor":"black","symbol":"plus","enabled":true}}],"step":"left","name":"x=Maintained","zIndex":1,"color":"Highcharts.getOptions().colors[\n0\n]","records":11,"n.max":11,"n.start":11,"events":7,"rmean":52.64545454545454,"se(rmean)":19.82860279556258,"median":31,"0.95LCL":18,"0.95UCL":null},{"data":[{"x":0,"y":1},{"x":5,"y":0.8333333333333334,"up":1,"low":0.6470369870133617,"group":"x=Nonmaintained"},{"x":8,"y":0.6666666666666667,"up":0.9946253602590915,"low":0.4468460811502639,"group":"x=Nonmaintained"},{"x":12,"y":0.5833333333333334,"up":0.9409980121738093,"low":0.3616137052103847,"group":"x=Nonmaintained"},{"x":16,"y":0.5833333333333334,"up":0.9409980121738093,"low":0.3616137052103847,"group":"x=Nonmaintained","marker":{"fillColor":"black","symbol":"plus","enabled":true}},{"x":23,"y":0.4861111111111112,"up":0.8833192253395576,"low":0.2675182488582666,"group":"x=Nonmaintained"},{"x":27,"y":0.388888888888889,"up":0.8157357156912723,"low":0.1853965260955568,"group":"x=Nonmaintained"},{"x":30,"y":0.2916666666666667,"up":0.7408220185158038,"low":0.1148311501524704,"group":"x=Nonmaintained"},{"x":33,"y":0.1944444444444445,"up":0.6642236598825634,"low":0.05692155257160417,"group":"x=Nonmaintained"},{"x":43,"y":0.09722222222222224,"up":0.6195486291190556,"low":0.0152565271708652,"group":"x=Nonmaintained"},{"x":45,"y":0,"up":null,"low":null,"group":"x=Nonmaintained"}],"step":"left","name":"x=Nonmaintained","zIndex":1,"color":"Highcharts.getOptions().colors[\n1\n]","records":12,"n.max":12,"n.start":12,"events":11,"rmean":22.70833333333334,"se(rmean)":4.180941981033154,"median":23,"0.95LCL":8,"0.95UCL":null}],"xAxis":{"title":{"text":"Time in days"},"crosshair":true}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":["hc_opts.series.0.color","hc_opts.series.1.color"],"jsHooks":[]}
```
