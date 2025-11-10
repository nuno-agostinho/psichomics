# Render boxplot

Render boxplot

## Usage

``` r
renderBoxplot(
  data,
  outliers = FALSE,
  sortByMedian = TRUE,
  showXlabels = TRUE,
  title = NULL,
  seriesName = "Gene expression"
)
```

## Arguments

- data:

  Data frame or matrix

- outliers:

  Boolean: draw outliers?

- sortByMedian:

  Boolean: sort box plots based on ascending median?

- showXlabels:

  Boolean: show labels in X axis?

## Value

Box plot

## Examples

``` r
psichomics:::renderBoxplot(data.frame(a=1:10, b=10:19, c=45:54))

{"x":{"hc_opts":{"chart":{"reflow":true,"zoomType":"x","type":"column"},"title":[],"yAxis":{"title":{"text":null},"min":0},"credits":{"enabled":false},"exporting":{"enabled":true,"formAttributes":{"target":"_blank"},"buttons":{"contextButton":{"text":"Export","menuItems":[{"text":"PNG image","onclick":"function () { this.exportChart({ type: 'image/png' }); }"},{"text":"JPEG image","onclick":"function () { this.exportChart({ type: 'image/jpeg' }); }"},{"text":"SVG vector image","onclick":"function () { this.exportChart({ type: 'image/svg+xml' }); }"},{"text":"PDF document","onclick":"function () { this.exportChart({ type: 'application/pdf' }); }"},{"separator":true},{"text":"CSV document","onclick":"function () { this.downloadCSV(); }"},{"text":"XLS document","onclick":"function () { this.downloadXLS(); }"}],"theme":{"fill":"transparent"}}}},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0},"treemap":{"layoutAlgorithm":"squarified"},"boxplot":{"color":"gray","fillColor":"orange"}},"series":[{"name":"Gene expression","data":[{"name":"a","low":1,"q1":3,"median":5.5,"q3":8,"high":10},{"name":"b","low":10,"q1":12,"median":14.5,"q3":17,"high":19},{"name":"c","low":45,"q1":47,"median":49.5,"q3":52,"high":54}],"id":null,"type":"boxplot"}],"xAxis":{"type":"category","labels":{"enabled":true},"visible":true}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":["hc_opts.exporting.buttons.contextButton.menuItems.0.onclick","hc_opts.exporting.buttons.contextButton.menuItems.1.onclick","hc_opts.exporting.buttons.contextButton.menuItems.2.onclick","hc_opts.exporting.buttons.contextButton.menuItems.3.onclick","hc_opts.exporting.buttons.contextButton.menuItems.5.onclick","hc_opts.exporting.buttons.contextButton.menuItems.6.onclick"],"jsHooks":[]}
```
