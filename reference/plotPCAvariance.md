# Create the explained variance plot from a PCA

Create the explained variance plot from a PCA

## Usage

``` r
plotPCAvariance(pca)
```

## Arguments

- pca:

  `prcomp` object

## Value

Plot variance as an `highchart` object

## See also

Other functions to analyse principal components:
[`calculateLoadingsContribution()`](https://nuno-agostinho.github.io/psichomics/reference/calculateLoadingsContribution.md),
[`performPCA()`](https://nuno-agostinho.github.io/psichomics/reference/performPCA.md),
[`plotPCA()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCA.md)

## Examples

``` r
pca <- prcomp(USArrests)
plotPCAvariance(pca)

{"x":{"hc_opts":{"chart":{"reflow":true,"zoomType":"xy"},"title":{"text":"Variance explained by each Principal Component (PC)"},"yAxis":{"title":{"text":"Percentage of variance"},"min":0,"max":100},"credits":{"enabled":false},"exporting":{"enabled":true,"formAttributes":{"target":"_blank"},"buttons":{"contextButton":{"text":"Export","menuItems":[{"text":"PNG image","onclick":"function () { this.exportChart({ type: 'image/png' }); }"},{"text":"JPEG image","onclick":"function () { this.exportChart({ type: 'image/jpeg' }); }"},{"text":"SVG vector image","onclick":"function () { this.exportChart({ type: 'image/svg+xml' }); }"},{"text":"PDF document","onclick":"function () { this.exportChart({ type: 'application/pdf' }); }"},{"separator":true},{"text":"CSV document","onclick":"function () { this.downloadCSV(); }"},{"text":"XLS document","onclick":"function () { this.downloadXLS(); }"}],"theme":{"fill":"transparent"}}}},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0,"dataLabels":{"format":"{point.eigenvalue:.2f}<br/>{point.y:.2f}%","align":"center","verticalAlign":"top","enabled":true}},"treemap":{"layoutAlgorithm":"squarified"}},"series":[{"data":[{"y":96.55342205668825,"eigenvalue":7011.114851023603,"cumvar":96.55342205668825},{"y":2.781733663217495,"eigenvalue":201.9923663226134,"cumvar":99.33515571990574},{"y":0.579953492234191,"eigenvalue":42.11265075533881,"cumvar":99.91510921213994},{"y":0.08489078786007123,"eigenvalue":6.164246184163202,"cumvar":100}],"type":"waterfall","cumvar":[96.55342205668825,99.33515571990574,99.91510921213994,100]}],"xAxis":{"title":{"text":"Principal Components"},"categories":[1,2,3,4],"crosshair":true},"legend":{"enabled":false},"tooltip":{"headerFormat":"<b>Principal component {point.x}<\/b> <br/>","pointFormat":"Eigenvalue: {point.eigenvalue:.2f}<br/>Variance: {point.y:.2f}%<br/>Cumulative variance: {point.cumvar:.2f}%"}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":["hc_opts.exporting.buttons.contextButton.menuItems.0.onclick","hc_opts.exporting.buttons.contextButton.menuItems.1.onclick","hc_opts.exporting.buttons.contextButton.menuItems.2.onclick","hc_opts.exporting.buttons.contextButton.menuItems.3.onclick","hc_opts.exporting.buttons.contextButton.menuItems.5.onclick","hc_opts.exporting.buttons.contextButton.menuItems.6.onclick"],"jsHooks":[]}
```
