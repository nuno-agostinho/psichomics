/* Extend functions available for Highcharts */

/**
 * Highcharts X-range series plugin (used for protein and transcript plot)
 * Modified from: http://jsfiddle.net/Lddp23qj/
*/
(function(H) {
    var defaultPlotOptions = H.getOptions().plotOptions,
        columnType = H.seriesTypes.column,
        each = H.each,
        extendClass = H.extendClass,
        pick = H.pick,
        Point = H.Point,
        Series = H.Series;
    defaultPlotOptions.xrange = H.merge(defaultPlotOptions.column, {
        tooltip: {
            headerFormat: "",
            pointFormat:
                "chr{series.options.chr}, {series.options.strand} strand" +
                "<br/>{point.name} ({point.length} nts): " +
                "{point.x:,.0f} to {point.x2:,.0f}<br/>" + 
                "<span style=\"color:{point.color}\">\u25CF</span>  " +
                "<b>{series.name} ({series.options.display})</b><br/>" +
                "Transcript type: <i>{series.options.biotype}</i>"
        }
    });
    
    H.seriesTypes.xrange = H.extendClass(columnType, {
        pointClass: extendClass(Point, {
            // Add x2 and yCategory to the available properties for tooltip formats
            getLabelConfig: function() {
                var cfg = Point.prototype.getLabelConfig.call(this);
                    cfg.x2 = this.x2;
                    cfg.yCategory = this.yCategory = this.series.yAxis
                        .categories && this.series.yAxis.categories[this.y];
                return cfg;
            }
        }),
        type: 'xrange',
        forceDL: true,
        parallelArrays: ['x', 'x2', 'y', 'width', 'name', 'chr'],
        requireSorting: false,
        animate: H.seriesTypes.line.prototype.animate,
        /* Borrow the column series metrics, but with swapped axes. This gives 
         * access to features like groupPadding, grouping, pointWidth, etc. */
        getColumnMetrics: function() {
            var metrics, chart = this.chart;

            function swapAxes() {
                each(chart.series, function(s) {
                    var xAxis = s.xAxis;
                        s.xAxis = s.yAxis;
                        s.yAxis = xAxis;
                });
            }
            swapAxes();
            this.yAxis.closestPointRange = 1;
            metrics = columnType.prototype.getColumnMetrics.call(this);
            swapAxes();
            return metrics;
        },
        translate: function() {
            columnType.prototype.translate.apply(this, arguments);
            var series = this,
                xAxis = series.xAxis,
                metrics = series.columnMetrics,
                minPointLength = series.options.minPointLength || 0;
            H.each(series.points, function(point) {
                var plotX = point.plotX,
                    plotX2 = xAxis.toPixels(H.pick(
                        point.x2, point.x + (point.len || 0)), true),
                    width = plotX2 - plotX,
                    widthDifference;
                if (minPointLength) {
                    widthDifference = width < minPointLength ? minPointLength -
                        width : 0;
                    plotX -= widthDifference / 2;
                    plotX2 += widthDifference / 2;
                }
                plotX = Math.max(plotX, -10);
                plotX2 = Math.min(Math.max(plotX2, -10), xAxis.len + 10);
                point.shapeArgs = {
                    x: plotX,
                    y: point.plotY + metrics.offset * 1.5 + 40 / point.width,
                    width: plotX2 - plotX > 0 ? plotX2 - plotX : 0,
                    height: point.width
                };
                point.tooltipPos[0] += width / 2 + plotX / 2;
                point.tooltipPos[1] -= metrics.width - metrics.offset;
            });
        }
    });
}(Highcharts));
    
/** 
 * Highcharts sparkline constructor
 */
$(function() {
    Highcharts.SparkLine = function(a, b, c) {
        var hasRenderToArg = typeof a === 'string' || a.nodeName,
            options = arguments[hasRenderToArg ? 1 : 0],
            defaultOptions = {
                chart: {
                    renderTo: (options.chart && options.chart.renderTo) || this
                }
            };

        options = Highcharts.merge(defaultOptions, options);

        return hasRenderToArg ?
            new Highcharts.Chart(a, options, c) :
            new Highcharts.Chart(options, b);
    };
});

/** Position the Highcharts tooltip to a relative point
 * @param {Number} width Tooltip's width
 * @param {Number} height Tooltip's height
 * @param {Object} point Position of interest
 * @return {Object} X and Y coordinates
 */
function tooltipPos (width, height, point) {
    return {
        x: point.plotX - width / 2,
        y: point.plotY - height * 1.2
    };
}

/** Draw sparklines based on the JSON code from the Sparkline HTML element
 */ 
function drawSparklines() {
    var $data = $('sparkline'), sparkline, obj;

    for (var i = 0; i < $data.length; i += 1) {
        sparkline = $($data[i]);
        obj = sparkline.data('sparkline'); // Obtain the JSON code
        obj.tooltip.positioner = tooltipPos; // Correctly position the tooltip
        sparkline.highcharts('SparkLine', obj);
    }
}