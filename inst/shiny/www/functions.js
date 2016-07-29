/** Check all groups when clicking the checkbox in the group table's header */
function checkAllGroups() {
    $("input[name='checkAllGroups']").change(function () {
        $("input[name='checkGroups']").prop('checked', $(this).prop("checked"));
    });
}

/**
 * Change active tab to the Data panel and expand the panel with the given value
 * @param {String} panelVal Value of the panel to open
 */
function showDataPanel(panelVal) {
    // Open Data tab
    $("ul[id='nav'] > li > a[data-value*='Data']").tab("show");

    // Expand panel of interest
    $("div[value*='" + panelVal + "'] > div[role='tabpanel']").collapse("show");
}

/**
 * Navigate user to survival analysis by quantification cut-off
 */
function showSurvCutoff() {
    var surv = "Survival curves";
    $("select[id*='selectizeAnalysis']").selectize()[0].selectize.setValue(surv);
    $("input[value='psiCutoff']").click();
}

/**
 * Navigate user to differential splicing of a given alternative splicing event
 * @param {String} event Alternative splicing event
 */
function showDiffSplicing (event) {
    var diff = "Differential analysis (per splicing event)";
    $("select[id*='selectizeAnalysis']").selectize()[0].selectize.setValue(diff);
    $("select[id*='selectizeEvent']").selectize()[0].selectize.setValue(event);
}

/**
 * Modify row of a table to include links to navigate user to the differential
 * splicing of the respective event
 * 
 * @param {Numeric} row Row to introduce links
 * @param {Object} data Table of interest
 * 
 * @return Links
 */
function createDiffSplicingLinks(row, data, index) {
    $('td:eq(0)', row).html("<a onclick='showDiffSplicing(\"" + data[0] + 
        "\")' href='javascript:void(0);'>" + data[0] + "</a>");
    return row;
}

/* Get which checkboxes are checked in the groups section */
Shiny.addCustomMessageHandler('getCheckedBoxes', function(variable) {   
    var selected = [];
    $("input[name='checkGroups']:checked").each(function() {
        selected.push($(this).attr('number'));
        
    });
    // Add value to variable in R
    Shiny.onInputChange(variable, selected);
});

/* Set given variable to zero */
Shiny.addCustomMessageHandler('setZero', function(variable) {
    Shiny.onInputChange(variable, 0);
});

/* Change document title to reflect if app is busy */
setInterval(function() {
    document.title = ($('html').hasClass('shiny-busy')) ?
        '[R] PSΨchomics' : 'PSΨchomics';
    }, 500);

/* Highcharts sparkline constructor */
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