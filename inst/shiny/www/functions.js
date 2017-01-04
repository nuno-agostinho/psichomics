/**
 * Update location according to browser navigation
 * 
 * Inspired by code available at
 * https://github.com/daattali/advanced-shiny/blob/master/navigate-history
 */
shinyjs.init = function() {
    window.onpopstate = function (event) {
        Shiny.onInputChange('appLocation', location.search);
    }
}

/**
 * Update URL to reflect the current browser navigation
 * @param {Object} params Pair key-value object containing URL parameters
 * 
 * Inspired by code available at
 * https://github.com/daattali/advanced-shiny/blob/master/navigate-history
 */
updateHistory = function(params) {
    var queryString = [];
    for (var key in params) {
        queryString.push(encodeURIComponent(key) + '=' + 
            encodeURIComponent(params[key]));
    }
    queryString = '?' + queryString.join('&');
    history.pushState(null, null, queryString)
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
 * @param {String} event Alternative splicing event
 */
function showSurvCutoff(event) {
    if (event !== null) changeEvent(event);
    
    var surv = "Survival analysis";
    $("ul[id='nav'] > li > ul > li > a[data-value*='" + surv + "']")
        .tab("show");
    $("input[value='psiCutoff']").click();
}

/**
 * Change selected event
 * @param {String} event Alternative splicing event
 */
function changeEvent (event) {
    $("select[id*='selectizeEvent']").selectize()[0].selectize.setValue(event);
}

/**
 * Navigate user to differential splicing of a given alternative splicing event
 * @param {String} event Alternative splicing event
 */
function showDiffSplicing (event) {
    var tabName = "Differential splicing analysis";
    $("ul[id='nav'] > li > ul > li > a[data-value*='" + tabName + "']")
        .tab("show");
    
    var diff = "Single event";
    $("a[data-value*='" + diff + "']").tab("show");
    changeEvent(event);
}

/**
 * Navigate user to the selectize input where to change the clinical groups used
 * for differential splicing analysis
 */
function changeDiffSplicingGroup () {
    var diff = "All events (table)";
    $("a[data-value*='" + diff + "']").tab("show");
    $("#analyses-diffSplicing-diffSplicingTable-diffGroups")[0].selectize
        .focus();
}

/**
 * Modify row of a table to include links to navigate user to the differential
 * splicing of the respective event
 * 
 * @param {Numeric} row Row to introduce links
 * @param {Object} data Table of interest
 * 
 * @return Same rows from input with links
 */
function createDiffSplicingLinks(row, data, index) {
    $('td:eq(0)', row).html("<a onclick='showDiffSplicing(\"" + data[0] + 
        "\")' href='javascript:void(0);'>" + data[0] + "</a>");
    return row;
}

/* Get which checkboxes are checked in the groups section */
Shiny.addCustomMessageHandler('getCheckedBoxes', function(variable) {   
    var selected = getSelectedCheckboxes();
    
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

$.fn.extend({
    /**
     * Play animation from Animate.css in a selected element
     * 
     * @param {String} animationName Name of animation of interest
     * @example $('div').animateCss("bounce");
     */
    animateCss: function (animationName) {
        var animationEnd = 'webkitAnimationEnd mozAnimationEnd MSAnimationEnd' +
            'oanimationend animationend';
        $(this).addClass('animated ' + animationName)
            .one(animationEnd, function() {
                $(this).removeClass('animated ' + animationName);
            });
    }
});