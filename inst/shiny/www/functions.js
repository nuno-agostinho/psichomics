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
 * Change active tab to the Data panel and collapse data panels
 * @param {String} modal Identifier of the modal to close (optional)
 */
function showDataPanel(modal) {
    if (typeof myVariable === 'undefined' && modal !== null) {
        $(modal).modal("hide"); // Close modal
    }
    
    // Open Data tab
    $("ul[id='nav'] > li > a[data-value*='Data']").click();
    
    // Collapse data panels
    $("div[id='data-accordion'] > div > div[class*='panel-collapse']")
        .collapse('hide');
}

/**
 * Change selected event
 * @param {String} event Alternative splicing event
 */
function changeEvent (event) {
    $("select[id*='selectizeEvent']").selectize()[0].selectize.setValue(event);
}

/**
 * Set selected transcript
 * @param {String} transcript Transcript identifier
 */
function setTranscript (transcript) {
    $("select[id*='selectedTranscript']").selectize()[0].selectize
        .setValue(transcript);
}

/**
 * Navigate user to differential splicing of a given alternative splicing event
 * @param {String} event Alternative splicing event
 * @param {Boolean} autoParams Automatically set expected parameters
 */
function showDiffSplicing (event, autoParams = false) {
    // Navigate to differential splicing analyses for a single event
    var tabName = "Differential splicing analysis";
    $("ul[id='nav'] > li > ul > li > a[data-value*='" + tabName + "']").click();
    var diff = "Single event";
    $("a[data-value*='" + diff + "']").click();
    
    // Change currently selected splicing event
    changeEvent(event);
    
    if (!autoParams) { return; }
    
    // Set whether using groups or not
    allEventsPage = "analyses-diffSplicing-diffSplicingTable";
    singleEventPage = "analyses-diffSplicing-diffSplicingEvent";
    groups = $("input[type='radio'][name='" + allEventsPage +
        "-diffGroupsSelection']:checked")[0].value;
    $("input[type='radio'][name='" + singleEventPage +
        "-diffGroupsSelection'][value=" + groups + "]").click();
        
    if (groups == "groups") {
        // Set selected groups
        items = $("select[id='" + allEventsPage + "-diffGroups']")[0].selectize
            .items;
        $("select[id='" + singleEventPage + "-diffGroups']")[0].selectize
            .setValue(items);
    }
    
    // Perform statistical analyses for a single event
    $("button[id='" + singleEventPage + "-analyse']")[0].click();
}

/**
 * Navigate user to survival analysis by quantification cut-off
 * @param {String} event Alternative splicing event
 * @param {Boolean} autoParams Automatically set expected parameters
 */
function showSurvCutoff(event, autoParams = false) {
    // Change currently selected splicing event
    if (event !== null) changeEvent(event);
    
    // Navigate to survival analyses
    var surv = "Survival analysis";
    $("ul[id='nav'] > li > ul > li > a[data-value*='" + surv + "']").click();
    
    if (!autoParams) { return; }
    
    if (event !== null) {
        allEventsPage = "analyses-diffSplicing-diffSplicingTable";
        survivalPage = "analyses-survival";
        
        // Set censoring interval
        censoring = $("input[type='radio'][name='" + allEventsPage +
            "-censoring']:checked")[0].value;
        $("input[type='radio'][name='" + survivalPage + "-censoring'][value=" +
            censoring + "]").click();
        
        // Set follow up time or starting time
        timeStart = $("select[id='" + allEventsPage + "-timeStart']")[0]
            .selectize.items;
        $("select[id='" + survivalPage + "-timeStart']")[0].selectize
            .setValue(timeStart);
            
        // Set event of interest
        event = $("select[id='" + allEventsPage + "-event']")[0].selectize
            .items;
        $("select[id='" + survivalPage + "-event']")[0].selectize
            .setValue(event);
            
        if (censoring == "interval" || censoring == "interval2") {
            // Set ending time
            timeStop = $("select[id='" + allEventsPage + "-timeStop']")[0]
                .selectize.items;
            $("select[id='" + survivalPage + "-timeStop']")[0].selectize
                .setValue(event);
        }
    }
    $("input[value='psiCutoff']").click();
    // Perform survival analyses
    setTimeout(function() {
        $("button[id='" + survivalPage + "-survivalCurves']")[0].click();
    }, 2000);
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
    var event = data[0];
    var eventID = event.replace(/ /g, "_");
    
    $('td:eq(0)', row).html("<a onclick='showDiffSplicing(\"" + eventID +
        "\")' href='javascript:void(0);' " + 
        "title='Differential splicing analyses for " + event + "'>" + event +
        "</a>");
    return row;
}

/* Change document title to reflect whether the app is busy */
setInterval(function() {
    document.title = ($('html').hasClass('shiny-busy')) ?
        '[R] PSΨchomics' : 'PSΨchomics';
    }, 500);

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