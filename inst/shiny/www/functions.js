/**
 * Update location according to browser navigation
 * 
 * Inspired by code available at
 * https://github.com/daattali/advanced-shiny/blob/master/navigate-history
 */
shinyjs.init = function() {
    window.onpopstate = function (event) {
        Shiny.onInputChange('appLocation', location.search);
    };
};

/**
 * Update URL to reflect the current browser navigation
 * @param {Object} params Pair key-value object containing URL parameters
 * 
 * Inspired by code available at
 * https://github.com/daattali/advanced-shiny/blob/master/navigate-history
 */
function updateHistory(params) {
    var queryString = [];
    for (var key in params) {
        queryString.push(encodeURIComponent(key) + '=' + 
            encodeURIComponent(params[key]));
    }
    queryString = '?' + queryString.join('&');
    history.pushState(null, null, queryString);
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
    $("#nav > li > a[data-value*='Data']").click();
    
    // Collapse data panels
    $("#data-accordion > div > div[class*='panel-collapse']").collapse('hide');
}

/**
 * Change selected event
 * @param {String} event Alternative splicing event
 */
function changeEvent (event) {
    $("#selectizeEventElem").selectize()[0].selectize.setValue(event);
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
 * @param {Boolean} autoParams Automatically set expected parameters based on
 * the choices for the exploratory differential analyses
 * @param {String} groupSelectize Identifier of the group selection element to
 * automatically set expected parameters based on its values
 */
function showDiffSplicing (event, autoParams = false, groupSelectize = null) {
    // Navigate to differential splicing analyses for a single event
    var tabName = "Differential splicing analysis";
    $("#nav > li > ul > li > a[data-value*='" + tabName + "']").click();
    var diff = "Single event";
    $("a[data-value*='" + diff + "']").click();
    
    // Change currently selected splicing event
    changeEvent(event);
    
    if (autoParams) {
        groupSelectize = "analyses-diffSplicing-diffSplicingTable-diffGroups";
    } else if (groupSelectize === null) {
        $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
        return;
    }
    
    // Set whether using groups or not
    singleEventPage = "analyses-diffSplicing-diffSplicingEvent";
    groups = $("input[type='radio'][name='" + groupSelectize +
        "Selection']:checked")[0].value;
    $("input[type='radio'][name='" + singleEventPage +
        "-diffGroupsSelection'][value=" + groups + "]").click();
        
    if (groups == "groups") {
        // Set selected groups
        items = $("#" + groupSelectize)[0].selectize.items;
        $("#" + singleEventPage + "-diffGroups")[0].selectize.setValue(items);
    }
    
    // Perform statistical analyses for a single event
    $("#" + singleEventPage + "-analyse")[0].click();
    $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
}

/**
 * Navigate user to survival analysis by quantification cutoff
 * @param {String} event Alternative splicing event
 * @param {Boolean} autoParams Automatically set expected parameters
 */
function showSurvCutoff(event, autoParams = false) {
    // Change currently selected splicing event
    if (event !== null) changeEvent(event);
    
    // Navigate to survival analyses
    var surv = "Survival analysis";
    $("#nav > li > ul > li > a[data-value*='" + surv + "']").click();
    
    var survivalPage = "analyses-survival";
    if (autoParams) {
        // Perform survival analyses once the optimal PSI is calculated
        $("#" + survivalPage + "-psiCutoff").one('change', function(){
            setTimeout(function() {
                $("#" + survivalPage + "-survivalCurves")[0].click();
            }, 500);
        });
    }
    
    // Set PSI cutoff
    $("input[value='psiCutoff']").click();
    
    if (autoParams && event !== null) {
        var allEventsPage = "analyses-diffSplicing-diffSplicingTable";
        
        // Set censoring interval
        censoring = $("input[type='radio'][name='" + allEventsPage +
            "-censoring']:checked")[0].value;
        $("input[type='radio'][name='" + survivalPage + "-censoring'][value=" +
            censoring + "]").click();
        
        // Set follow up time or starting time
        timeStart = $("#" + allEventsPage + "-timeStart")[0].selectize.items;
        $("#" + survivalPage + "-timeStart")[0].selectize.setValue(timeStart);
            
        // Set event of interest
        event = $("#" + allEventsPage + "-event")[0].selectize.items;
        $("#" + survivalPage + "-event")[0].selectize.setValue(event);
            
        if (censoring == "interval" || censoring == "interval2") {
            // Set ending time
            timeStop = $("#" + allEventsPage + "-timeStop")[0].selectize.items;
            $("#" + survivalPage + "-timeStop")[0].selectize.setValue(event);
        }
    }
    $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
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
        "\", autoParams=true)' href='javascript:void(0);' " + 
        "title='Differential splicing analyses for " + event + "'>" + event +
        "</a>");
    return row;
}

/**
 * Update slider value for PSI cutoff
 * 
 * @param {Numeric} value Slider value
 */
function setPSIcutoffSlider(value) {
    $("input[id*='psiCutoff']").data("ionRangeSlider").update({from: value});
}

/**
 * Get tooltip for p-value plot in survival analysis
 * 
 * @param {Object} object Tooltip object
 * 
 * @return HTML object that renders the appropriate tooltip
 */
function getPvaluePlotTooltip(object) {
    head = 'PSI cutoff: ' + object.point.x + '<br/>';
    minuslog10pvalue = '-log₁₀(p-value): ' + object.point.y.toFixed(3) +
        ' <span style=\"color:' + object.point.color + '\">' +
        '\u25CF</span><br/>';
    pvalue = 'p-value: ' + Math.pow(10, -object.point.y).toFixed(3) + "<br/>";
    
    patients = '';
    if(object.point.patients2 !== null) {
        patients = object.point.patients1 + ' vs ' + object.point.patients2 +
            ' patients';
    }
    
    return head + minuslog10pvalue + pvalue + patients;
}

/* Change document title to reflect whether the app is busy */
setInterval(function() {
    document.title = ($('html').hasClass('shiny-busy')) ?
        '🤔 psichomics' : 'psichomics';
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

(function() {
 /* Shiny Registration */
    var fileBrowserInputBinding = new Shiny.InputBinding();
    $.extend(fileBrowserInputBinding, {
        find: function(scope) {
            return( $(scope).find(".fileBrowser-input") );
        },
        getId: function(el) { return($(el).attr('id')); },
        getValue: function(el) {
 return($(el).data('val') || 0); },
        setValue: function(el, value) { $(el).data('val', value); },
        receiveMessage: function(el, data) {
            // This is used for receiving messages that tell the input object to do
            // things, such as setting values (including min, max, and others).
            var $widget = $(el).parentsUntil('.fileBrowser-input-container')
                .parent();
            var $path = $widget.find('input.fileBrowser-input-chosen-dir');

            if (data.path) {                $path.val(data.path);
                $path.trigger('change');
            }
        },
        subscribe: function(el, callback) {
            $(el).on("click.fileBrowserInputBinding", function(e) {
                var $el = $(this);
                var val = $el.data('val') || 0;
                $el.data('val', val + 1);
                callback();
            });
        },
        unsubscribe: function(el) { $(el).off(".fileBrowserInputBinding"); }
    });
    Shiny.inputBindings
        .register(fileBrowserInputBinding,
            "oddhypothesis.fileBrowsernputBinding");
})();