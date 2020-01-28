/** psichomics.js */

/* Ensure code escaping */
window.escape = window.escape || window.encodeURI;

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
 * Change active tab to the given tab of the Groups panel
 * @param {String} tab Tab to open
 */
function showGroups(type) {
    $("a[data-value='Groups']")[0].click();
    
    var mode;
    if (type === "Samples" | type === "Patients") {
        mode = 0;
    } else if (type === "ASevents" || type === "Genes") {
        mode = 1;
    }
    $("#groupsTypeTab a")[mode].click();
}

/**
 * Render group DataTable
 * 
 * @param table DataTable
 */
function renderGroupTable(table) {
    var getIcon = function(symbol) {
        return '<i class="fa fa-' + symbol + '" aria-hidden="true"></i>';
    };
    var plusIcon  = getIcon("plus-circle");
    var minusIcon = getIcon("minus-circle");
    
    var cols  = table.columns()[0].slice(-2);
    table.columns(cols).visible(false, false);
    table.columns.adjust().draw(false);
    
    var format = function(data) {
        return '<table class="table table-details" border="0">' + '<tr>' + 
        '<td>Subset:</td>' + '<td>'+ data[data.length - 2] + '</td>' + '</tr>' + 
        '<tr>' + '<td>Input:</td>' + '<td>' + data[data.length - 1] + '</td>' + 
        '</tr>' + '</table>';
    };
    
    table.on('click', 'td.details-control', function() {
        var td = $(this), row = table.row(td.closest('tr'));
        if (row.child.isShown()) {
            // Hide extra information
            row.child.hide();
            td.html(plusIcon);
        } else {
            // Show extra information
            row.child( format(row.data()), 'no-padding' ).show();
            td.html(minusIcon);
        }
    });
}

/**
 * Prepare interface for group selection
 */
function renderGroupSelection (item, escape) {
    var description =  item.label.split(" #")[0];
    var colour = "#" + item.label.split(" #")[1];
    return "<div><b><font color='" + colour + "'>\u25CF</font> " + 
        escape(item.value) + "</b><small> " + escape(description) + 
        "</small></div>";
}

/**
 * Change selected alternative splicing event
 * @param {String} event Alternative splicing event
 */
function changeEvent (event) {
    $("#selectizeEventElem").selectize()[0].selectize.setValue(event);
}

/**
 * Render alternative splicing event
 */
function renderEvent (item, escape) {
    var tmp    = item.label.split(" __ "),
        parsed = tmp[0].split("_"),
        type   = parsed[0],
        pos    = parsed[1],
        gene   = parsed.slice(-1)[0],
        svg    = tmp[1];
    if (svg === undefined) svg = "";
    return `<div>${svg}<b>${gene}</b> ${pos}<br>${type}</div>`;
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
 * Render gene expression normalisation options
 */
function renderGEnormOptions (item, escape) {
    var description;
    switch(item.value) {
        case "TMM":
            description = "This method is recommended" +
                          " for most RNAseq data where more " +
                          "than half of the genes are believed " +
                          "not differentially expressed " + 
                          "between any pair of the samples.";
            break;
        case "RLE":
            description = "The median library is " +
                          "calculated from the geometric mean " +
                          "of all columns and the median ratio " +
                          "of each sample to the median library" +
                          " is taken as the scale factor.";
            break;
        case "upperquartile":
            description = "The scale factors are " +
                          "calculated from a given quantile of " +
                          "the counts for each library, after " +
                          "removing genes with zero counts in " +
                          "all libraries.";
            break;
        case "none":
            description = "";
            break;
        case "quantile":
            description = "Forces the entire empirical distribution of each " +
                          "column to be identical (only performed if " +
                          "<i>voom</i> is selected).";
    }
    return "<div><b>" + escape(item.label) + "</b></br>" +
        "<small>" + description + "</small></div>"; 
}

/**
 * Navigate user to differential splicing of a given alternative splicing event
 * @param {String} event Alternative splicing event
 * @param {String} groups List of groups used for differential analysis
 */
function showDiffSplicing (event, groups = null) {
    var singleEventPage = "analyses-diffSplicing-diffSplicingEvent";
    
    // Navigate to differential splicing analyses for a single event
    var diff = "Individual alternative splicing event";
    $("a[data-value*='" + diff + "']").click();
    
    // Change currently selected splicing event
    changeEvent(event);
        
    if (groups !== null) {
        // Set selected groups
        $("input[type='radio'][name='" + singleEventPage +
            "-diffGroupsSelection'][value=groups]").click();
        $("#" + singleEventPage + "-diffGroups")[0].selectize.setValue(groups);
    } else {
        $("input[type='radio'][name='" + singleEventPage +
            "-diffGroupsSelection'][value=noGroups]").click();
    }
    
    // Perform statistical analyses for a single event
    $("#" + singleEventPage + "-analyse")[0].click();
    $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
}

/**
 * Navigate user to differential expression of a given gene
 * @param {String} gene Gene symbol
 * @param {String} groups List of groups used for differential analysis
 * @param {String} geneExpr Gene expression dataset name
 */
function showDiffExpression (gene, groups = null, geneExpr = null) {
    var singleEventPage = "analyses-diffExpression-diffExpressionEvent";
    
    // Navigate to differential analyses for a single gene
    var page = "Individual gene";
    $("a[data-value*='" + page + "']").click();
    
    // Set selected gene expression
    if (geneExpr !== null) {
        var geneExprSel = $("#" + singleEventPage + "-geneExpr")[0].selectize;
        geneExprSel.addOption({label: geneExpr, value: geneExpr});
        geneExprSel.refreshOptions(false);
        geneExprSel.addItem(geneExpr);
    }
    
    // Set selected gene
    if (gene !== null) {
        var geneSel = $("#" + singleEventPage + "-gene")[0].selectize;
        geneSel.addOption({label: gene, value: gene});
        geneSel.refreshOptions(false);
        geneSel.addItem(gene);
    }
    
    // Set selected groups
    if (groups !== null) {
        $("input[type='radio'][name='" + singleEventPage +
            "-diffGroupsSelection'][value=groups]").click();
        $("#" + singleEventPage + "-diffGroups")[0].selectize.setValue(groups);
    }
    
    // Perform statistical analyses for a single event
    $("#" + singleEventPage + "-analyse")[0].click();
    $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
}

/**
 * Prepare a cell from DataTable with a link to survival
 * 
 * @param data Cell data from DataTable
 * @param type Type of data from DataTable
 * @param row Row from DataTable
 */
function linkToShowSurv(data, type, row) {
    if (type === 'display') {
        console.log(row[0]);
        var event = row[0].replace(new RegExp(" ", 'g'), "_");
        event = '\'' + event + '\'';
        data  = '<a href="javascript:void(0)" onclick="showSurvCutoff(' +
            event + ')">' + data + '</a>';
    }
    return data;
}

/**
 * Navigate user to survival analysis by a value cutoff
 * @param {String} event Alternative splicing event
 * @param {String} groups List of groups used for survival analysis
 * @param {Boolean} autoParams Automatically set expected parameters
 * @param {Boolean} psi Prepare splicing quantification (true) or gene 
 * expression cutoff (false)?
 */
function showSurvCutoff(event, groups = null, autoParams = true, psi = true) {
    // Change currently selected splicing event
    if (event !== null && psi) changeEvent(event);
    
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
    
    if ( psi ) {
        $("input[value='psiCutoff']").click();
    } else {
        $("input[value='geCutoff']").click();
        
        // Set selected gene
        var singleEventPage = "analyses-diffExpression-diffExpressionEvent";
        var gene = $("#" + singleEventPage + "-gene")[0].selectize.items[0];
        var geneSel = $("#" + survivalPage + "-gene")[0].selectize;
        geneSel.addOption({label: gene, value: gene});
        geneSel.refreshOptions(false);
        geneSel.addItem(gene);
    }
    
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
        
        // Set selected groups
        $("#" + survivalPage + "-sampleFiltering")[0].selectize.setValue(
            groups);
    }
    $('html, body').animate({ scrollTop: 0 }, 'slow'); // Scroll to top
}

/**
 * Update slider value for GE cutoff
 * 
 * @param {Numeric} value Slider value
 */
function setGEcutoffSlider(value) {
    $("input[id*='geCutoff']").data("ionRangeSlider").update({from: value});
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
    pvalue = 'p-value: ' + Math.pow(10, -object.point.y).toFixed(3) + 
        ' <span style=\"color:' + object.point.color + '\">' +
        '\u25CF</span><br/>';
    minuslog10pvalue = '-log₁₀(p-value): ' + object.point.y.toFixed(3) + 
        "<br/>";
    patients = '';
    if (object.point.patients2 !== null) {
        patients = object.point.patients1 + ' vs ' + object.point.patients2 +
            ' patients';
    }
    
    return head + pvalue + minuslog10pvalue + patients;
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
        getValue: function(el) { return($(el).data('val') || 0); },
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