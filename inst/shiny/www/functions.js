/* Check all groups when clicking the checkbox in the group table's header*/
/* TODO(NunoA): Make this work... */
$("input[name='checkAllGroups']").change(function () {
    $("input[name='checkGroups']").prop('checked', $(this).prop("checked"));
});

/* Get which checkboxes are checked in the groups section */
Shiny.addCustomMessageHandler('getCheckedBoxes', function(variable) {   
    var selected = [];
    $("input[name='checkGroups']:checked").each(function() {
        selected.push($(this).attr('number'));
        
    });
    /* Add value to variable in R */
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