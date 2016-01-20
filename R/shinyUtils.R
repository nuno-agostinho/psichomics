badge <- function (inputId, label)
    span(class = "badge", id = inputId, label)

buttonGroups <- function(...) {
    div(class="btn-group", role="group", ...)
}

progressbar <- function(value, min = 0, max = 100, label = paste0(value, "%"),
                        striped = T) {
    stripedText <- ifelse(striped, "progress-bar-striped", "")
    div(class="progress",
        div(class=paste("progress-bar", stripedText),
            role="progressbar",
            "aria-valuenow"=value,
            "aria-valuemin"=min,
            "aria-valuemax"=max,
            style=paste0("width: ", value, "%;"),
            label
            )
    )
}

alertNew <- function(...) {
    div(class="alert alert-warning alert-dismissible", role="alert", ...)
}

dropdown <- function(inputId) {
    div(class="dropdown",
        tags$button(class="btn btn-default dropdown-toggle",
                    type="button", id=inputId,
                    "data-toggle"="dropdown",
                    "aria-haspopup"="true",
                    "aria-expanded"="true",
                    list("Dropdown", span(class="caret"))),
        tags$ul(class="dropdown-menu", "aria-labelledby"="dropdownMenu1",
                tags$li(a(href="#", "Action")),
                tags$li(a(href="#", "Another action")),
                tags$li(role="separator", class="divider"),
                tags$li(a(href="#", "Something else")))
        )
}