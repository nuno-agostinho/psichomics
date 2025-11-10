# User interface

The user interface (UI) controls the layout and appearance of the app.
All CSS modifications are in the file `shiny/www/styles.css`

## Usage

``` r
appUI()

analysesUI(id, tab)

diffEventUI(id, ns, psi = TRUE)

correlationUI(id)

diffExpressionUI(id, tab)

diffExpressionEventUI(id)

diffExpressionTableUI(id)

diffSplicingUI(id, tab)

diffSplicingEventUI(id)

diffSplicingTableUI(id)

dimReductionUI(id, tab)

icaUI(id)

pcaUI(id)

infoUI(id)

survivalUI(id)

templateUI(id)

dataUI(id, tab)

firebrowseUI(id, panel)

geNormalisationFilteringUI(id, panel)

gtexDataUI(id, panel)

inclusionLevelsUI(id, panel)

inclusionLevelsFilterUI(id, panel)

localDataUI(id, panel)

recountDataUI(id, panel)

groupsUI(id, tab)

helpUI(id, tab)
```

## Arguments

- id:

  Character: identifier

- tab:

  Function to process HTML elements

- panel:

  Function to enclose interface

## Value

HTML elements
