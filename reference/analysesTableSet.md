# Set of functions to render differential analyses (plot and table)

Set of functions to render differential analyses (plot and table)

Set up environment and redirect user to a page based on click
information

## Usage

``` r
analysesTableSet(
  session,
  input,
  output,
  analysesType,
  analysesID,
  getAnalysesData,
  getAnalysesFiltered,
  setAnalysesFiltered,
  getAnalysesSurvival,
  getAnalysesColumns,
  setAnalysesColumns,
  getResetPaging,
  setResetPaging
)

processClickRedirection(click, psi = NULL, survival = FALSE)

analysesPlotSet(
  session,
  input,
  output,
  analysesType,
  analysesID,
  getAnalysesData,
  getAnalysesFiltered,
  getAnalysesSurvival
)
```

## Arguments

- session:

  Shiny session

- input:

  Shiny input

- output:

  Shiny output

- analysesType:

  Character: type of analyses (`GE` or `PSI`)

- analysesID:

  Character: identifier

- getAnalysesData:

  Function: get analyses data

- getAnalysesFiltered:

  Function: get filtered analyses data

- setAnalysesFiltered:

  Function: set filtered analyses data

- getAnalysesSurvival:

  Function: get survival data

- getAnalysesColumns:

  Function: get columns

- setAnalysesColumns:

  Function: set columns

- getResetPaging:

  Function: get toggle of reset paging

- setResetPaging:

  Function: set toggle of reset paging

- click:

  List: click information

- psi:

  Data frame or matrix: alternative splicing quantification

- survival:

  Boolean: redirect to survival page?

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
