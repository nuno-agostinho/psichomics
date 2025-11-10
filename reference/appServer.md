# Server logic

Instructions to build the Shiny app

## Usage

``` r
appServer(input, output, session)

analysesServer(input, output, session)

diffEventServer(ns, input, output, session, psi)

correlationServer(input, output, session)

diffExpressionServer(input, output, session)

diffExpressionEventServer(input, output, session)

diffExpressionTableServer(input, output, session)

diffSplicingServer(input, output, session)

diffSplicingEventServer(input, output, session)

diffSplicingTableServer(input, output, session)

dimReductionServer(input, output, session)

icaServer(input, output, session)

pcaServer(input, output, session)

infoServer(input, output, session)

survivalServer(input, output, session)

templateServer(input, output, session)

dataServer(input, output, session)

firebrowseServer(input, output, session)

geNormalisationFilteringServer(input, output, session)

gtexDataServer(input, output, session)

inclusionLevelsServer(input, output, session)

inclusionLevelsFilterServer(input, output, session)

localDataServer(input, output, session)

recountDataServer(input, output, session)

groupsServer(input, output, session)

helpServer(input, output, session)
```

## Arguments

- input:

  Shiny input

- output:

  Shiny output

- session:

  Shiny session

## Value

`NULL` (function is only used to modify the Shiny session's state or
internal variables)
