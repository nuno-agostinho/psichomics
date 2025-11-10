# Get alternative splicing events from genes or vice-versa

Get alternative splicing events from genes or vice-versa

## Usage

``` r
getSplicingEventFromGenes(genes, ASevents, data = NULL)

getGenesFromSplicingEvents(ASevents, data = NULL)
```

## Arguments

- genes:

  Character: gene symbols (or TCGA-styled gene symbols)

- ASevents:

  Character: alternative splicing events

- data:

  Matrix or data frame: alternative splicing information

## Value

Named character containing alternative splicing events or genes and
their respective genes or alternative splicing events as names
(depending on the function in use)

## Details

A list of alternative splicing events is required to run
`getSplicingEventFromGenes`

## Examples

``` r
ASevents <- c("SE_1_+_201763003_201763300_201763374_201763594_NAV1",
              "SE_1_+_183515472_183516238_183516387_183518343_SMG7",
              "SE_1_+_183441784_183471388_183471526_183481972_SMG7",
              "SE_1_+_181019422_181022709_181022813_181024361_MR1",
              "SE_1_+_181695298_181700311_181700367_181701520_CACNA1E")
genes <- c("NAV1", "SMG7", "MR1", "HELLO")

# Get splicing events from genes
matchedASevents <- getSplicingEventFromGenes(genes, ASevents)

# Names of matched events are the matching input genes
names(matchedASevents)
#> [1] "NAV1" "SMG7" "SMG7" "MR1" 
matchedASevents
#>                                                  NAV1 
#> "SE_1_+_201763003_201763300_201763374_201763594_NAV1" 
#>                                                  SMG7 
#> "SE_1_+_183515472_183516238_183516387_183518343_SMG7" 
#>                                                  SMG7 
#> "SE_1_+_183441784_183471388_183471526_183481972_SMG7" 
#>                                                   MR1 
#>  "SE_1_+_181019422_181022709_181022813_181024361_MR1" 

# Get genes from splicing events
matchedGenes <- getGenesFromSplicingEvents (ASevents)

# Names of matched genes are the matching input alternative splicing events
names(matchedGenes)
#> [1] "SE_1_+_201763003_201763300_201763374_201763594_NAV1"   
#> [2] "SE_1_+_183515472_183516238_183516387_183518343_SMG7"   
#> [3] "SE_1_+_183441784_183471388_183471526_183481972_SMG7"   
#> [4] "SE_1_+_181019422_181022709_181022813_181024361_MR1"    
#> [5] "SE_1_+_181695298_181700311_181700367_181701520_CACNA1E"
matchedGenes
#>    SE_1_+_201763003_201763300_201763374_201763594_NAV1 
#>                                                 "NAV1" 
#>    SE_1_+_183515472_183516238_183516387_183518343_SMG7 
#>                                                 "SMG7" 
#>    SE_1_+_183441784_183471388_183471526_183481972_SMG7 
#>                                                 "SMG7" 
#>     SE_1_+_181019422_181022709_181022813_181024361_MR1 
#>                                                  "MR1" 
#> SE_1_+_181695298_181700311_181700367_181701520_CACNA1E 
#>                                              "CACNA1E" 
```
