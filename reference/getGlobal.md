# Get or set globally accessible elements

Get or set globally accessible elements

## Usage

``` r
getGlobal(category = getCategory(), ..., sep = "_")

setGlobal(category = getCategory(), ..., value, sep = "_")

setData(data)

setDataTable(name, value, category = getCategory())

getAutoNavigation()

setAutoNavigation(auto)

getCores()

setCores(integer)

getSignificant()

setSignificant(integer)

getPrecision()

setPrecision(integer)

getASevents()

getAnnotationHub()

setAnnotationHub(ah)

getASevent()

setASevent(event, data = NULL)

getEvent()

setEvent(event, data = NULL)

getGenes()

getCategories()

getCategory()

setCategory(category)

getCategoryData()

getActiveDataset()

setActiveDataset(dataset)

getClinicalData(attrs = NULL)

getSubjectId()

getSubjectAttributes()

getSampleInfo()

setSampleInfo(value, category = getCategory())

getSampleId()

getSampleAttributes()

getJunctionQuantification(category = getCategory())

getGeneExpression(item = NULL, category = getCategory(), EList = FALSE)

setNormalisedGeneExpression(geneExpr, category = getCategory())

getInclusionLevels()

setInclusionLevels(incLevels, category = getCategory())

getInclusionLevelsSummaryStatsCache(category = getCategory())

setInclusionLevelsSummaryStatsCache(cache, category = getCategory())

getPCA(category = getCategory())

setPCA(pca, category = getCategory())

getICA(category = getCategory())

setICA(ica, category = getCategory())

getCorrelation(category = getCategory())

setCorrelation(correlation, category = getCategory())

getGroupIndependenceTesting(category = getCategory())

setGroupIndependenceTesting(groupIndependenceTesting, category = getCategory())

getSpecies(category = getCategory())

setSpecies(species, category = getCategory())

getAssemblyVersion(category = getCategory())

setAssemblyVersion(assembly, category = getCategory())

getAnnotationName(category = getCategory())

setAnnotationName(annotName, category = getCategory())

getURLtoDownload()

setURLtoDownload(url)
```

## Arguments

- category:

  Character: data category

- ...:

  Arguments to identify a variable

- sep:

  Character to separate identifiers

- value:

  Value to attribute to an element

- data:

  Matrix or data frame: alternative splicing information

- name:

  Character: data table name

- auto:

  Boolean: enable automatic navigation of browser history?

- integer:

  Integer: value of the setting

- ah:

  AnnotationHub

- event:

  Character: alternative splicing event

- dataset:

  Character: dataset name

- attrs:

  Character: name of attributes to retrieve (if `NULL`, the whole
  dataset is returned)

- item:

  Character: name of specific item to retrieve (if `NULL`, the whole
  list is returned)

- EList:

  Boolean: return gene expression datasets as `EList` if possible or as
  data frames?

- geneExpr:

  Data frame or matrix: normalised gene expression

- incLevels:

  Data frame or matrix: inclusion levels

- cache:

  List of summary statistics

- pca:

  `prcomp` object (principal component analysis)

- ica:

  Object containing independent component analysis

- correlation:

  `prcomp` object (correlation analyses)

- groupIndependenceTesting:

  Object containing group independence testing results

- species:

  Character: species

- assembly:

  Character: assembly version

- annotName:

  Character: annotation name

- url:

  Character: URL links to download

## Value

Getters return globally accessible data, whereas setters return `NULL`
as they are only used to modify the Shiny session's state

## Note

Needs to be called inside a reactive function

## See also

Other functions to get and set global variables:
[`getClinicalMatchFrom()`](https://nuno-agostinho.github.io/psichomics/reference/getClinicalMatchFrom.md),
[`getDifferentialExpression()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialExpression.md),
[`getDifferentialSplicing()`](https://nuno-agostinho.github.io/psichomics/reference/getDifferentialSplicing.md),
[`getGroups()`](https://nuno-agostinho.github.io/psichomics/reference/getGroups.md),
[`getHighlightedPoints()`](https://nuno-agostinho.github.io/psichomics/reference/getHighlightedPoints.md),
[`getSelectedDataPanel()`](https://nuno-agostinho.github.io/psichomics/reference/getSelectedDataPanel.md)
