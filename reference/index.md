# Package index

## General functions

- [`psichomics()`](https://nuno-agostinho.github.io/psichomics/reference/psichomics.md)
  : Start graphical interface of psichomics
- [`parseSplicingEvent()`](https://nuno-agostinho.github.io/psichomics/reference/parseSplicingEvent.md)
  : Parse alternative splicing event identifier
- [`getSplicingEventData()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventData.md)
  : Get splicing event information for given alternative splicing
  quantification data
- [`getSplicingEventFromGenes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventFromGenes.md)
  [`getGenesFromSplicingEvents()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventFromGenes.md)
  : Get alternative splicing events from genes or vice-versa
- [`plotSplicingEvent()`](https://nuno-agostinho.github.io/psichomics/reference/plotSplicingEvent.md)
  : Plot diagram of alternative splicing events
- [`parseCategoricalGroups()`](https://nuno-agostinho.github.io/psichomics/reference/parseCategoricalGroups.md)
  : Parse categorical columns in a data frame
- [`getDownloadsFolder()`](https://nuno-agostinho.github.io/psichomics/reference/getDownloadsFolder.md)
  : Get the path to the Downloads folder

## TCGA data retrieval

Retrieve TCGA data using [Firebrowse
API](http://firebrowse.org/api-docs/)

- [`isFirebrowseUp()`](https://nuno-agostinho.github.io/psichomics/reference/isFirebrowseUp.md)
  :

  Check if [FireBrowse API](http://firebrowse.org/api-docs/) is running

- [`getTCGAdataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md)
  [`getTCGAdates()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md)
  [`getTCGAcohorts()`](https://nuno-agostinho.github.io/psichomics/reference/getTCGAdataTypes.md)
  : Get available parameters for TCGA data

- [`loadTCGAdata()`](https://nuno-agostinho.github.io/psichomics/reference/loadTCGAdata.md)
  : Download and process TCGA data

- [`parseTCGAsampleTypes()`](https://nuno-agostinho.github.io/psichomics/reference/parseTcgaSampleInfo.md)
  [`parseTCGAsampleInfo()`](https://nuno-agostinho.github.io/psichomics/reference/parseTcgaSampleInfo.md)
  : Parse sample information from TCGA sample identifiers

## GTEx data retrieval

- [`getGtexDataTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexDataTypes.md)
  [`getGtexReleases()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexDataTypes.md)
  : Get GTEx data information
- [`getGtexTissues()`](https://nuno-agostinho.github.io/psichomics/reference/getGtexTissues.md)
  : Get GTEx tissues from given GTEx sample attributes
- [`loadGtexData()`](https://nuno-agostinho.github.io/psichomics/reference/loadGtexData.md)
  : Download and load GTEx data

## SRA data retrieval

- [`loadSRAproject()`](https://nuno-agostinho.github.io/psichomics/reference/loadSRAproject.md)
  :

  Download and load SRA projects via
  [recount2](https://jhubiostatistics.shinyapps.io/recount/)

## User-provided data loading

- [`loadLocalFiles()`](https://nuno-agostinho.github.io/psichomics/reference/loadLocalFiles.md)
  : Load local files
- [`prepareSRAmetadata()`](https://nuno-agostinho.github.io/psichomics/reference/prepareSRAmetadata.md)
  [`prepareJunctionQuant()`](https://nuno-agostinho.github.io/psichomics/reference/prepareSRAmetadata.md)
  [`prepareGeneQuant()`](https://nuno-agostinho.github.io/psichomics/reference/prepareSRAmetadata.md)
  : Prepare user-provided files to be loaded into psichomics

## Gene expression pre-processing

- [`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)
  : Plot row-wise statistics
- [`plotGeneExprPerSample()`](https://nuno-agostinho.github.io/psichomics/reference/plotGeneExprPerSample.md)
  : Plot distribution of gene expression per sample
- [`plotLibrarySize()`](https://nuno-agostinho.github.io/psichomics/reference/plotLibrarySize.md)
  : Plot library size
- [`filterGeneExpr()`](https://nuno-agostinho.github.io/psichomics/reference/filterGeneExpr.md)
  : Filter genes based on their expression
- [`normaliseGeneExpression()`](https://nuno-agostinho.github.io/psichomics/reference/normaliseGeneExpression.md)
  [`normalizeGeneExpression()`](https://nuno-agostinho.github.io/psichomics/reference/normaliseGeneExpression.md)
  : Filter and normalise gene expression
- [`convertGeneIdentifiers()`](https://nuno-agostinho.github.io/psichomics/reference/convertGeneIdentifiers.md)
  : Convert gene identifiers

## PSI quantification and filtering

- [`getSplicingEventTypes()`](https://nuno-agostinho.github.io/psichomics/reference/getSplicingEventTypes.md)
  : Get supported splicing event types

- [`listSplicingAnnotations()`](https://nuno-agostinho.github.io/psichomics/reference/listSplicingAnnotations.md)
  : List alternative splicing annotations

- [`loadAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/loadAnnotation.md)
  :

  Load alternative splicing annotation from `AnnotationHub`

- [`quantifySplicing()`](https://nuno-agostinho.github.io/psichomics/reference/quantifySplicing.md)
  : Quantify alternative splicing events

- [`discardLowCoveragePSIvalues()`](https://nuno-agostinho.github.io/psichomics/reference/discardLowCoveragePSIvalues.md)
  : Remove alternative splicing quantification values based on coverage

- [`filterPSI()`](https://nuno-agostinho.github.io/psichomics/reference/filterPSI.md)
  : Filter alternative splicing quantification

- [`plotRowStats()`](https://nuno-agostinho.github.io/psichomics/reference/plotRowStats.md)
  : Plot row-wise statistics

## Custom alternative splicing annotation

- [`parseSuppaAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md)
  [`parseVastToolsAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md)
  [`parseMisoAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md)
  [`parseMatsAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md)
  : Parse events from alternative splicing annotation
- [`prepareAnnotationFromEvents()`](https://nuno-agostinho.github.io/psichomics/reference/prepareAnnotationFromEvents.md)
  : Prepare annotation from alternative splicing events

## Data Grouping

- [`createGroupByAttribute()`](https://nuno-agostinho.github.io/psichomics/reference/createGroupByAttribute.md)
  : Split elements into groups based on a given column of a dataset

- [`filterGroups()`](https://nuno-agostinho.github.io/psichomics/reference/filterGroups.md)
  : Filter groups with less data points than the threshold

- [`groupPerElem()`](https://nuno-agostinho.github.io/psichomics/reference/groupPerElem.md)
  : Assign one group to each element

- [`getGeneList()`](https://nuno-agostinho.github.io/psichomics/reference/getGeneList.md)
  : Get curated, literature-based gene lists

- [`getSampleFromSubject()`](https://nuno-agostinho.github.io/psichomics/reference/getSampleFromSubject.md)
  : Get samples matching the given subjects

- [`getSubjectFromSample()`](https://nuno-agostinho.github.io/psichomics/reference/getSubjectFromSample.md)
  : Get subjects from given samples

- [`testGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/testGroupIndependence.md)
  : Multiple independence tests between reference groups and list of
  groups

- [`plotGroupIndependence()`](https://nuno-agostinho.github.io/psichomics/reference/plotGroupIndependence.md)
  :

  Plot `-log10(p-values)` of the results obtained after multiple group
  independence testing

## Principal component analysis (PCA)

- [`performPCA()`](https://nuno-agostinho.github.io/psichomics/reference/performPCA.md)
  : Perform principal component analysis after processing missing values
- [`plotPCAvariance()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCAvariance.md)
  : Create the explained variance plot from a PCA
- [`calculateLoadingsContribution()`](https://nuno-agostinho.github.io/psichomics/reference/calculateLoadingsContribution.md)
  : Calculate the contribution of PCA loadings to the selected principal
  components
- [`plotPCA()`](https://nuno-agostinho.github.io/psichomics/reference/plotPCA.md)
  : Create a scatterplot from a PCA object

## Independent component analysis (ICA)

- [`performICA()`](https://nuno-agostinho.github.io/psichomics/reference/performICA.md)
  : Perform independent component analysis after processing missing
  values
- [`plotICA()`](https://nuno-agostinho.github.io/psichomics/reference/plotICA.md)
  : Create multiple scatterplots from ICA

## Differential analyses

- [`diffAnalyses()`](https://nuno-agostinho.github.io/psichomics/reference/diffAnalyses.md)
  : Perform statistical analyses
- [`plotDistribution()`](https://nuno-agostinho.github.io/psichomics/reference/plotDistribution.md)
  : Plot sample distribution

## Gene expression and alternative splicing correlation

- [`correlateGEandAS()`](https://nuno-agostinho.github.io/psichomics/reference/correlateGEandAS.md)
  : Correlate gene expression data against alternative splicing
  quantification
- [`` `[`( ``*`<GEandAScorrelation>`*`)`](https://nuno-agostinho.github.io/psichomics/reference/plot.GEandAScorrelation.md)
  [`plot(`*`<GEandAScorrelation>`*`)`](https://nuno-agostinho.github.io/psichomics/reference/plot.GEandAScorrelation.md)
  [`print(`*`<GEandAScorrelation>`*`)`](https://nuno-agostinho.github.io/psichomics/reference/plot.GEandAScorrelation.md)
  [`as.table(`*`<GEandAScorrelation>`*`)`](https://nuno-agostinho.github.io/psichomics/reference/plot.GEandAScorrelation.md)
  : Display results of correlation analyses

## Survival analysis

- [`getAttributesTime()`](https://nuno-agostinho.github.io/psichomics/reference/getAttributesTime.md)
  : Get time values for given columns in a clinical dataset
- [`processSurvTerms()`](https://nuno-agostinho.github.io/psichomics/reference/processSurvTerms.md)
  : Process survival curves terms to calculate survival curves
- [`survfit(`*`<survTerms>`*`)`](https://nuno-agostinho.github.io/psichomics/reference/survfit.survTerms.md)
  : Create survival curves
- [`survdiffTerms()`](https://nuno-agostinho.github.io/psichomics/reference/survdiffTerms.md)
  : Test Survival Curve Differences
- [`plotSurvivalCurves()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalCurves.md)
  : Plot survival curves
- [`testSurvival()`](https://nuno-agostinho.github.io/psichomics/reference/testSurvival.md)
  : Test the survival difference between groups of subjects
- [`assignValuePerSubject()`](https://nuno-agostinho.github.io/psichomics/reference/assignValuePerSubject.md)
  : Assign average sample values to their corresponding subjects
- [`labelBasedOnCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/labelBasedOnCutoff.md)
  : Label groups based on a given cutoff
- [`optimalSurvivalCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/optimalSurvivalCutoff.md)
  : Calculate optimal data cutoff that best separates survival curves
- [`plotSurvivalPvaluesByCutoff()`](https://nuno-agostinho.github.io/psichomics/reference/plotSurvivalPvaluesByCutoff.md)
  : Plot p-values of survival difference between groups based on
  multiple cutoffs

## Gene, transcript and protein annotation retrieval

- [`queryEnsemblByGene()`](https://nuno-agostinho.github.io/psichomics/reference/queryEnsemblByGene.md)
  [`queryEnsemblByEvent()`](https://nuno-agostinho.github.io/psichomics/reference/queryEnsemblByGene.md)
  : Query information from Ensembl
- [`ensemblToUniprot()`](https://nuno-agostinho.github.io/psichomics/reference/ensemblToUniprot.md)
  : Convert from Ensembl to UniProt identifier
- [`plotProtein()`](https://nuno-agostinho.github.io/psichomics/reference/plotProtein.md)
  : Plot protein features
- [`plotTranscripts()`](https://nuno-agostinho.github.io/psichomics/reference/plotTranscripts.md)
  : Plot transcripts
