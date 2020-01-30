# psichomics 1.12.1 (29 January, 2020)

* Alternative splicing events can now be represented via diagrams:
    - Redesign of alternative splicing event selection (graphical interface)
    - `plotSplicingEvent()` plots diagram representation of alternative
    splicing events
    - In the visual interface, alternative splicing event diagrams were added
    below distribution plots (to quickly illustrate higher and lower values of
    alternative splicing quantification) and in annotation page
* User-provided junction quantification loading:
    - Support junction coordinates from mitochondrial, Z and W chromosomes
    - Fix issues with files containing splice junctions within `random`, `alt`
    and `unknown` chromosomes by discarding those rows (a warning is raised)
* Alternative splicing annotation:
    - `listSplicingAnnotations()` can now be filtered by species, assembly and
    data of available annotations
* Improve import/export of data groups from/to a file, including colour support
* Copy-edit and improve all tutorials, welcome message and help tab
* Include link to article in Methods in Molecular Biology:
https://doi.org/10.1007/978-1-0716-0301-7_10

## More bug fixes and minor changes

* TCGA and SRA data loading:
    - `loadFirebrowseData()` now returns expected data when asking for multiple
    datasets (such as in the case of performing a pan-cancer analysis) in both
    visual and command-line interface
    - SRA projects containing only one column of extra information in sample
    metadata are now correctly loaded instead of raising an error
    (`loadSRAproject()`)
* Data loading and manipulation:
    - Copy-edit information on the format of user-provided files
    - Warn when discarding rows with duplicated rownames after loading
    user-provided files
    - `getGtexDataTypes()` is now exported, as expected
    - `parseSplicingEvent()` now returns the coordinates as numeric if 
    `coords = TRUE` and `char = FALSE`
    - Improve dialog when trying to load a local folder without any supported 
    files available
* Alternative splicing annotation:
    - Confirmation dialog is not displayed any more when creating a folder
    (specially useful while running the visual interface)
    - Allow to select cache directory of `AnnotationHub` (command-line
    interface)
    - Fix `prepareGeneQuant()` discarding the argument `strandedness` if
    either `stranded` or `stranded (reverse)`
* Data grouping:
    - Samples not associated with any subject are now kept when exporting groups
    to a file
    - Fix issues related with importing groups based on a file with groups of
    splicing events and/or genes
    - Inform user when groups are successfully loaded from a file and whether
    any group elements are discarded in the process
    - Show number of genes contained in pre-made gene lists
    - Discard unavailable genes when creating group of genes based on pre-made
    gene lists (unless these are automatically created at startup)
* Density plot (`plotDistribution()`):
    - Fix visual bug when plotting a group with only one sample (if one data
    point is available for a group, only the rug plot is drawn)
    - After hiding all plot series, rug plots of the different groups can be 
    distinguished based on the Y axis (different arbitrary Y values are given to
    each rug plot series)
    - Rug plot labels now show data values if sample names are not provided
    - Rug plot labels can now be rotated (rotation is not enabled by default
    given `Highcharts` issues that may occur at different zoom levels and
    depending on proximity between different sample values)
    - Fix misguiding example in function documentation
    - When hovering the values in the rug plot, the colour of the tooltip is now
    the same used for the rug points as expected
* Survival analysis (p-value plot):
    - Fix alternative splicing quantification cutoff being selected based on the
    one whose difference has the highest (instead of the lowest) p-value
    - Fix plot line label presenting "p < 0.05" independently of the threshold 
    used for significance
* Gene, transcript and protein annotation:
    - Show available annotation information and query PubMed even if Ensembl is
    down
    - Avoid app crash when searching for PubMed articles too many times
* Improve warning/error alerts:
    - Fix alerts crashing the visual interface
    - Improve message formatting regarding HTML-containing alerts
* Improve unit tests and function documentation

# psichomics 1.10.2 (4 Oct, 2019)

Fix warnings in Windows-specific CI tools and Bioconductor

# psichomics 1.10.1 (19 Sep, 2019)

Avoid timeouts in CI tools and Bioconductor builds by skipping a specific unit
test

# psichomics 1.8.3 (21 April, 2019)

Replace deprecated `R.utils::evalWithTimeout()` with `R.utils::withTimeout()`

# psichomics 1.8.2 (26 March, 2019)

* Data loading:
    - GTEx data can now be automatically downloaded and loaded on-demand
    - By default, data table now only displays at most the first 100 columns 
    for performance reasons
* Alternative splicing quantification:
    - Allow to discard samples before alternative splicing quantification
    - In alternative splicing quantification dataset summary, plot
    quantification based on median, variance and range per splicing event across
    samples to provide tools to filter quantification (`plotRowStats()`)
* Gene expression filtering and normalisation:
    - Allow to discard samples before filtering and normalisation
    - Filter low read counts using `edgeR::filterByExpr()`
    - Allow to perform `limma::voom()` on gene expression data (without design
    matrix) and to apply its normalisation methods
    - In gene expression dataset summary, plot distribution of gene expression 
    per sample, distribution of library sizes and gene-wise mean and variance of
    gene expression across samples to provide the user tools to assess gene 
    expression normalisation (`plotGeneExprPerSample()`, `plotDistribution()` 
    and  `plotRowStats()`, respectively)
    - Convert between different gene identifiers (the original identifier is 
    kept in some conditions, read `convertGeneIdentifiers()`); in the visual
    interface, when filtering and normalising gene expression, ENSEMBL 
    identifiers are converted to gene symbols, by default
* Groups:
    - By default, load pre-made lists of genes when loading gene expression or
    loading/performing alternative splicing quantification
    - Added pre-made list of genes that encode for RNA-binding proteins 
    (Sebestyen et al. 2016), useful to postulate about the regulatory role of
    those proteins based on gene expression and PSI correlation analyses
* Correlation analyses:
    - Allow to use groups of genes and alternative splicing events in 
    correlation analyses
    - Plot specific combinations of gene and alternative splicing events 
    (`[.GEandAScorrelation`())
    - Display progress when performing correlation analyses
    - Display correlation results in a table (`as.table()`)
* Survival:
    - Render p-value plot by cutoff in command-line interface 
    (`plotSurvivalPvaluesByCutoff()`)
* Gene, transcript and protein information:
    - Modify keywords used to search for PubMed articles

## Bug fixes and minor changes

* Improve console logging of error and warning alerts
* Fix crash when loading psichomics with test data that is not locally available
(by automatically downloading said data if not found)
* Allow to edit file path in file/directory browser elements
* Documentation:
    - Export functions mentioned in the documentation
    - Hide documentation of internal functions from the PDF reference manual
* Loading SRA data:
    - Accept a vector of files as the first argument (easier to use with 
    `list.files()`)
    - Ask to overwrite file if one exists with the same name as the output file
* Groups:
    - Automatically set dropdown width for group attribute selection
    - Minor improvements to the group creation interface
    - Fix error when creating groups containing only samples and no matching 
    subjects
    - Fix warning when displaying group preview only based on subjects or 
    samples
* Alternative splicing quantification:
    - Automatically set the human genome version after loading data from TCGA
    (hg19), GTEx (hg19) or recount2 (hg38)
    - Fix progress bar
    - Decrease loading time after quantifying alternative splicing
    - By default, quantify skipped exons, mutually exclusive exons, alternative 
    3' and 5' splice sites, and alternative first and last exons; the default
    option is now consistent across both the visual and command-line interfaces)
* Differential analyses:
    - Allow distribution plots to show the name of the samples when hovering or
    when a rug plot is rendered if `rugLabels = TRUE` (function
    `plotDistribution()`)
    - Fix distribution plots requiring all samples in a group when using 
    function `plotDistribution()`
    - Fix group colours and opacity for rug plot points within the distribution 
    plot (occurred in the command-line version; function `plotDistribution()`)
    - Fix wrong information in the table of differential splicing results (only
    occurs when the first splicing event is one for which there is not enough
    information to calculate statistical tests)
    - Fix inconsistency when presenting median and variance differences between
    gene expression and alternative splicing quantification
    - Fix error when groups contain samples outside the data being analysed
* Gene, transcript and protein information:
    - Fix article title formatting (e.g. bold and italics)
* Update psichomics citation with journal publication date

# psichomics 1.8.1 (3 December, 2018)

* Add tag `ImmunoOncology` to BiocViews

# psichomics 1.6.2 (2 October, 2018)

* Update citations to link to article in Nucleic Acids Research:
https://doi.org/10.1093/nar/gky888
* Copy-edit README, tutorials, DESCRIPTION, NEWS and code
* Update screen shot
* Add metadata regarding loaded datasets (for instance, path and name of files 
used for dataset loading, options set upon dataset creation, etc.)
* Differential splicing and expression analysis:
    - Allow to input formulas when highlighting points; this is an improvement
    when plotting -log10(q-values), for instance, where it is now possible to
    highlight values >= -log10(0.05), i.e. q-values <= 0.5, instead of inputting
    an approximate value
* Correlation analyses:
    - Colour samples by group
    - Render contour based on a density estimate
    - Display exportable table with correlation analyses
* Groups:
    - Show group colour in group selection element
* Retrieve articles from PubMed instead of PMC

## Bug fixes and minor changes

* Improve file browser input in remote servers:
    - Allow to directly input a path without using the file browser dialog
    - Display warning when calling the file browser dialog in RStudio Server or
    in an unsupported system
* TCGA data:
    - Inform of sequencing technology used to obtain gene expression data
* GTEx data:
    - Fix issues in downstream analyses when using GTEx data subset by tissues
* Differential analysis:
    - Avoid limiting X axis of density plot when values outside 0 and 1 are
    provided, as in the case of gene expression (`plotDistribution()` function)
    - Fix filenames of exported tables from differential expression analysis
    incorrectly mentioning differential splicing analysis instead
* Correlation analysis:
    - When plotting correlation analysis, loess curve was performed based on
    gaussian fitting, independent of the "family" argument of
    `plotCorrelation()`
* Improve console logging of error and warning messages

# psichomics 1.6.1 (5 July, 2018)

* Improve support for analysing SRA and user-provided data:
	- New tutorial on performing alternative splicing analysis using SRA and 
	user-provided data from FASTQ files (other tutorials were also updated)
	- Process gene and splice junction counts from STAR output for subsequent
	loading into psichomics
	- Automatically download and load data from select SRA projects using the
	recount R package
* Data grouping:
    - Preview groups based on selected attribute before creation
* Gene, transcript and protein information:
    - Display alternative splicing events in transcript plot
* Correlation between gene expression and alternative splicing:
    - Allow to change plot height and font size

## Bug fixes and minor changes

* Support new version of Human (hg38 assembly) alternative splicing annotation 
(fixes wrong coordinates for many minus-strand splicing events)
* Fix issues with TCGA sample metadata:
    - No sample metadata loaded when loading TCGA files from local folder
    - Not retrieving sample metadata for all TCGA samples when multiple junction
    quantification files are loaded
* Data grouping:
    - Correctly suggest how to select variables in subset expression and give 
    more examples on how to use this feature
    - Fix subject/sample attributes not being updated for selection for the 
    `GREP` interface
    - Fix error alerts not appearing when trying to create groups based on 
    invalid parameters
    - Add an error alert when trying to extract group data from an user-provided
    file without such information
* Alternative splicing quantification:
    - Fix crash when quantifying splicing for a single sample
* Differential splicing analysis (exploratory):
    - Fix crash when inputting non-numeric values or minimum larger than the
    maximum limit as the axes for the splicing scatter plot

# psichomics 1.4.5 (4 Apr, 2018)

* Mention psichomics manuscript throughout psichomics
* Copy-edit graphical user interface and respective tutorial
* Fix warnings and errors in Bioconductor

# psichomics 1.4.4 (12 Feb, 2018)

* Update CITATION file to show citation to article in bioRxiv:
https://www.biorxiv.org/content/early/2018/02/07/261180
* Update vignettes to include a case study based on the aforementioned article
* Gene expression pipeline:
    - Perform `limma-trend` by default
* Alternative splicing quantification:
    - Quantification of alternative first and last exons: following more 
    thorough testing, the new exon-centred method was considered to be less 
    relevant to exploratory analysis (specially when compared with other types 
    of events); as such, both methods are now available for quantification
* Dimensionality reduction:
    - Change tolerated missing values per event to 5% by default for both PCA
    and ICA (in both visual and command-line interfaces)
    - PCA: In the table that shows the events that most contribute to the 
    selected principal components, show the rank
* Groups:
    - Automatically set sample type groups (i.e. normal, primary solid tumour,
    metastatic, etc.) for TCGA samples (visual interface only)
* Differential analysis:
    - Use numeric fields instead of sliders to precisely filter data
    - Fix table for differential expression not being filtered based on
    highlighted genes
    - Filter splicing events or genes to use when performing exploratory
    differential analyses
* Survival analysis:
    - Allow to stratify patients based on optimal gene expression cutoff
    - Select samples to be used for survival analysis

## Bug fixes

* Fix problems related with `DT` versions >= 0.3:
    - Groups displayed as having the same attributes as last created group
    - Table not being updated in differential analyses according to event
    filtering based on the volcano plot

# psichomics 1.4.3 (12 Jan, 2018)

* Alternative splicing quantification:
    - Improve speed and memory usage when dealing with larger datasets
    - Improve quantification of alternative first and last exons: quantify 
    alternative first and last exons based on all exon-exon junction reads that 
    support each of the alternative exons
    - Print progress bar in R console
* Principal component analysis:
    - Change tolerated missing values per event to 5% by default
    - Show/save table with the contribution of events (for alternative splicing 
    quantification) or genes (for gene expression) to the selected principal
    components
    - Add extra information when hovering variables in loading plot
    - Allow to plot top 100 variables that most contribute to the selected 
    principal components (faster rendering of and interaction with loading 
    plots)
* Differential analysis:
    - Allow to input a list of groups for the "group" argument of the functions
    `diffAnalyses()` and `plotDistribution()` (command-line interface)
    - Allow to input a non-numeric vector or a row of a matrix/data frame in the
    "data" argument of the function `plotDistribution()` (command-line
    interface)
* Correlation between gene expression and alternative splicing:
    - Perform correlation between gene expression of multiple genes and 
    quantification of multiple alternative splicing events

## Bug fixes

* Fix unnamed events when only one event for a event type is returned
* Minor copy-editing and overall improvements

# psichomics 1.4.2 (19 Dec, 2017)

* Fix error when trying to load alternative splicing annotation (given updated
hg19 and hg38 annotation that is now available for use with psichomics)

# psichomics 1.4.1 (14 Dec, 2017)

* GTEx data loading:
    - Add input elements to allow GTEx gene expression loading in the graphical
    interface
    - Fix bug that did not allow to select tissues to load GTEx v7 data
    (graphical interface)
    - Fix splicing events not being quantified based on GTEx v7 junction reads
* Gene expression normalisation:
    - Fix misleading gene expression (non-)normalisation by converting reads to 
    counts per million (CPM) using `edgeR::cpm()` after normalisation using
    `edgeR::calcNormFactor()`
* Alternative splicing quantification:
    - Updated support to properly parse new notation of alternative splicing 
    annotation from Bioconductor (backwards compatible with older notation)
    - Raise error when no splicing events after quantification
    - Fix warning following the quantification of splicing events or its loading
    (incorrect parsing of gene information from splicing events)
* Dimensionality reduction:
    - Use the number (instead of the percentage) of tolerated missing values per
    sample as the argument to impute data from the remaining samples for those
    values before performing dimensionality reduction; by default, missing 
    values are tolerated for 10 samples
* Update file description and README
* Minor bug fixes and improvements

# psichomics 1.4.0 (22 Oct, 2017)

* Support for loading new GTEx V7 data
* Support gene expression data:
    - Load, filter, normalise and perform log2-transformation on gene expression
    data from TCGA
    - Perform principal component analysis based on gene expression data, 
    survival analysis by gene expression cutoff and pairwise differential 
    gene expression analysis
    - Correlate gene expression of a given gene against PSI values of multiple
    alternative splicing events
* Data loading:
    - Add step-wise instructions about loading of user-provided files
    - Filter GTEx junction quantification based on tissues of interest (all 
    tissues are loaded by default)
    - Quantify splicing based on a list of genes (splicing events within all 
    genes are quantified by default)
    - Parse sample information from TCGA samples using `parseTcgaSampleInfo()`
    - Generate TCGA sample metadata when loading TCGA junction quantification
    - Present data summary after loading the data
* Data grouping:
    - Redesigned group creation and selection
    - Create groups based on genes and alternative splicing events
    - Assign a customisable colour per data group
    - Export or import patient and sample identifiers of data groups
    - Add new set operations when grouping (such as complement, subtraction and
    symmetric difference)
    - Suggest attributes of interest when creating groups
    - Allow to retrieve the universe of patient and sample identifiers by
    performing the complement group without any group selected
    - Statistically analyse group independence (useful to assess the overlap 
    between a PCA cluster and groups derived from clinical and sample 
    attributes, for instance)
* Differential analysis:
    - Label points based on top differentially spliced events or genes, selected 
    alternative splicing events and/or selected genes
    - Create AS event and gene groups based on filtered or selected AS events 
    and genes in the tables
* Dimensionality reduction techniques:
    - Subset data based on groups of AS events and genes before performing 
    dimensionality reduction
    - Create data groups based on the partitioning clustering of PCA scores
    - Perform independent component analysis (ICA) on alternative splicing 
    quantification and gene expression data
* Survival analysis:
    - Add p-value plot to visually infer the significance of survival analyses
    based on multiple alternative splicing quantification cutoffs
* Gene, transcript and protein information:
    - Information retrieval is now only dependent on a user-defined gene,
    instead of requiring alternative splicing quantification data to be loaded
    
## Bug fixes and other improvements

* Show progress bar when running in the command-line interface
* Fix inconsistent browser history navigation
* Updated the CLI vignette with information on analysing gene expression data
and a quick reference for functions
* Update minimum version required of shiny (1.0.3)
* Avoid replacing selected groups when manipulating new ones
* Differential splicing analysis:
    - Fix data not being rendered in the table when zooming in the plot after 
    data transformation was applied
    - Return p-value of NA instead of 0 when the value of Fligner-Killeen's Test
    for Homogeneity of Variance is infinite
    - Discard value transformations that may return invalid data for the values
    chosen for the X and Y axes
    - Fix point that remains highlighted in the plot after deselecting the only
    selected row of the table
    - Improve readability of plot's tooltip
    - Improve survival curves based on the optimal alternative splicing
    quantification cutoff:
        - Include the survival curve previews in 3 new columns within the
        differential splicing analyses table, instead of below that table; those
        columns consist of the survival curves, the optimal PSI cutoff and the 
        respective log-rank test's p-value
        - Allow to use survival data when plotting and table sorting
        - Include the optimal PSI cutoff and the respective log-rank test's 
        p-value in exported tables
        - Fix link to survival analyses using the previously calculated PSI
        cutoff
* Principal component analysis:
    - When clicking on a alternative splicing event in the loadings plot, the
    appropriate differential splicing analyses will now be automatically
    rendered with the respective options, as expected
* Survival analysis:
    - Properly set the title of survival curves based on the selected splicing
    event's quantification
    - Improve readability of Cox PH models
    - When performing survival analyses by alternative splicing cutoff, each
    patient is assigned the PSI value from the respective sample; for patients
    with more than one sample, the assigned sample is chosen based on the most
    frequent sample type across all patients (before, the first matched 
    non-normal or non-control samples were used)
* Multiple other bug fixes and visual improvements

# psichomics 1.2.1 (24 Apr, 2017)

* Gene, protein and transcript information:
    - Fix missing file required for transcript plots
* Update command-line interface tutorial to render a transcript plot

# psichomics 1.2.0 (22 Apr, 2017)

* Gene, protein and transcript information:
    - Fix tooltip text presentation in transcript plot
    - Fix JavaScript issues when zooming the transcript plot
    - Fix error when plotting events associated with multiple genes
    - Fix error when plotting single-exon transcripts
    - Protein name, length and function are now presented when available
    - Improved general presentation of the information
* Differential splicing analyses:
    - Click and drag in the plot to zoom in and subsequently filter events
    shown in the table
    - Decreased step of sliders
    - Improve interface of previewed survival curves
    - When clicking on a table link to navigate to differential splicing
    analyses of a single event, the appropriate analyses will now be
    automatically rendered with the respective options, as expected
* Settings (renamed to "Help"):
    - Add links to tutorials and user feedback
    - Add app information and acknowledgements
    - Remove unused option for choosing cores (all performed operations are 
    still single-core, given the difficulty of working with multi-processes in
    Shiny)
* Improve dialogues regarding missing data and other minor interface elements
* Update documentation with volcano plot

# psichomics 1.0.9 (10 Apr, 2017)

* Differential splicing analyses:
    - Add volcano plot to represent events through selected attributes, such as
    p-values and descriptive statistics (e.g. median and variance) between 
    groups of interest
    - Transform values of the X and Y axis in the plot using log transformed, 
    inverted and absolute values, for instance
    - Highlight events in the plot based on values of the X and Y axis
    - Table of differential analyses per alternative splicing event is filtered
    according to highlighted and selected events in the plot
* Gene, protein and transcript information:
    - Transcript plot is now interactive and zoomable
    - Protein are now rendered based on selected transcript alone
    - Faster parsing of UniProt's web API response
    - Improve display of article information when data is missing
* Principal component analysis:
    - Improve presentation of available options
* When clicking on previews of differential splicing and survival analyses, the
appropriate analyses will now be automatically rendered with the respective
options
* Fix buggy browser history when the user is directed to a different tab
* Consistently use Firebrowse and Firehose across the package
* Update documentation

# psichomics 1.0.8 (21 Feb, 2017)

* Support GTEx data loading and analysis
* Fix clinical data dependency:
    - Fix error when trying to load a file containing alternative splicing
    quantification without first loading clinical data
    - Fix error where samples from junction quantification were matched to
    clinical information even if clinical data were not loaded
    - Inform user when clinical data is not loaded while trying to plot survival
    curves
* Improve data grouping:
    - Create sample groups like patient groups and perform set operations
    between any created groups
    - Create groups using patient and sample identifiers
    - Check number of patients and samples per group
    - Rename selected groups
    - Alert user when groups cannot be created due to missing data
* Differential splicing analysis:
    - Analyse all samples as one group
* Survival analysis:
    - Select any clinical attribute for starting/follow up and ending times
* Create table containing TCGA sample metadata when calculating or loading
alternative splicing quantification
* Minor UI improvements

# psichomics 1.0.7 (22 Jan, 2017)

* Survival analysis:
    - Fix error caused by some non-matched patients not being in the 
    patient-sample matching matrix

# psichomics 1.0.6 (17 Jan, 2017)

* Update tutorials with more relevant and complex examples
* Update minimum versions required of `highcharter` (0.5.0) and shiny (1.0.0):
    - Fix function usage as according to new version of `highcharter`
    - More options available when exporting plots (PNG, JPEG, SVG, XLS and CSV)
* Faster alternative splicing quantification
* Differential splicing analysis:
    - Fix major bug where samples could be placed in the wrong groups
    - Shorten speed of the calculation for the optimal PSI cutoff that minimises
    the survival difference
    - Fix not performing statistical tests for two selected sample types while
    analysing a single event with three or more sample types
    - Fix differential analysis on one splicing event not working when using
    `diffAnalyses()`
    - Fix differential analysis not showing for individual events before
    navigating to the page where the analysis is performed for all events
    - Improve readability and information of statistical tests for single events
* Principal component analysis:
    - Shorten time taken to calculate principal components and to render the
    loadings plot
    - Fix loadings plot error when rendering some principal components
* Survival analysis:
    - Fix incorrect number of patients from the survival groups in the 
    contextual information for the selected cutoff (below the slider)
    - Improve how alternative splicing quantification is assigned to patients
    based on their samples
* Protein annotation:
    - Warn user when trying to render proteins with no domains

# psichomics 1.0.5 (7 Jan, 2017)

* Navigate history using the browser forward and back buttons
* Fix delay when displaying large data by removing columns containing missing
values exclusively
* Principal component analysis:
    - Improve speed when calculating total contribution of each variable to the
    principal components
* Survival analysis:
    - Shorten calculation of optimal PSI that minimises the survival difference
    - Improve visual cues of optimal PSI cutoff and present p-value of selected
    PSI cutoff
    - Fix ambiguous error messages
    - Fix incorrect Cox model results for formula-based calculations
    - Fix null Cox models crashing the program
* Differential splicing analysis:
    - Select sample types for differential splicing analysis
    - Fix statistical tests not displaying for individual events after
    differentially analysing all events using the other statistical tests

# psichomics 1.0.4 (18 Dec, 2016)

* Correctly load files and quantify alternative splicing for PRAD, OV and PAAD
tumour types from The Cancer Genome Atlas (TCGA)
* Fix session disconnecting when exporting plots in Firefox
* Improve text and behaviour of fields to select datasets and AS events
* Fix author names and add contributor

# psichomics 1.0.3 (13 Dec, 2016)
    
* Bug fixes regarding gene annotation:
    - Fix disabled gene selection when choosing a splicing event associated with
    a single gene after selecting an event related to multiple genes
    - Fix display of PubMed articles related to previously selected gene when
    selecting a single-gene-associated event after selecting an event related to
    multiple genes
* Bug fixes regarding groups:
    - Fix groups by rows not working
    - Fix group selection not working when only one group exists
    - Improve argument name of `getGroupsFrom()`
* Other minor improvements

# psichomics 1.0.2 (3 Dec, 2016)

* Fix UTF-8 encoding in author list

# psichomics 1.0.1 (1 Dec, 2016)

* Improve metadata (title, description, authors and vignette titles)

# psichomics 1.0.0 (5 Oct, 2016)
    
* First release of psichomics
* Quantify, analyse and visualise alternative splicing data
