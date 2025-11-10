# Preparing an Alternative Splicing Annotation for psichomics

## Creating custom alternative splicing annotation

psichomics quantifies alternative splicing based on alternative splicing
event annotations from [MISO](http://genes.mit.edu/burgelab/miso/),
[SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa),
[VAST-TOOLS](https://github.com/vastgroup/vast-tools) and
[rMATS](http://rnaseq-mats.sourceforge.net). New alternative splicing
annotation may be prepared and used in psichomics by parsing alternative
splicing events from those tools. Please contact me if you would like to
see support for other tools.

This tutorial will guide you on how to parse alternative splicing events
from different tools. To do so, start by loading the following packages:

``` r
library(psichomics)
library(plyr)
```

### SUPPA annotation

[SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa) generates
alternative splicing events based on a transcript annotation. Start by
running SUPPA’s `generateEvents` script with a transcript file (GTF
format) for all event types, if desired. See [SUPPA’s
page](https://bitbucket.org/regulatorygenomicsupf/suppa) for more
information.

The resulting output will include a directory containing tab-delimited
files with alternative splicing events (one file for each event type).
Hand the path of this directory to the function
[`parseSuppaAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md),
prepare the annotation using
[`prepareAnnotationFromEvents()`](https://nuno-agostinho.github.io/psichomics/reference/prepareAnnotationFromEvents.md)
and save the output to a RDS file:

``` r
# suppaOutput <- "path/to/SUPPA/output"

# Replace `genome` for the string with the identifier before the first
# underscore in the filenames of that directory (for instance, if one of your 
# filenames of interest is "hg19_A3.ioe", the string would be "hg19")
suppa <- parseSuppaAnnotation(suppaOutput, genome="hg19")
annot <- prepareAnnotationFromEvents(suppa)

# suppaFile <- "suppa_hg19_annotation.RDS"
saveRDS(annot, file=suppaFile)
```

### rMATS annotation

Just like [SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa),
[rMATS](http://rnaseq-mats.sourceforge.net) also allows to generate
alternative splicing events based on a transcript annotation, although
two BAM or FASTQ files are required to generate alternative splicing
events. Read [rMATS’ page](http://rnaseq-mats.sourceforge.net) for more
information.

The resulting output of [rMATS](http://rnaseq-mats.sourceforge.net) is
then handed out to the function
[`parseMatsAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md):

``` r
# matsOutput <- "path/to/rMATS/output"
mats <- parseMatsAnnotation(
    matsOutput,         # Output directory from rMATS
    genome = "fromGTF", # Identifier of the filenames
    novelEvents=TRUE)   # Parse novel events?
annot <- prepareAnnotationFromEvents(mats)

# matsFile <- "mats_hg19_annotation.RDS"
saveRDS(annot, file=matsFile)
```

### MISO annotation

Simply retrieve [MISO’s alternative splicing
annotation](https://miso.readthedocs.io/en/fastmiso/annotation.html) and
give the path to the downloaded folder as input.

``` r
# misoAnnotation <- "path/to/MISO/annotation"
miso <- parseMisoAnnotation(misoAnnotation)
annot <- prepareAnnotationFromEvents(miso)

# misoFile <- "miso_AS_annotation_hg19.RDS"
saveRDS(annot, file=misoFile)
```

### VAST-TOOLS annotation

Download and extract [VAST-TOOLS’ alternative splicing
annotation](https://github.com/vastgroup/vast-tools#vastdb-libraries)
and use the path to the `TEMPLATES` subfolder as the input of
[`parseVastToolsAnnotation()`](https://nuno-agostinho.github.io/psichomics/reference/parseMisoAnnotation.md).
Complex events (i.e. alternative coordinates for the exon ends) are not
currently supported.

``` r
# vastAnnotation <- "path/to/VASTDB/libs/TEMPLATES"
vast <- parseVastToolsAnnotation(vastAnnotation, genome="Hsa")
annot <- prepareAnnotationFromEvents(vast)

# vastFile <- "vast_AS_annotation_hg19.RDS"
saveRDS(annot, file=vastFile)
```

### Combining annotation from different sources

To combine the annotation from different sources, provide the parsed
annotations of interest simultaneously to the function
`prepareAnnotationFromEvents`:

``` r
# Combine the annotation from SUPPA, MISO, rMATS and VAST-TOOLS
annot <- prepareAnnotationFromEvents(suppa, vast, mats, miso)

# annotFile <- "AS_annotation_hg19.RDS"
saveRDS(annot, file=annotFile)
```

## Quantifying alternative splicing using the custom annotation

The created alternative splicing annotation can be used in psichomics
for alternative splicing quantification. To do so, when using the GUI
version of psichomics, be sure to select the **Load annotation from
file…** option, click the button that appears below and select the
recently created RDS file.

Otherwise, if you are using the CLI version, perform the following
steps:

``` r
annot <- readRDS(annotFile) # "annotFile" is the path to the annotation file
junctionQuant <- readFile("ex_junctionQuant.RDS") # example set

psi <- quantifySplicing(annot, junctionQuant)
```

    ## Using 1 of 60 events (2%) whose junctions are present in junction quantification data...

``` r
psi # may have 0 rows because of the small junction quantification set
```

    ##                                                   Normal 1  Normal 2 Normal 3
    ## SE_1_+_23385660_23385840_23385851_23395032_KDM1A 0.7777778 0.4444444      0.6
    ##                                                  Cancer 1  Cancer 2  Cancer 3
    ## SE_1_+_23385660_23385840_23385851_23395032_KDM1A      0.4 0.4285714 0.5555556

## Feedback

All feedback on the program, documentation and associated material
(including this tutorial) is welcome. Please send any suggestions and
comments to:

> Nuno Saraiva-Agostinho (<nunodanielagostinho@gmail.com>)
>
> [Disease Transcriptomics Lab, Instituto de Medicina Molecular
> (Portugal)](http://imm.medicina.ulisboa.pt/group/distrans/)
