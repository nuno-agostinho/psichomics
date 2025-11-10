# Match MISO's splicing event IDs with the IDs present in the alternative splicing annotation file and get events in a data frame

Match MISO's splicing event IDs with the IDs present in the alternative
splicing annotation file and get events in a data frame

## Usage

``` r
parseMisoEventID(eventID, annotation, IDcolumn)
```

## Arguments

- eventID:

  Character: alternative event IDs

- annotation:

  Data.frame: alternative event annotation file

- IDcolumn:

  Integer: index of the column with the event ID's in the alternative
  event annotation file

## Value

Data frame of the matching events (or `NA` when nothing matches)

## Details

For faster execution times, provide a vector of event IDs.

For more information about MISO, see <http://miso.readthedocs.org>.

## Note

If possible, it's recommend to use smaller subsets of the alternative
events' annotation instead of all data for faster runs. For example,
when trying to match only skipped exons event IDs, only use the
annotation of skipped exons instead of using a mega annotation with all
event types.

## Examples

``` r
eventID <- c("114785@uc001sok.1@uc001soj.1", "114784@uc001bxm.1@uc001bxn.1")
# the annotation is one of the GFF3 files needed to run MISO
gff3 <- system.file("extdata", "miso_AS_annot_example.gff3", 
                    package="psichomics")
annotation <- read.delim(gff3, header=FALSE, comment.char="#")
IDcolumn <- 9
psichomics:::parseMisoEventID(eventID, annotation, IDcolumn)
#> [[1]]
#>      V1  V2   V3       V4       V5 V6 V7 V8
#> 1 chr12 AFE gene 57916659 57920171  .  +  .
#> 2 chr12 AFE mRNA 57919131 57920171  .  +  .
#> 3 chr12 AFE exon 57919131 57920171  .  +  .
#> 4 chr12 AFE mRNA 57916659 57918199  .  +  .
#> 5 chr12 AFE exon 57916659 57916794  .  +  .
#> 6 chr12 AFE exon 57917812 57917875  .  +  .
#> 7 chr12 AFE exon 57918063 57918199  .  +  .
#> 
#> [[2]]
#>      V1  V2   V3       V4       V5 V6 V7 V8
#> 13 chr1 AFE gene 34630512 34631443  .  -  .
#> 14 chr1 AFE mRNA 34630512 34630875  .  -  .
#> 15 chr1 AFE exon 34630512 34630875  .  -  .
#> 16 chr1 AFE mRNA 34631348 34631443  .  -  .
#> 17 chr1 AFE exon 34631348 34631443  .  -  .
#> 
```
