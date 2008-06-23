##--------------------------------------------------------------------------------

## Note: This script has been adopted from the 'davidTiling' package.
## However, since ChIP-chip is not strand-specific,
##  the strand details were omitted.

## a. path to the result table after running the condense*Output.pl script

## FORMAT of that position file should be in NimbleGen's pos format:
# here's an example of their typical head (comment char and one whitespace at start of each line have been added only just here)

# SEQ_ID   CHROMOSOME     PROBE_ID        POSITION        LENGTH  COUNT
# chr1:4483138-4484138    CHR01   CHR01P004483140 4483140 50      1
# chr1:4483138-4484138    CHR01   CHR01P004483245 4483245 50      1
# ...


### THE FUNCTION posToProbeAnno HAS SUPERCEDED THIS SCRIPT


##--------------------------------------------------
## gff

 ## see script 'retrieveGenomicFeatureAnnotation.R'

