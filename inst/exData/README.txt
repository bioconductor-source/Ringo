August 2007		                                   Joern Toedling

The example data in this directory are some snippets of the demo ChIP-chip
data that NimbleGen provide on their FTP server for free download. 
ftp://ftp.nimblegen.com/pub/demo_data/

Of the original files, only those parts concerning reporters/probes
matching the forward strand (Container: FORWARD1) of chromosome 9 
are provided here to keep these example files reasonably small.

* spottypes.txt
This file has not been provided by NimbleGen but was created by hand. 
Basically, its job is to tell the Bioconductor
functions what kind of reporters(probes) on your array are of biological
interest and which ones are controls.
The file says that only reporters whose type (indicated by
GENE_EXPR_OPTION in the pair files, or as CONTAINER in the NimbleGen
design file *.NDF) starts with FORWARD or REVERSE (i.e. FORWARD1,
FORWARD2, REVERSE4 etc.) are considered as interesting "Probe" later on,
while the other types are different forms of controls. 
Have a look at the start of the example pair
file MOD_20551_PMT1_pair.txt
The first line is a comment, contains details on the hybridization and
is ignored by Ringo. The second one are the column headers. Notice that
the second column "GENE_EXPR_OPTION" holds the reporter type. In these
shortened example files there are only probes on interest included.
Normal pair files include lines which say , e.g., "NGS_CONTROLS" or
"H_CODE" here. With the supplied spottypes.txt file these lines are
interpreted as "Negative" or "H_Code" type reporters then.
The file spottypes.txt is used by the limma function "controlStatus"
when reading in the pair files. 
In many cases, you can probably use this spottypes file for your pair files
as well. In some cases, however, you may need to edit it to account for
additional reporter types, kinds of GENE_EXPR_OPTION entries in your files,
that are not covered by the spottypes.txt file yet.

These snippets are used as an example data within Ringo with permission
of NimbleGen Systems, Inc.
