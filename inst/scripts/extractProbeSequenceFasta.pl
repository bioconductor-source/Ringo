#!/usr/bin/perl -w

# February 2006                           Joern Toedling

# script to extract the probe sequence from Nimblegen '*.ndf' chip design files

# where to write the extracted Probes
$outfile = "probeSequences.fa";
open(OUT,">$outfile") || die "can't open output file: $!";

my %seenProbeIDs;

foreach $file (@ARGV) {  # @ARGV array of arguments: specify all files with wildcards

    open(DATA, $file) || die "can't open data file: $file $!";

    $line = <DATA>; # skip first line
    chomp($line);
    @headers = split(/\t/,$line);
    print "Number of columns in file: $#headers \n";

    while (defined($line=<DATA>)){

        chomp($line); #remove trailing newline
	#$line =~ s/^\s+//; # remove leading spaces
	#$line =~ s/\s+$//; # remove trailing spaces
	@fields = split(/\t/,$line); # split one tab
	$probeID = $fields[12];  # unique Nimblegen Probe ID
	$assignID = $fields[4]; # ID the probe is assigned to (specified by person ordering the chip, not unique
	$probeSequence = $fields[5]; # base sequence of probe
        if (exists($seenProbeIDs{$probeID})){next;} #skip already written ProbeIDs
        else {
	    $seenProbeIDs{$probeID}=1;
	    print OUT ">$probeID\n$probeSequence\n";
	}
	
    } #while there are still lines in the file

    close (DATA);

} #foreach file

close (OUT); # stop output
