#!/usr/bin/perl -w

# script to condense the MUMmer output files per chromosome into 1 file:
# files are given as argument:
# for example: condenseMUMmerOutput.pl chr*out 

# where to write the condensed stuff?
$outfile = "allChrom_Output.txt";
$discardIncomplete = 1;

open(OUT,">$outfile") || die "can't open output file: $!";
# start output:
print OUT "Probe\tChromosome\tStrand\tStart\tEnd\n";


foreach $file (@ARGV) {  # @ARGV array of arguments: specify all files with wildcards

    open(DATA, $file) || die "can't open data file: $file $!";

    while (defined($line=<DATA>)){
        chomp($line); #remove trailing newline
        if ($line =~ /^>.+/ ){ # it's again an ID, so no hit info
	    $previousLine = $line;
	    next;
	} else { # it's not an ID, but hit information
            $previousLine =~ s/^>\s*//; # remove leading '>' and spaces
 	    @previousFields = split(/\s+/,$previousLine); # split one any whitespaces
	    $probeID = $previousFields[0];
            $direction = $previousFields[1]; # stated as second column if there
            if ($direction eq 'Reverse'){
		$strand = "-";
	    } else {
		$strand = "+";
	    }
            $queryLen = $previousFields[$#previousFields];
            $line =~ s/^\s+//; # remove leading spaces
            $line =~ s/\s+$//; # remove trailing spaces
  	    @fields = split(/\s+/,$line); # split one any whitespaces
	    $chrom = $fields[0];  # chromosome number
            $start = $fields[1];  # genomic start site
            $length = $fields[3]; # what is third number?
            # discard incomplete matches:
            if ($discardIncomplete eq 1){
		if ($length < $queryLen){ next;}}
	    $end = $start+$length-1;
	    print OUT "$probeID\t$chrom\t$strand\t$start\t$end\n";
	} # else
    } #while there are still lines in the file

    close (DATA);

} #foreach file

close (OUT); # stop output
