#!/usr/bin/perl -w

# script to condense the BLAT output files per chromosome into 1 file:
# files are given as argument:
# for example: condenseBlatOutput.pl *psl

# where to write the condensed stuff?
my $outfile = "allChromBlatOut.txt";
my $matchlength = 24;

open(OUT,">$outfile") || die "can't open output file: $!";
# start output:
print OUT "Probe\tChromosome\tStrand\tStart\tEnd\tMismatches\n";

foreach $file (@ARGV) {  # @ARGV array of arguments: specify all files with wildcards

    open(DATA, $file) || die "can't open data file: $file $!";

    for ($j=0; $j<5; $j++){<DATA>} # skip head of out file, 5 lines
	
    #my $chrom = $file;
    #$chrom =~ s/out.psl//;   # remove rest of file name and just keep chromosome number

    while (defined($line=<DATA>)){
        chomp($line); #remove trailing newline
	#$line =~ s/^\s+//; # remove leading spaces
	#$line =~ s/\s+$//; # remove trailing spaces
	@fields = split(/\t/,$line); # split one tabs
	$match = $fields[0];  # chromosome number
	$mismatch = $fields[1];  # genomic start site
	if (($match < $matchlength) || ($mismatch>1)) 
	{next;} #skip imperfect matches
	$strand = $fields[8]; #
	$probe = $fields[9]; # 
	$chrom = $fields[13]; # 
	$chrom =~ s/(chr)//; # remove trailing 'chr'
	$hitstart = $fields[15] + 1; # +1 since the query match always starts at position 0
	$hitend = $fields[16]; # 
	print OUT "$probe\t$chrom\t$strand\t$hitstart\t$hitend\t$mismatch\n";
    } #while there are still lines in the file
    
    close (DATA);
    
} #foreach file

close (OUT); # stop output
