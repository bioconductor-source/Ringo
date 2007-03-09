#!/usr/bin/perl -w

# script to condense the Exonerate output files per chromosome into 1 file:
# files are given as argument:
# for example: condenseExonerateOutput.pl *out

# where to write the condensed stuff?
my $outfile = "allChromExonerateOut.txt";
my $matchlength = 23;
# minimum number of nucleotides in the probe that should match, provide a 
#  second chance (setting the Exonerate score in the programm call being the
#  first) to only obtain "good" matches of probes to the genome

open(OUT,">$outfile") || die "can't open output file: $!";
# start output:
print OUT "SEQ_ID\tPROBE_ID\tCHROMOSOME\tPOSITION\tLENGTH\tMISMATCHES\n";
# position is always the most 5' position of the match on the Plus Strand 
## our ChIP enrichment is not strand-specific

foreach $file (@ARGV) {  # @ARGV array of arguments: specify all files with wildcards

    open(DATA, $file) || die "can't open data file: $file $!";

    for ($j=0; $j<2; $j++){<DATA>} # skip head of out file, 2 lines
	
    #my $chrom = $file;
    #$chrom =~ s/out.psl//;   # remove rest of file name and just keep chromosome number

    while (defined($line=<DATA>)){
        chomp($line); #remove trailing newline
	#$line =~ s/^\s+//; # remove leading spaces
	#$line =~ s/\s+$//; # remove trailing spaces
	if ($line =~ /^--.completed.*/){last;}
	@fields = split(/\t/,$line); # split one tabs
	$match = $fields[8];  # number of matching bases
        $probelen = $fields[2]; # length of probe
	$mismatch = $probelen - $match;  # 
	if (($match < $matchlength) || ($mismatch>1)) 
	{next;} #skip imperfect matches
	$strand = $fields[3]; # 
	$probe = $fields[1]; # probe id
	$chrom = $fields[4]; # 
	$chrom =~ s/(chr)//; # remove trailing 'chr'
	$hitstart = $fields[5] + 1; # +1 since the query match always starts at position 0
	$hitend = $fields[6]; # 
        $hitlength = $hitend - $hitstart + 1;
	$seqID = "chr".$chrom.":".$hitstart."-".$hitend;
	print OUT "$seqID\t$probe\t$chrom\t$hitstart\t$hitlength\t$mismatch\n";
    } #while there are still lines in the file
    
    close (DATA);
    
} #foreach file

close (OUT); # stop output
