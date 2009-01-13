#!/usr/bin/perl -w

# script to condense the Exonerate output files per chromosome into 1 file:
# files are given as argument:
# for example: condenseExonerateOutput.pl *out

## 13 January 2009:
## changes to handle new Exonerate (version 2.2.0) output, in which the
## match start coordinate can be higher than the end coordinate if the
## match is on the minus strand (previously in this case the output file
## would contain the indication that the reverse complement of the query
## matches at these positions but for each match the start would be
## smaller than the end coordinate of the match.

# where to write the condensed stuff?
my $outfile = "allChromExonerateOut.txt";
my $matchlength = 23;
# minimum number of nucleotides in the probe that should match, provide a 
#  second chance (setting the Exonerate score in the programm call being the
#  first) to only obtain "good" matches of probes to the genome

## other variables:
my ($line, $file, $match, $probelen, $mismatch);
my ($querystrand, $probe, $chrom);
my ($hitstart, $hitend, $hitstrand, $seqID);
my @fields;

open(OUT,">$outfile") || die "can't open output file: $!";
# start output:
print OUT "SEQ_ID\tPROBE_ID\tCHROMOSOME\tPOSITION\tLENGTH\tSTRAND\tMISMATCHES\n";
# position is always the most 5' position of the match on the Plus Strand 
## ChIP-chip is not strand-specific

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
	$querystrand = $fields[3]; # 
	$probe = $fields[1]; # probe id
	$chrom = $fields[4]; # 
	$chrom =~ s/(chr)//; # remove trailing 'chr'
	$hitstart = $fields[5]; 
	$hitend = $fields[6]; # 
        # coordinates are reported somewhat funny by Exonerate, 
        # as positions between bases starting at 0, see this example from the
        # Exonerate web page:
        # sequence:          T A G A C
        # position:         0 1 2 3 4 5
        # a match to all 5 bases would be reported to start at 0 and end at 5
        $hitstrand = $fields[7];
	# if match is to the reverse complement, the end will be smaller
        #  than the start, which will lead to negative length matches
        if ($hitstrand eq "-"){
	    my $temppos = $hitend;
            $hitend = $hitstart;
            $hitstart = $temppos +1;
	} else {
	    $hitstart = $hitstart+1;
	}
        $hitlength = $hitend - $hitstart + 1;
	$seqID = "chr".$chrom.":".$hitstart."-".$hitend;
	print OUT "$seqID\t$probe\t$chrom\t$hitstart\t$hitlength\t$hitstrand\t$mismatch\n";
    } #while there are still lines in the file
    close (DATA);
} #foreach file

close (OUT); # stop output
