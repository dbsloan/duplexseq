#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use sloan;

my $usage =
"\nUsage: perl $0 [options/arguments]
   
   This script takes output fasta files from duplexseq_process_reads.pl
   and applies filters and end-trimming
   
   REQUIRED ARGUMENTS
   
   Fasta File Input
         --input
         A fasta output file from duplexseq_process_reads.pl
   
   Output Name
         --output
         Name of new fasta file where output will be written
   
   
   OPTIONAL ARGUMENTS
   
   Unmerged Input
         --unmerged [default: off] 
         Add this flag if the input fasta is an unmmerged file.

   Trim Non-Overlap Regions
         --overlap_only [default: off]      
         Delete sequnece regions that were not in the overlap region of the
         merged read. This option is ignored if --unmerged is specified.

   Discard Reads with Ns
         --discard_n_reads [default: off] 
         Add this flag to exlude any read that contains an N.

   Trim Ligation Ends
         --endtrim [default: 0]      
         Remove this many bases from start of reads. For merged reads,
         this will be removed from both ends. For unmerged reads, only
         the 5' end will be trimmed.
        
   Minimum Length
         --min_len [default: 1]   
         Exclude reads shorter than this threshold after trimming steps are
         applied.
      
";

our $INPUT;
our $OUTPUT;
our $UNMERGED;
our $OVERLAP_ONLY;
our $DISCARD_N_READS;
our $ENDTRIM = 0;
our $MIN_LEN = 1;

GetOptions(
    'input=s'  => \$INPUT,
    'output=s'  => \$OUTPUT,
    'unmerged'  => \$UNMERGED,    
    'overlap_only'  => \$OVERLAP_ONLY,    
    'discard_n_reads'  => \$DISCARD_N_READS,
    'endtrim=i'  => \$ENDTRIM,
    'min_len=i'  => \$MIN_LEN
);

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";

##Check parameter inputs
$INPUT or die ("\n$usage\n\nERROR: Must provide a fasta input filename with --input.\n\n");
-e $INPUT or die ("\n$usage\n\nERROR: File specified with --input does not exist: $INPUT.\n\n");
(int($ENDTRIM) == $ENDTRIM and $ENDTRIM >= 0) or die ("\n$usage\n\nERROR: endtrim must be a non-negative integer\.\n\n");
(int($MIN_LEN) == $MIN_LEN and $MIN_LEN > 0) or die ("\n$usage\n\nERROR: min_len must be a positive integer\.\n\n");
$UNMERGED and $OVERLAP_ONLY and (print STDERR ("\nWARNING: Ignoring --overlap_only because --unmerged was specified.\n\n") and $OVERLAP_ONLY = 0);

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nReading in input file $INPUT and applying specified filtering and trimming steps.\n\n";

my %fasta = fasta2hash($INPUT);

my $family_count = 0;
my $out_family_count = 0;
my $n_exclude_count = 0;
my $min_len_count = 0;
my $FHO = open_output ($OUTPUT);

#loop through input sequences individually
foreach (sort keys %fasta){
	++$family_count;
	my $seq = $fasta{$_};
	chomp $_;
	my @sl = split (/\ /, $_);
	
	#discard reads containing an N if the --discard_n_reads option is specified
	$DISCARD_N_READS and $seq =~ /N/ and ++$n_exclude_count and next;
	
	#trim read ends if --endtrim is specified
	if ($ENDTRIM){
		if ($UNMERGED){
			length ($seq) > $ENDTRIM or (++$min_len_count and next);
			$seq = substr ($seq, $ENDTRIM);
		}else{
			length ($seq) > 2*$ENDTRIM or (++$min_len_count and next);
			$seq = substr ($seq, $ENDTRIM, -1*$ENDTRIM);
			$sl[2] = max(1, $sl[2] - $ENDTRIM);
			$sl[3] = min (length($seq), $sl[3] - $ENDTRIM);	
		}
	}
	
	#trim to overlapping region within merged reads if --overlap_only is specified
	if ($OVERLAP_ONLY){
		$seq = substr ($seq, $sl[2] - 1, $sl[3] - $sl[2] + 1);
		$sl[2] = 1;
		$sl[3] = length ($seq);
	}
	
	#discard trimmed reads with length less than min_len
	length ($seq) > $MIN_LEN or (++$min_len_count and next);

	++$out_family_count;
	#print trimmed seqs to output fasta file
	print $FHO ">", join (' ', @sl), "\n$seq\n";
}

print "\n" . (localtime) . "\nProcessing Complete\n\n";
print "Read counts and filtering:\n";
print "$family_count\tNumber of input seq families\n";
print "$out_family_count\tNumber of families written to output after processing\n";
$DISCARD_N_READS and print "$n_exclude_count\tExcluded families that contained N in sequence\n";
print "$min_len_count\tExcluded families with a sequence length < $MIN_LEN after trimming\n";