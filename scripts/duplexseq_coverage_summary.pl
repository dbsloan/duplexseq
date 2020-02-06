#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes cov_stats.txt files produced by the duplexseq pipeline
   and summarizes coverage by reference sequence and reference base (A, C,
   G, or T). It allows for trimming of positions near read ends in case a
   similar filter is being applied to variant calls.
   
   REQUIRED ARGUMENTS
   
   Input File List
         --input
         Name (and path if not in local directory) of file containing a list
         of cov_stats output files from duplexseq_sam_to_snps.pl. There should
         be one file name per line, and the path should be included if the
         files are not in the same directory as the SNP files. Alternatively,
         the --input_path option can be used to pre-append the same path to all
         file names.
   
   Output Name
         --output
         Name of output file
   
      
   OPTIONAL ARGUMENTS
 
   Path to Input Coverage Files
         --input_path   
         If the input cov files are not in the local directory and their full
         path is not specified in the file provided with the --input argument,
         then provide the path with this option. This should not be used if
         the list of file names already includes path information. It will not
         be applied to the location of the input list of filenames itself (that
         path should be provided as part of the --input option if necessary).

   Minimum Read Position
         --min_read_pos [default: 1]      
         This filter will exclude coverage found near the ends of reads. For
         example, setting this parameter to 11 will effectively remove any 
         coverage contribtued by the first or last 10 bps of a read. The
         deafult setting of 1 will not filter anything. This should be set to
         the same value used when filtering SNPs with the script
         duplexseq_snps_post_processing.pl.

                 
";

our $INPUT;
our $OUTPUT;
our $INPUT_PATH;
our $MIN_READ_POS = 1;

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";


GetOptions(
    'input=s'  => \$INPUT,
    'output=s'  => \$OUTPUT,
    'input_path=s'  => \$INPUT_PATH,
    'min_read_pos=i'  => \$MIN_READ_POS
);

##check for proper options specification
$INPUT or die ("\n$usage\n\nERROR: Must provide a file with input filenames with --input.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide a output filename with --output.\n\n");

##read in input file names from the provide list file
my @input_file_lines = file_to_array ($INPUT);
my @files;
if ($INPUT_PATH){
	substr ($INPUT_PATH,-1) eq '/' or $INPUT_PATH .= '/';
}

foreach (@input_file_lines){
	$_ =~ /^\s+$/ and next;
	chomp $_;
	if ($INPUT_PATH){
		push (@files, $INPUT_PATH . $_);
	}else{
		push (@files, $_);	
	}
}

my $FHO = open_output ($OUTPUT);

print $FHO "File\tReference\tTotal\tGC\tAT\tA\tC\tG\tT\n";

##loop through every provided coverage file
foreach (@files){

	my @coverage_lines = file_to_array ($_);

	my %cov_HoH; #hash of hashes with sequence name as top key and base (A, C, G, or T) as secondary key. Value is total coverage count;
	my $in_header_lines = 1;	
	foreach my $line (@coverage_lines){
		
		if ($in_header_lines){
			if ($line =~ /^Sequence\tRead\ Position/){
				$in_header_lines = 0;
			}
		}else{
			chomp $line;
			my @sl = split (/\t/, $line);
			$sl[1] >= $MIN_READ_POS or next;
			$sl[3] =~ /^[ACGT]$/ or die ("\nERROR: non-ACGT character in the following line in $_:\n$line\n\n");
			$cov_HoH{$sl[0]}->{$sl[3]} += $sl[4];
		}		
	}

	foreach my $seq_name (sort keys %cov_HoH){
	
		my $a = $cov_HoH{$seq_name}->{"A"}; 
		my $c = $cov_HoH{$seq_name}->{"C"}; 
		my $g = $cov_HoH{$seq_name}->{"G"}; 
		my $t = $cov_HoH{$seq_name}->{"T"};
		my $gc = $c + $g;
		my $at = $a + $t;
		my $total = $gc + $at;
		
		print $FHO "$_\t$seq_name\t$total\t$gc\t$at\t$a\t$c\t$g\t$t\n";			
	}	


}

print "\n" . (localtime) . "\nAnalysis Complete\n\n";
