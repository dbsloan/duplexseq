#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use POSIX qw(floor);
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes a SAM file summarizing the result of reference mapping with
   SSCS or DCS output from duplex sequencing. It returns a summary of short indel
   calls.
   
   REQUIRED ARGUMENTS
   
   SAM File
         --sam
         File name of SAM output from bowtie2 mapping of DCS/SSCS sequences
         against a reference genome.
   
   Consensus Sequences Fasta File
         --seqs
         Fasta file containing DCS or SSCS reads from duplex sequencing.

   Output Name
         --output
         Base name for all output files (additional extensions will be added)


   OPTIONAL ARGUMENTS
  
   Excluded Reference Sequences
         --excluded_refs      
         Comma-delimited string of sequences within the reference fasta to exclude

   k-mer Database
         --kmer_db
         Database of k-mer counts as produced kmc_dump or similar program. Should
         be formated as a text file with one k-mer per line. Line format should
         be k-mer sequence followed by one white-space character followed by k-mer
         count. Counts of indentified variants will be reported. This feature may
         be useful for identifying contaminating sequences, such as apparent 
         mitochondrial variants that are actually derived from nuclear sequences 
         (i.e., numts). In this scenario the k-mer database might be generated
         from a standard shotgun library of total-cellular DNA.

   k-mer Size
         --kmer_size [default = 39]
         Size of k-mers. If a k-mer database is specified with --kmer_db to generate
         counts, this k-mer size must match the size used to generate the database
         with kmc or similar program.
";

our $SAM;
our $SEQS;
our $OUTPUT;
our $EXCLUDE_REFS;
our $KMER_DB;
our $KMER_SIZE=39;


GetOptions(
    'sam=s'  => \$SAM,
    'seqs=s'  => \$SEQS,
    'output=s'  => \$OUTPUT,
    'exclude_refs=s'  => \$EXCLUDE_REFS,
    'kmer_db=s'  => \$KMER_DB,
    'kmer_size=i'  => \$KMER_SIZE
);


print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";

$SAM or die ("\n$usage\n\nERROR: Must provide a SAM input filename with --sam.\n\n");
$SEQS or die ("\n$usage\n\nERROR: Must provide a DCS or SSCS fasta input filename with --seqs.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide a base name for output files with --output.\n\n");

#populate hash of kmer counts if a kmer database was provided
my %kmer_hash;
if ($KMER_DB){
	-e $KMER_DB or die ("\n$usage\n\nERROR: Specified k-mer database ($KMER_DB\) does not exist.\n\n");
	my $FH_KMER = open_file ($KMER_DB);
	my $kmer_first_line = <$FH_KMER>;
	chomp $kmer_first_line;
	my @sl = split (/\s/, $kmer_first_line);
	my $kmer_test_length = length($sl[0]);
	$kmer_test_length == $KMER_SIZE or die ("\n$usage\n\nERROR: The k-mer size in the specified database ($KMER_DB\) is $kmer_test_length\, which does not match specified k-mer size ($KMER_SIZE). Change the k-mer size to $kmer_test_length with --kmer_size or regenerate the k-mer database with the desired size.\n\n");
	$kmer_hash{$sl[0]} = $sl[1];
	while (<$FH_KMER>){
		chomp $_;
		my @sl2 = split (/\s/, $_);
		$kmer_hash{$sl2[0]} = $sl2[1];
	}
}

#this is different than for SNPs because cannot assume there is always a single bp in the middle
my $kmer_left_flank = floor ($KMER_SIZE / 2);
my $kmer_right_flank;
if ($kmer_left_flank == $KMER_SIZE / 2){
	$kmer_right_flank = $kmer_left_flank;
}else{
	$kmer_right_flank = $kmer_left_flank + 1;
}


##read in all full headers from the input DCS or SSCS fasta file
my ($headersRef, $seqsRef) = get_fasta_names_and_seqs($SEQS);
my @headers = @{$headersRef};

##build a hash of arrays to store info on each read family from fasta headers
my %readfamily_HoA;
foreach (@headers){
	chomp $_;
	my @sl = split (/\ /, $_);
	$readfamily_HoA{$sl[0]}[0] = $sl[1];
}

##store reference names to be excluded
my %excludes;
if ($EXCLUDE_REFS){
	my @exclude_array = split (/\,/, $EXCLUDE_REFS);
	foreach (@exclude_array){
		$excludes{$_} = 1;
	}
}

##Open output files and print header lines
my $FHO = open_output ("$OUTPUT\.indels.txt");
print $FHO "ID\tFamily Size\tReference\tRef Position\tRead Position\tMap Orientation\tIndel Type\tIndel Length\tIndel Seq\tk-mer";
if ($KMER_DB){
print $FHO "\tk-mer Count";
}
print $FHO "\tRead Sequence\n";


##loop through all lines in SAM file
my $FH = open_file ($SAM);
while (<$FH>){

	$_ =~ /^\@/ and next;
	my @sl = split (/\t/, $_);

	##skip lines that do not map
	$sl[2] eq '*' and next;

	##skip lines that map to reference sequences on the exclude list
	exists ($excludes{$sl[2]}) and next;

	##look at reads that show a single indel
	if ($sl[5] =~ /^(\d+)M(\d+)([ID])(\d+)M$/){
		
		my $left_match = $1;
		my $indel_length = $2;
		my $indel_type = $3;
		my $right_match = $4;
		
		## exclude reads with SNPs
		if ($_ =~ /XM\:i\:0\s/){
			
			my $read_pos = min ($left_match, $right_match) + 1;			
			my $ref_pos = $sl[3] + $left_match;
			my $map_orientation;
			if ($sl[1] == 0){
				$map_orientation = "F";
			}elsif ($sl[1] == 16){
				$map_orientation = "R";
			}else{
				die ("\nERROR: bit string not equal to 0 or 16:\n $_\n\n");
			}
			
			my $indel_seq;
			
			if ($indel_type eq 'D'){
				if ($_ =~ /MD\:Z\:\d+\^([ACGTN]+)\d/){
					$indel_seq = $1;
				}else{
					die ("\nERROR: Could not parse deletion sequence:\n $_\n\n");
				}
			}elsif($indel_type eq 'I'){
				$indel_seq = substr($sl[9], $left_match, $indel_length);
			}else{
				die ("\nERROR: Did not properly identify indel type:\n $_\n\n");		
			}
			
			my $kmer = "";
			
			if (length($sl[9]) < $kmer_left_flank + $kmer_right_flank + 1){
				$kmer = "NA";
			}
			
			#this extracts k-mer. Assuming there is enough sequence, it always pulls the pull left flank even if there is an insertion, so the kmer may not be centered around the variant.
			unless ($kmer eq "NA"){
				my $right_length_check = $right_match;
				if ($indel_type eq 'I'){
					$right_length_check += $indel_length;
				}
				if ($left_match >= $kmer_left_flank and $right_length_check >= $kmer_right_flank){
					$kmer = substr ($sl[9], $left_match - $kmer_left_flank, $KMER_SIZE)
				}elsif ($left_match < $kmer_left_flank){
					$kmer = substr ($sl[9], 0, $KMER_SIZE);
				}else{
					$kmer = substr ($sl[9], length($sl[9]) - $KMER_SIZE);
				}
			}
			
			print $FHO "$sl[0]\t", $readfamily_HoA{$sl[0]}[0], "\t$sl[2]\t$ref_pos\t$read_pos\t$map_orientation\t$indel_type\t$indel_length\t$indel_seq\t$kmer";
			
			if ($kmer eq "NA" or (not defined ($KMER_DB))){
				print $FHO "\t$sl[9]\n";
			}else{
				my $kmer_count = 0;
				if (exists ($kmer_hash{$kmer})){
					$kmer_count += $kmer_hash{$kmer};
				}
				if (exists ($kmer_hash{revcom($kmer)})){
					$kmer_count += $kmer_hash{revcom($kmer)};
				}
				print $FHO "\t$kmer_count\t$sl[9]\n";
			}		
		}	
	}	
}

print "\n" . (localtime) . "\nAnalysis with $0 complete.\n";