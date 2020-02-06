#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IPC::Cmd qw(can_run);
use List::Util qw(min max);
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes raw Illumina fastq input files from a duplex sequencing 
   library and returns single-standed and double-stranded consensus sequences.
   
   REQUIRED ARGUMENTS
   
   R1 Fastq File
         --r1_fastq
         File containing Illumina read 1 sequences. Can be gzipped.
   
   R2 Fastq File
         --r2_fastq
         File containing Illumina read 2  sequences. Can be gzipped.

   Output Name
         --output
         Base name for all output files (additional extensions will be added)
   
   
   OPTIONAL ARGUMENTS
 
   Minimum SSCS Family Size
         --min_sscs [default: 3]      
         Minimum number of reads in a family to generate a single-stranded
         consensus sequence

   Minimum agreement to call SSCS base
         --min_agree [default: 0.8]      
         Proportion of reads within a read family that must agree to make an 
         SSCS base call.

   Minimum SSCS applied to only one strand
         --min_sscs_one_strand [default: off]      
         Add this flag if min SSCS read number is only required for one of the
         two complementary read families (i.e., a double-stranded consensus can 
         be built with only a single read in the complementary family).

   Cutadapt Executable
         --cutadapt_exe [default: cutadapt]      
         Full path to cutadapt if not in default PATH

   Cutadapt R1 Adapter
         --cutadapt_adap1 [default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]   
         Adapter sequences to trim from read 1

   Cutadapt R2 Adapter
         --cutadapt_adap2 [default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]   
         Adapter sequences to trim from read 2

   Cutadapt Error Tolerance
         --cutadapt_err_tol [default: 0.15]   
         Error tolerance for cutadapt trimming

   Cutadapt Cores
         --cutadapt_cores [default: 4]   
         Number of cores to request for cutadapt step. This may need to be set
         to 1 if using a cutadapt version with Python v2.

   Cutadapt Minimum Length
         --cutadapt_min_len [default: 75]   
         Value for minimum_length parameter in cutadapt

   Cutadapt Quality Score for Trimming
         --cutadapt_trimq [default: 20]   
         Quality score used for q parameter in cutadapt
         
   Turn Off Cutadapt Trimming
         --disable_cutadapt 
         Add this flag to turn off primer trimming with cutadapt.

   BBmerge Executable
         --bbmerge_exe [default: bbmerge.sh]      
         Full path to bbmerge.sh if not in default PATH
         
   BBmerge Minimum Overlap
         --bbmerge_minoverlap [default: 30]   
         Sequence length for minoverlap parameter in bbmerge.sh

   BBmerge Mismatches
         --bbmerge_mismatches [default: 5]   
         Maximum number of mismathces for mismatches parameter in bbmerge.sh
   
   Random Barcode Length
         --barcode_len [default: 12] 
         Number of Ns in the random barcode (unique molecular identifier) at
         beginning of each read.

   Linker Sequence
         --linker [default: TGACT] 
         Number of Ns in the random barcode (unique molecular identifier) at
         beginning of each read.
 
   Disable Repetitive Barcode Filter
    	--disable_rep_filter
    	Add this flag to turn off the default filter that excludes reads with
    	barcodes that are just one long homopolymer (e.g. AAAAAAAAAAAA).

   Disable Filtering of Barcodes with Ns
    	--disable_n_filter
    	Add this flag to turn off the default filter that excludes reads with
    	barcodes that contain an N.

   Exclude Reads with Low Barcode Quality   
    	--min_barcode_qual [default: 20]
    	Discard reads that have basecalls in barcode with quality value lower
    	than this.
    	
   Illumina Phred Quality Score Version  
    	--phred_offset [default: 33]
    	ASCII offset value for quality score encoding. 33 is standard for
    	current Illumina runs.

   Delete Fastq Files
    	--delete_intermediate_fastqs
    	Add this flag to delete the large intermediate fastq files and only
    	save the final fasta consensus sequence files.

   Combine DCS Fasta Output Files
    	--combine_dcs
    	Add this flag to print all DCS output to a single fasta file
    	regardless of whether they were merged by bbmerge.

   Suppress SSCS Fasta Output Files
    	--suppress_sscs
    	Add this flag to avoid producing SSCS fasta output.

   Filtering Optical Duplicates    
    	--min_optical_dist [default: 0]
    	Minimum pixel distance for retaining a read with the same barcode as 
    	another read on the same tile within the same Illumina flow cell.
    	Specifying non-zero values here will help eliminate optical/clustering
    	duplicates.
                 
";

our $R1_FASTQ;
our $R2_FASTQ;
our $OUTPUT;
our $MIN_SSCS = 3;
our $MIN_AGREE = 0.8;
our $MIN_SSCS_ONE_STRAND = 0;
our $BBMERGE_EXE = "bbmerge.sh";
our $BBMERGE_MINOVERLAP = 30;
our $BBMERGE_MISMATCHES = 5;
our $DISABLE_CUTADAPT;
our $CUTADAPT_EXE = "cutadapt";
our $CUTADAPT_ERR_TOL = 0.15;
our $CUTADAPT_CORES = 4;
our $CUTADAPT_MIN_LEN = 75;
our $CUTADAPT_ADAP1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
our $CUTADAPT_ADAP2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
our $CUTADAPT_TRIMQ = 20;
our $BARCODE_LEN = 12;
our $LINKER = "TGACT";
our $DISABLE_REP_FILTER;
our $MIN_OPTICAL_DIST = 0;
our $MIN_TILE_EDGE_DIST = 0;
our $MIN_BAROCODE_QUAL = 20;
our $PHRED_OFFSET = 33;
our $DISABLE_N_FILTER;
our $DELETE_INTERMEDIATE_FASTQS;
our $COMBINE_DCS;
our $SUPPRESS_SSCS;

##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";


GetOptions(
    'r1_fastq=s'  => \$R1_FASTQ,
    'r2_fastq=s'  => \$R2_FASTQ,
    'output=s'  => \$OUTPUT,
    'min_sscs=i'  => \$MIN_SSCS,    
    'min_agree=i'  => \$MIN_AGREE,    
    'min_sscs_one_strand'  => \$MIN_SSCS_ONE_STRAND,    
    'bbmerge_exe=s'  => \$BBMERGE_EXE,
    'bbmerge_minoverlap=i'  => \$BBMERGE_MINOVERLAP,
    'bbmerge_mismatches=i'  => \$BBMERGE_MISMATCHES,
    'disable_cutadapt'  => \$DISABLE_CUTADAPT,
    'cutadapt_exe=s'  => \$CUTADAPT_EXE,
    'cutadapt_err_tol=i'  => \$CUTADAPT_ERR_TOL,
    'cutadapt_cores=i'  => \$CUTADAPT_CORES,
    'cutadapt_min_len=i'  => \$CUTADAPT_MIN_LEN,
    'cutadapt_adap1=s'  => \$CUTADAPT_ADAP1,
    'cutadapt_adap2=s'  => \$CUTADAPT_ADAP2,
    'cutadapt_trimq=i'  => \$CUTADAPT_TRIMQ,
    'barcode_len=i'  => \$BARCODE_LEN,
    'linker=s'  => \$LINKER,
    'disable_rep_filter'  => \$DISABLE_REP_FILTER,
	'min_optical_dist=i'  => \$MIN_OPTICAL_DIST,
    'min_tile_edge_dist=i'  => \$MIN_TILE_EDGE_DIST,
    'min_barcode_qual=i'  => \$MIN_BAROCODE_QUAL,
    'phred_offset=i'  => \$PHRED_OFFSET,
    'disable_n_filter'  => \$DISABLE_N_FILTER,
    'delete_intermediate_fastqs'  => \$DELETE_INTERMEDIATE_FASTQS,
    'combine_dcs'  => \$COMBINE_DCS,
    'suppress_sscs'  => \$SUPPRESS_SSCS
);

##Check dependencies
unless ($DISABLE_CUTADAPT){
	can_run ($CUTADAPT_EXE) or die ("\n$usage\n\nERROR: cutadapt is enabled, but could not find following executable in PATH: $CUTADAPT_EXE. Provide a path and executable name with --cutadapt_exe or disable with --disable_cutadapt.\n\n");
}
can_run ($BBMERGE_EXE) or die ("\n$usage\n\nERROR: Could not find following executable in PATH: $BBMERGE_EXE. Provide a path and executable name with --bbmerge_exe.\n\n");

##Check parameter inputs
(int($BBMERGE_MINOVERLAP) == $BBMERGE_MINOVERLAP and $BBMERGE_MINOVERLAP >= 0) or die ("\n$usage\n\nERROR: bbmerge_minoverlap must be a non-negative integer\.\n\n");
(int($CUTADAPT_TRIMQ) == $CUTADAPT_TRIMQ and $CUTADAPT_TRIMQ >= 0) or die ("\n$usage\n\nERROR: cutadapt_trimq must be a non-negative integer\.\n\n");
(int($BBMERGE_MISMATCHES) == $BBMERGE_MISMATCHES and $BBMERGE_MISMATCHES >= 0) or die ("\n$usage\n\nERROR: bbmerge_mismatches must be a non-negative integer\.\n\n");
(int($BARCODE_LEN) == $BARCODE_LEN and $BARCODE_LEN > 0) or die ("\n$usage\n\nERROR: barcode_len must be a positive integer\.\n\n");
(int($MIN_BAROCODE_QUAL) == $MIN_BAROCODE_QUAL and $MIN_BAROCODE_QUAL >= 0) or die ("\n$usage\n\nERROR: min_barcode_qual must be a non-negative integer\.\n\n");
(int($PHRED_OFFSET) == $PHRED_OFFSET and $PHRED_OFFSET >= 0) or die ("\n$usage\n\nERROR: phred_offset must be a non-negative integer\.\n\n");
($CUTADAPT_ERR_TOL >= 0 and $CUTADAPT_ERR_TOL <= 1) or die ("\n$usage\n\nERROR: cutadapt_err_tol must be a number from 0 to 1\.\n\n");
(int($CUTADAPT_CORES) == $CUTADAPT_CORES and $CUTADAPT_CORES >= 1) or die ("\n$usage\n\nERROR: cutadapt_cores must be a positive integer\.\n\n");
(int($CUTADAPT_MIN_LEN) == $CUTADAPT_MIN_LEN and $CUTADAPT_MIN_LEN >= 0) or die ("\n$usage\n\nERROR: cutadapt_min_len must be a non-negative integer\.\n\n");
(int($MIN_OPTICAL_DIST) == $MIN_OPTICAL_DIST and $MIN_OPTICAL_DIST >= 0) or die ("\n$usage\n\nERROR: min_optical_dist must be a non-negative integer\.\n\n");


##Check file inputs
$R1_FASTQ or die ("\n$usage\n\nERROR: Must provide a fastq input filename with --r1_fastq.\n\n");
-e $R1_FASTQ or die ("\n$usage\n\nERROR: File specified with --r1_fastq does not exist: $R1_FASTQ.\n\n");
$R2_FASTQ or die ("\n$usage\n\nERROR: Must provide a fastq input filename with --r2_fastq.\n\n");
-e $R2_FASTQ or die ("\n$usage\n\nERROR: File specified with --r2_fastq does not exist: $R2_FASTQ.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide output name with --output.\n\n");

##Run cutadapt
unless ($DISABLE_CUTADAPT){
	my $cutadapt_string = "$CUTADAPT_EXE -a $CUTADAPT_ADAP1 -A $CUTADAPT_ADAP2 -q $CUTADAPT_TRIMQ --minimum-length $CUTADAPT_MIN_LEN --cores $CUTADAPT_CORES -e $CUTADAPT_ERR_TOL -o $OUTPUT\.trim1.fq -p $OUTPUT\.trim2.fq $R1_FASTQ $R2_FASTQ > $OUTPUT\.cutadapt_log.txt";
	print "\n" . (localtime) . "\nRunning $CUTADAPT_EXE to trim for Illumina adapters and low-quality ends, using the following command:\n\n$cutadapt_string\n\n";
	system ($cutadapt_string);
}

##Run bbmerge
my $bb_input1 = $R1_FASTQ;
my $bb_input2 = $R2_FASTQ;
unless ($DISABLE_CUTADAPT){
	$bb_input1 = "$OUTPUT\.trim1.fq";
	$bb_input2 = "$OUTPUT\.trim2.fq";
}
my $bbmerge_string = "$BBMERGE_EXE in1=$bb_input1 in2=$bb_input2 out=$OUTPUT\.merged.fq outu=$OUTPUT\.unmerged1.fq outu2=$OUTPUT\.unmerged2.fq ihist=$OUTPUT\.bbmerge_hist.txt ordered=t qtrim=r minoverlap=$BBMERGE_MINOVERLAP mismatches=$BBMERGE_MISMATCHES 2> $OUTPUT\.bbmerge_log.txt";
print "\n" . (localtime) . "\nRunning $BBMERGE_EXE to collapse trimmed reads with the following command:\n\n$bbmerge_string\n\n";
system ($bbmerge_string);

#make HoA to store R1 and R2 lengths before bbmerge occurred. Note storing MAX lengths for each family and not distinguishing between alpha-beta and beta-alpha
my %premerge_length_HoA;
my $FHT1 = open_file("$OUTPUT\.trim1.fq");
my $FHT2 = open_file("$OUTPUT\.trim2.fq");
print "\n" . (localtime) . "\nStoring pre-merged lengths of trimmed reads to determine length of merge overlap...\n\n";

while (my $line = <$FHT1>){
	$line = <$FHT1>;
	my $line2 = <$FHT2>;
	$line2 = <$FHT2>;

	chomp $line;
	chomp $line2;
	
	my $bc1 = substr($line, 0, $BARCODE_LEN);	
	my $bc2 = substr($line2, 0, $BARCODE_LEN);
	
	my $cmp_test = $bc1 cmp $bc2;
	
	my $duplex_umi;	
	if ($cmp_test == -1){
		$duplex_umi = $bc1 . $bc2;
	}elsif ($cmp_test == 1){
		$duplex_umi = $bc2 . $bc1;
	}else{
		next;	
	}
	
	#for storing read lengths, flip positions for beta-alpha
	if (exists $premerge_length_HoA{$duplex_umi}){
		if ($cmp_test == -1){
			if (length ($line) > $premerge_length_HoA{$duplex_umi}[0]){
				$premerge_length_HoA{$duplex_umi}[0] = length ($line);
			}
			if (length ($line2) > $premerge_length_HoA{$duplex_umi}[1]){
				$premerge_length_HoA{$duplex_umi}[1] = length ($line2);
			}
		}else{
			if (length ($line) > $premerge_length_HoA{$duplex_umi}[1]){
				$premerge_length_HoA{$duplex_umi}[1] = length ($line);
			}
			if (length ($line2) > $premerge_length_HoA{$duplex_umi}[0]){
				$premerge_length_HoA{$duplex_umi}[0] = length ($line2);
			}		
		}
	}else{	
		if ($cmp_test == -1){
			$premerge_length_HoA{$duplex_umi}[0] = length ($line);
			$premerge_length_HoA{$duplex_umi}[1] = length ($line2);
		}else{
			$premerge_length_HoA{$duplex_umi}[1] = length ($line);
			$premerge_length_HoA{$duplex_umi}[0] = length ($line2);		
		}
	}
	$line = <$FHT1>;
	$line = <$FHT1>;
	$line = <$FHT2>;
	$line = <$FHT2>;
}
close $FHT1;
close $FHT2;

if ($DELETE_INTERMEDIATE_FASTQS){
	unless ($DISABLE_CUTADAPT){
		unlink "$OUTPUT\.trim1.fq";
		unlink "$OUTPUT\.trim2.fq";
	}
} 

##Create pseudo read2 seqs for the merged reads and filter based on presence of linker sequences. Trim 3' barcodes/linkers off of merged reads
print "\n" . (localtime) . "\nProcessing bbmerge output to screen for expected linker sequences and to trim 3' barcodes...\n\n";
my $FHO1 = open_output("$OUTPUT\_linker_filt1.fq");
my $FHO2 = open_output("$OUTPUT\_linker_filt2.fq");
my $FHM = open_file ("$OUTPUT\.merged.fq");
my $FHUM1 = open_file ("$OUTPUT\.unmerged1.fq");
my $FHUM2 = open_file ("$OUTPUT\.unmerged2.fq");

my $merged_read_count = 0;
my $merged_linker_filter = 0;

while (my $header_line = <$FHM>){#loop over read blocks in merged fastq file from bbmerge
	++$merged_read_count;
	my $seq_line = <$FHM>;
	my $plus_line = <$FHM>;
	my $qual_line = <$FHM>;
	chomp $seq_line;
	
	#confirm expected linker seqs on both ends
	if (substr($seq_line, $BARCODE_LEN, length($LINKER)) eq $LINKER and substr($seq_line, -1 * ($BARCODE_LEN + length($LINKER)), length($LINKER)) eq revcom($LINKER)){
		chomp $header_line;
	
		print $FHO1 "$header_line\n";
		my @sl = split (/\s+/, $header_line);
		#flag pseudo read 2 headers with "merge" to distinguish from unmerged reads
		print $FHO2 "$sl[0] 2", substr($sl[1], 1), "\tmerge\n";
			
		print $FHO1 substr($seq_line, 0, length($seq_line) - ($BARCODE_LEN + length($LINKER))), "\n";
		print $FHO2 revcom(substr($seq_line, $BARCODE_LEN + length($LINKER))), "\n";
		
		print $FHO1 $plus_line;
		print $FHO2 $plus_line;
	
		chomp $qual_line;
		
		print $FHO1 substr($qual_line, 0, length($qual_line) - ($BARCODE_LEN + length($LINKER))), "\n";
		my $rev_qual_line = reverse (substr($qual_line, $BARCODE_LEN + length($LINKER)));
		print $FHO2 "$rev_qual_line\n";
	}else{
		++$merged_linker_filter;
	}
}

##Filter unmerged reads based on presence of linker sequences
my $unmerged_read_count = 0;
my $unmerged_linker_filter = 0;

while (my $header_line1 = <$FHUM1>){#loop over read blocks in unmerged fastq files from bbmerge
	++$unmerged_read_count;
	my $seq_line1 = <$FHUM1>;
	my $plus_line1 = <$FHUM1>;
	my $qual_line1 = <$FHUM1>;
	my $header_line2 = <$FHUM2>;
	my $seq_line2 = <$FHUM2>;
	my $plus_line2 = <$FHUM2>;
	my $qual_line2 = <$FHUM2>;

	#confirm expected linker seqs near 5' ends on both reads
	if (substr($seq_line1, $BARCODE_LEN, length($LINKER)) eq $LINKER and substr($seq_line2, $BARCODE_LEN, length($LINKER)) eq $LINKER ){		
		print $FHO1 "$header_line1$seq_line1$plus_line1$qual_line1";
		print $FHO2 "$header_line2$seq_line2$plus_line2$qual_line2";
	}else{
		++$unmerged_linker_filter;
	}
}

close $FHO1;
close $FHO2;
close $FHM;
close $FHUM1;
close $FHUM2;

if ($DELETE_INTERMEDIATE_FASTQS){
	unlink "$OUTPUT\.merged.fq";
	unlink "$OUTPUT\.unmerged1.fq";
	unlink "$OUTPUT\.unmerged2.fq";
} 

print "Total merged reads:\t $merged_read_count\n";
print "Total unmerged read pairs:\t $unmerged_read_count\n";
print "Merged reads removed because of missing linker sequence ($LINKER):\t $merged_linker_filter\n";
print "Unmerged read pairs removed because of missing linker sequence ($LINKER):\t $unmerged_linker_filter\n";

print "\n" . (localtime) . "\nProcessing filtered reads to aggregate read families...\n\n";

my $FHR1 = open_file ("$OUTPUT\_linker_filt1.fq");
my $FHR2 = open_file ("$OUTPUT\_linker_filt2.fq");

my $identical_bc_count = 0;
my $n_bc_count = 0;
my $homopolymer_bc_count = 0;
my $lowqual_bc_count = 0;
my $postlinker_readcount = 0;

#complex data structure. Top-level key for hash is the 24 bp ID generated by concatenating the two 12 bp UMIs
#Next level is a 2 element array: 1st element is alpha-beta group [array], second element is beta-alpha group [array]
#Next level is a 2 element array: 1st element is R1 seqs, second element is R2 seqs.
##But for beta-alpha fams read 1 and read 2 are swapped so they are directly comparable.
#Bottom level is array of seqs for that read family and read "type"
my %seq_HoAoAoA;

#hash keep track of whether reads in a UMI family were unmerged (0) or merged (1). Exclude unmerged reads when other family members are merged
my %merge_history;
my $hetero_merge_excludes = 0;

#keep track of read headers to parse optical distance data if a filter has been applied
my %optical_headers_HoA;
my $opt_count = 0;

while (my $line1_1 = <$FHR1>){
	my $line1_2 = <$FHR1>;
	my $line1_3 = <$FHR1>;
	my $line1_4 = <$FHR1>;
	my $line2_1 = <$FHR2>;
	my $line2_2 = <$FHR2>;
	my $line2_3 = <$FHR2>;
	my $line2_4 = <$FHR2>;

	++$postlinker_readcount;

	chomp $line1_2;
	chomp $line2_2;
	
	my $bc1 = substr($line1_2, 0, $BARCODE_LEN);	
	my $bc2 = substr($line2_2, 0, $BARCODE_LEN);
	
	#Remove barcodes that are just one long homopolymer
	unless ($DISABLE_REP_FILTER){
		$bc1 =~ /^(.)\1+$/ and ++$homopolymer_bc_count and next;
		$bc2 =~ /^(.)\1+$/ and ++$homopolymer_bc_count and next;
	}

	#Remove barcodes that contain an N
	unless ($DISABLE_N_FILTER){
		$bc1 =~ /N/ and ++$n_bc_count and next;
		$bc2 =~ /N/ and ++$n_bc_count and next;
	}

	#Exclude if barcode contains any sites with quality score below threshold
	my $low_qual = 0;
	if ($MIN_BAROCODE_QUAL){
		for (my $i=0; $i < $BARCODE_LEN; ++$i){
			ord(substr($line1_4, $i, 1)) - $PHRED_OFFSET < $MIN_BAROCODE_QUAL and $low_qual=1 and ++$lowqual_bc_count and last;
			ord(substr($line2_4, $i, 1)) - $PHRED_OFFSET < $MIN_BAROCODE_QUAL and $low_qual=1 and ++$lowqual_bc_count and last;
		}
		$low_qual and next;
	}	

	#extract sequence after barcode and linker
	my $subseq1 = substr ($line1_2, $BARCODE_LEN + length($LINKER));
	my $subseq2 = substr ($line2_2, $BARCODE_LEN + length($LINKER));

	#concatenate barcodes based on alphabetical order
	my $cmp_test = $bc1 cmp $bc2;
	my $duplex_umi;	
	if ($cmp_test == -1){
		$duplex_umi = $bc1 . $bc2;
	}elsif ($cmp_test == 1){
		$duplex_umi = $bc2 . $bc1;
	}else{
		#if alpha and beta barcodes are identical to each other, skip them
		++$identical_bc_count;
		next;	
	}

	#Apply optical distance filter by skipping reads that are two physically close to a previously processed read with same UMI.
	if ($MIN_OPTICAL_DIST){
		my $opt_flag = 0;
		if (exists ($optical_headers_HoA{$duplex_umi})){
			my @sh1 = split (/\s/, $line1_1);
			my @split_read1 = split (/\:/, $sh1[0]);
			foreach my $header (@{$optical_headers_HoA{$duplex_umi}}){
				my @sh2 = split (/\s/, $header);
				my @split_read2 = split (/\:/, $sh2[0]);
				my $tile1 = $split_read1[2] . "_" . $split_read1[3] . "_" . $split_read1[4];
				my $tile2 = $split_read2[2] . "_" . $split_read2[3] . "_" . $split_read2[4];
			
				if ($tile1 eq $tile2){
					if (sqrt(abs($split_read1[5] - $split_read2[5])**2 + abs($split_read1[6] - $split_read2[6])**2) < $MIN_OPTICAL_DIST){
						$opt_flag = 1;
						last;
					}
				}			
			}
		}
		$opt_flag and ++$opt_count and next;
		push (@{$optical_headers_HoA{$duplex_umi}}, $line1_1);
	}


	#check if the reads were merged by bbmerge (flag was added in R2 header) and store in merge history hash
	#skip unmerged reads if merged reads are already present for that family
	#deleted previous unmerged reads if this read is merged
	my $merge = 0;
	$line2_1 =~ /\tmerge\n/ and $merge = 1;
	if (exists ($merge_history{$duplex_umi})){
		if ($merge == 0 and $merge_history{$duplex_umi} == 1){
			++$hetero_merge_excludes;
			next;
		}elsif ($merge == 1 and $merge_history{$duplex_umi} == 0){
			if (exists ($seq_HoAoAoA{$duplex_umi}[0][0])){
				$hetero_merge_excludes += scalar (@{$seq_HoAoAoA{$duplex_umi}[0][0]});
			}
			if (exists ($seq_HoAoAoA{$duplex_umi}[1][0])){
				$hetero_merge_excludes += scalar (@{$seq_HoAoAoA{$duplex_umi}[1][0]});
			}
			delete $seq_HoAoAoA{$duplex_umi};
		}
	}
	$merge_history{$duplex_umi} = $merge;
	
	#store reads in complex data structure described above
	if ($cmp_test == -1){
		push (@{$seq_HoAoAoA{$duplex_umi}[0][0]}, $subseq1);
		push (@{$seq_HoAoAoA{$duplex_umi}[0][1]}, $subseq2);
	}elsif ($cmp_test == 1){
		push (@{$seq_HoAoAoA{$duplex_umi}[1][0]}, $subseq2);
		push (@{$seq_HoAoAoA{$duplex_umi}[1][1]}, $subseq1);	
	}
} 

##Loop through read families and generate consensus sequences and summary of read-family sizes
print "\n" . (localtime) . "\nGenerating consensus sequences and read-family summary stats...\n\n";

my $family_count = 0;
my $double_strands = 0;
my $single_strands = 0;
my %double_HoH;
my %single_hash;
my $dcs_length_mismatches = 0;
my $lost_start_lengths = 0;


$SUPPRESS_SSCS or my $FHSM0 = open_output ("$OUTPUT\.SSCS.merge.fas");
$SUPPRESS_SSCS or my $FHSM1 = open_output ("$OUTPUT\.SSCS.high.merge.fas");
$SUPPRESS_SSCS or my $FHSM2 = open_output ("$OUTPUT\.SSCS.low.merge.fas");

$SUPPRESS_SSCS or my $FHSU0_1 = open_output ("$OUTPUT\.SSCS.unmerge.R1.fas");
$SUPPRESS_SSCS or my $FHSU1_1 = open_output ("$OUTPUT\.SSCS.high.unmerge.R1.fas");
$SUPPRESS_SSCS or my $FHSU2_1 = open_output ("$OUTPUT\.SSCS.low.unmerge.R1.fas");

$SUPPRESS_SSCS or my $FHSU0_2 = open_output ("$OUTPUT\.SSCS.unmerge.R2.fas");
$SUPPRESS_SSCS or my $FHSU1_2 = open_output ("$OUTPUT\.SSCS.high.unmerge.R2.fas");
$SUPPRESS_SSCS or my $FHSU2_2 = open_output ("$OUTPUT\.SSCS.low.unmerge.R2.fas");


$COMBINE_DCS and my $FHDC = open_output ("$OUTPUT\.DCS.fas");
$COMBINE_DCS or my $FHDM = open_output ("$OUTPUT\.DCS.merge.fas");
$COMBINE_DCS or my $FHDU_1  = open_output ("$OUTPUT\.DCS.unmerge.R1.fas");
$COMBINE_DCS or my $FHDU_2  = open_output ("$OUTPUT\.DCS.unmerge.R2.fas");


foreach my $id (sort keys %seq_HoAoAoA){
	++$family_count;
	
	my $high_strand;
	my $low_strand;
	my $dcs = 0;
		
	#determine whether both strands are represented and store which strand has the higher number of reads (use alpha-beta if it's a tie)
	if (exists ($seq_HoAoAoA{$id}[0]) and exists ($seq_HoAoAoA{$id}[1])){
		++$double_strands;
		++$double_HoH{scalar (@{$seq_HoAoAoA{$id}[0][0]})}->{scalar (@{$seq_HoAoAoA{$id}[1][0]})};
		$dcs = 1;
		
		if (scalar (@{$seq_HoAoAoA{$id}[0][0]}) >=  scalar (@{$seq_HoAoAoA{$id}[1][0]})){
			$high_strand = 0;			
			$low_strand = 1;			
		}else{
			$high_strand = 1;			
			$low_strand = 0;			
		}
		
	}else{
		++$single_strands;
		if (exists $seq_HoAoAoA{$id}[0]){
			++$single_hash{scalar (@{$seq_HoAoAoA{$id}[0][0]})};
			$high_strand = 0;
		}else{
			++$single_hash{scalar (@{$seq_HoAoAoA{$id}[1][0]})};
			$high_strand = 1;
		}
	}
	
	##keep track of R1 vs R2 history, which will need to determine strand asymmetry for meaqures such as oxoQ or GIV.
	my $print_strand1;
	my $print_strand2;
	if ($high_strand == 0){
		$print_strand1 = "F";
		$print_strand2 = "R";
	}else{
		$print_strand1 = "R";
		$print_strand2 = "F";
	}

	#Store seq sets as arrays
	my @seqs_high = @{$seq_HoAoAoA{$id}[$high_strand][0]};
	my @seqs_low;
	my @seqs_high_R2;
	my @seqs_low_R2;

	scalar (@seqs_high) >= $MIN_SSCS or next;
	$dcs and @seqs_low = @{$seq_HoAoAoA{$id}[$low_strand][0]};
	
	#create R2 arrays if these reads were not merged
	unless ($merge_history{$id}){
		@seqs_high_R2 = @{$seq_HoAoAoA{$id}[$high_strand][1]};
		$dcs and @seqs_low_R2 = @{$seq_HoAoAoA{$id}[$low_strand][1]};
	}
	
	#check whether the the other strand has enough reads unless the min_sscs_one_strand flag is set
	if ($dcs){
		unless ($MIN_SSCS_ONE_STRAND){
			scalar (@seqs_low) >= $MIN_SSCS or $dcs=0;
		}
	}
	
	my $sscs_high;
	my $sscs_low;
	my $sscs_high_R2;
	my $sscs_low_R2;
	
	$sscs_high = generate_consensus (\@seqs_high, $MIN_SSCS, $MIN_AGREE);
	$dcs and $sscs_low = generate_consensus (\@seqs_low, $MIN_SSCS, $MIN_AGREE);
	unless ($merge_history{$id}){
		$sscs_high_R2 = generate_consensus (\@seqs_high_R2, $MIN_SSCS, $MIN_AGREE);
		$dcs and $sscs_low_R2 = generate_consensus (\@seqs_low_R2, $MIN_SSCS, $MIN_AGREE);
	}
	
	my $dcs_seq;
	my $dcs_seq_R2;
	my $len_mismatch = 0;
	if ($dcs){
		$dcs_seq = generate_dcs($sscs_high, $sscs_low);
		if ($merge_history{$id}){
			unless ($dcs_seq){
				++$dcs_length_mismatches;
				$len_mismatch = 1;
			}
		}else{
			$dcs_seq_R2 = generate_dcs($sscs_high_R2, $sscs_low_R2);
			unless ($dcs_seq and $dcs_seq_R2){
				++$dcs_length_mismatches;
				$len_mismatch = 1;
			}
		}
	}
	
	#calculate region of overlap if read was merged. Then print fasta output
	if ($merge_history{$id}){
		my $L1;
		my $L2;
		my $linker_length = length ($LINKER);
		
		#very rarely there is no matching entry in this hash. Perhaps because the merging process actually changed the barcode. Just skip these.
		unless (exists($premerge_length_HoA{$id})){
			++$lost_start_lengths and next; 
		}
		
		if ($premerge_length_HoA{$id}[0] - $BARCODE_LEN - $linker_length < length($sscs_high)){
			$L1 = $premerge_length_HoA{$id}[0] - $BARCODE_LEN - $linker_length;
		}else{
			$L1 = length($sscs_high);
		}
	
		if ($premerge_length_HoA{$id}[1] - $BARCODE_LEN - $linker_length < length($sscs_high)){
			$L2 = $premerge_length_HoA{$id}[1] - $BARCODE_LEN - $linker_length;
		}else{
			$L2 = length($sscs_high);
		}
	
		my $S1 = length($sscs_high) - $L2 + 1;
		my $S2 = $L1; 
	
		if ($dcs){
			#only print DCS and overlap coordinates for SSCS, if there was no DCS length mismatch.
			unless ($SUPPRESS_SSCS){
				print $FHSM1 ">$id\-$print_strand1 ", scalar (@seqs_high);
				$len_mismatch or print $FHSM1 " $S1 $S2";
				print $FHSM1 "\n$sscs_high\n";
				print $FHSM2 ">$id\-$print_strand2 ", scalar (@seqs_low);
				$len_mismatch or print $FHSM2 " $S1 $S2";
				print $FHSM2 "\n$sscs_low\n";
			}

			unless ($len_mismatch){
				if ($COMBINE_DCS){
					print $FHDC ">$id ", scalar (@seqs_high), "-", scalar (@seqs_low), " $S1 $S2\n$dcs_seq\n";
				}else{
					print $FHDM ">$id ", scalar (@seqs_high), "-", scalar (@seqs_low), " $S1 $S2\n$dcs_seq\n";
				}
			}
		}else{
			$SUPPRESS_SSCS or print $FHSM0 ">$id\-$print_strand1 ", scalar (@seqs_high), " $S1 $S2\n$sscs_high\n";
		}
	
	}else{ #if not merged, print to separate R1 and R2 files
		if ($dcs){
			unless ($SUPPRESS_SSCS){
				print $FHSU1_1 ">$id\-$print_strand1 ", scalar (@seqs_high), "\n$sscs_high\n";		
				print $FHSU1_2 ">$id\-$print_strand2 ", scalar (@seqs_high_R2), "\n$sscs_high_R2\n";		
				print $FHSU2_1 ">$id\-$print_strand2 ", scalar (@seqs_low), "\n$sscs_low\n";		
				print $FHSU2_2 ">$id\-$print_strand1 ", scalar (@seqs_low_R2), "\n$sscs_low_R2\n";
			}
			unless ($len_mismatch){
				if ($COMBINE_DCS){
					print $FHDC  ">$id\[R1\] ", scalar (@seqs_high), "-", scalar (@seqs_low), "\n$dcs_seq\n"; 		
					print $FHDC  ">$id\[R2\] ", scalar (@seqs_high_R2), "-", scalar (@seqs_low_R2), "\n$dcs_seq_R2\n";
				}else{
					print $FHDU_1  ">$id\[R1\] ", scalar (@seqs_high), "-", scalar (@seqs_low), "\n$dcs_seq\n"; 		
					print $FHDU_2  ">$id\[R2\] ", scalar (@seqs_high_R2), "-", scalar (@seqs_low_R2), "\n$dcs_seq_R2\n";				
				}
			}		
		}else{
			unless ($SUPPRESS_SSCS){
				print $FHSU0_1 ">$id\-$print_strand1 ", scalar (@seqs_high), "\n$sscs_high\n";
				print $FHSU0_2 ">$id\-$print_strand2 ", scalar (@seqs_high_R2), "\n$sscs_high_R2\n";
			}		
		}
	}
}


if ($DELETE_INTERMEDIATE_FASTQS){
	unlink "$OUTPUT\_linker_filt1.fq";
	unlink "$OUTPUT\_linker_filt2.fq";
} 

my %min_dcs_size_counts;

foreach my $val1 (keys %double_HoH){
	foreach my $val2 (keys %{$double_HoH{$val1}}){
		$min_dcs_size_counts{min($val1, $val2)} += $double_HoH{$val1}->{$val2};
	}
}

my $FHS = open_output ("$OUTPUT\.family_stats.txt");


print $FHS "Total families: $family_count\n";
print $FHS "Double-stranded families: $double_strands\n";
print $FHS "Single-stranded families: $single_strands\n";

print $FHS "Family Size (SS)\tCount\n";

foreach my $size (sort {$a <=> $b} keys %single_hash){
	print $FHS "$size\t$single_hash{$size}\n";
}


print $FHS "\nMin Family Size (DS)\tCount\n";
foreach my $size (sort {$a <=> $b} keys %min_dcs_size_counts){	
	print $FHS "$size\t$min_dcs_size_counts{$size}\n";
}

print $FHS "\nAB Family Size\tBA Family Size\tCount\n";

foreach my $size_ab (sort {$a <=> $b} keys %double_HoH){	
	foreach my $size_ba (sort {$a <=> $b} keys %{$double_HoH{$size_ab}}){
		print $FHS "$size_ab\t$size_ba\t$double_HoH{$size_ab}->{$size_ba}\n";
	}
}

print "\n" . (localtime) . "\nAnalysis Complete";
$DELETE_INTERMEDIATE_FASTQS and print ". Only the final fasta consensus sequence files have been retained because --delete_intermediate_fastqs was specified.";
print "\n\n";

print "Read counts and filtering:\n";
print "$postlinker_readcount\tInput reads after linker filtering\n";
$DISABLE_N_FILTER or print "$n_bc_count\tFiltered reads that contained one or more Ns in barcode\n";
$DISABLE_REP_FILTER or print "$homopolymer_bc_count\tFiltered reads that contained the same repetitive base for the their entire barcode\n";
$MIN_BAROCODE_QUAL and print "$lowqual_bc_count\tFiltered reads that contained one or more base quality scores below $MIN_BAROCODE_QUAL in barcode\n";
print "$identical_bc_count\tFiltered reads with alpha and beta barcodes that are identical to each other\n";
$MIN_OPTICAL_DIST and print "$opt_count\tFiltered reads that were within $MIN_OPTICAL_DIST pixels of another read with the same barcode\n";
print "$hetero_merge_excludes\tFltered reads that were unmerged when others in the same family were successfully merged by bbmerge\n";
print "$dcs_length_mismatches\tDCS sequences not reported because the two SSCSs differed in length\n";
print "$lost_start_lengths\tFamilies excluded because pre-merged read lengths could not be recovered\n";

#$MIN_TILE_EDGE_DIST and print "$edge_count\tFiltered reads that were within $MIN_TILE_EDGE_DIST pixels of left and/or lower border of tile\n";

#subroutine that takes an array of sequences and returns a consensus sequence
#Called as follows: generate_consensus (\@seqs, $MIN_SSCS, $MIN_AGREE)
sub generate_consensus{

	my $array_ref = shift @_;
	my $count_thresh = shift @_;
	my $agree_thresh = shift @_;
	my @seq_array = @{$array_ref};
	
	my $max_len = 0;
	
	foreach (@seq_array){
		length ($_) > $max_len and $max_len = length ($_);
	}
	
	my $consensus_seq;
	for (my $i = 0; $i < $max_len; ++$i){
		my $site_read_count = 0;
		my %site_hash;
		foreach (@seq_array){
			if (length ($_) > $i){
				my $nuc = substr($_, $i, 1);
				unless ($nuc eq "N"){
					++$site_read_count;
					++$site_hash{$nuc}			
				} 			
			}
		}	
		if ($site_read_count < $count_thresh){
			$consensus_seq .= "N";
		}else{
			my @nuc_keys = sort { $site_hash{$b} <=> $site_hash{$a} } keys %site_hash;
			if ($site_hash{$nuc_keys[0]} / $site_read_count >= $agree_thresh){
				$consensus_seq .= $nuc_keys[0];
			}else{
				$consensus_seq .= "N";		
			}
		}
	}
	return $consensus_seq;
}

#subroutine that takes an two sequences and returns a consensus
#Called as follows: generate_dcs ($seq1, $seq2)
sub generate_dcs{
	my $seq1 = shift @_;
	my $seq2 = shift @_;
	my $dcs_consensus_seq;
	if (length ($seq1) == length ($seq2)){
		if ($seq1 eq $seq2){
			$dcs_consensus_seq = $seq1;
		}else{
			for (my $i = 0; $i < length($seq1); ++$i){
				if (substr($seq1, $i, 1) eq substr($seq2, $i, 1)){
					$dcs_consensus_seq .= substr($seq1, $i, 1);
				}else{
					$dcs_consensus_seq .= "N";
				}
			}
		}
	}
	return $dcs_consensus_seq;
}