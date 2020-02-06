#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IPC::Cmd qw(can_run);
use List::Util qw (min max);
use Bio::SearchIO;
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes SNP output files from the duplexseq pipeline
   and combines them while applying various filters and adding
   information related to kmer counts, contamination databases,
   large repeat sequences, etc.
   
   REQUIRED ARGUMENTS
   
   Input File List
         --input
         Name (and path if not in local directory) of file containing a list
         of SNP output files from duplexseq_sam_to_snps.pl. There should be
         one file name per line, and the path should be included if the files
         are not in the same directory as the SNP files. Alternatively, the
         --input_path option can be used to pre-append the same path to all
         file names.
   
   Output Name
         --output
         Name of output file
   
 
   OPTIONAL ARGUMENTS

   Path to Input SNP Files
         --input_path   
         If the SNP input files are not in the local directory and their full
         path is not specified in the file provided with the --input argument,
         then provide the path with this option. This should not be used if
         the list of file names already includes path information. It will not
         be applied to the location of the input list of filenames itself (that
         path should be provided as part of the --input option if necessary).
              
   Minimum Read Position
         --min_read_pos [default: 1]      
         This filter will exclude variants found near the ends of reads. For
         example, setting this parameter to 11 will effectively remove any 
         variants in the first or last 10 bps of a read. The deafult setting
         of 1 will not filter any variants.

   k-mer Database
         --kmer_db
         Database of k-mer counts as produced kmc_dump or similar program.
         Should be formated as a text file with one k-mer per line. Line
         format should be k-mer sequence followed by one white-space character
         followed by k-mer count. Counts of indentified variants will be
         reported. More than one database can be specified with a comma-
         delimited list. This feature may be useful for identifying
         contaminating sequences, such as apparent mitochondrial variants that
         are actually derived from nuclear sequences (i.e., numts). In this
         scenario, the k-mer database might be generated from a standard
         shotgun library of total-cellular DNA. Calling this option requires
         python3 be installed and in your PATH with that executable name. A 
         PERL-based implementation of k-mer counting is available in the
         initial SAM parsing part of the pipeline, but it is much slower and
         impractical with large k-mer databases.

   k-mer Database Path
         --kmer_db_path [default = ./ (current working directory)]
         Path to location of the datebase file name(s) specified with --kmer_db	

   k-mer Script
         --kmer_script [default: duplexseq_snps_kmercount.py]
         The external script that is called if the --k-mer_db option is
         specified. Provide full path and file name if it is not already in
         your local directory.
         
   k-mer Size
         --kmer_size [default = 39]
         Size of k-mers. If a k-mer database is specified with --kmer_db to
         generate counts, this k-mer size must match the size used to generate
         the database with kmc or similar program.

   Recombination Check File
         --recomb_check
         Fasta file with genome to be used to check for variants that might
         have been created by recombinaion between short, non-identical
         repeats rather than de novo mutations. This would typically be the
         reference fasta file used for mapping. A blast database for this
         fasta must also exist in the same location (generated with
         make_blast_db in the NCBI BLAST+ package). The NCBI BLAST+
         executables must be installed and in your PATH if this option is
         used. It is also requires BioPerl.

   Recombination Check Flank
         --recomb_flank [default = 100]
         The amount of flanking sequence from the reference genome surrounding
         the variant that is used for searching for possible recombinants.
         This is only used if --recomb_check option is specified.

   Recombination Check E-value
         --recomb_evalue [default = 1e-10]
         The e-value threshold applied for BLAST searches in checking for
         variants resulting from recombination if the --recomb_check option is
         specified.

   Contamination Check File
         --contam_check
         Fasta file with sequences representing potential sources of
         contamination. If this option is specified, variants will be flagged
         if the entire duplex read has a perfect match to this database. A
         blast database for this fasta must also exist in the same location
         (generated with make_blast_db in the NCBI BLAST+ package). The NCBI
         BLAST+ executables must be installed and in your PATH if this option
         is used.

   Repeat File
         --repeat_file
         File describing repeat content in the genome to identify positions as
         repetitive and map them identical positions elsewhere in the genome.
         A tab-delimited list of start/end coordinates and orientation info.

   Minimum Repeat Length
         --min_rep_len [default = 100]
         The minimum length of a perfect repeat to be used in repeat mapping.
                
";

our $INPUT;
our $OUTPUT;
our $INPUT_PATH;
our $MIN_READ_POS = 1;
our $KMER_DB;
our $KMER_DB_PATH = './';
our $KMER_SCRIPT = "duplexseq_snps_kmercount.py";
our $KMER_SIZE = 39;
our $RECOMB_CHECK;
our $RECOMB_FLANK = 100;
our $RECOMB_EVALUE = 1e-10;
our $CONTAM_CHECK;
our $REPEAT_FILE;
our $MIN_REP_LEN = 100;


##Print out start time and command line call of Perl script
print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";


GetOptions(
    'input=s'  => \$INPUT,
    'output=s'  => \$OUTPUT,
    'input_path=s'  => \$INPUT_PATH,
    'min_read_pos=i'  => \$MIN_READ_POS,    
    'kmer_db=s'  => \$KMER_DB,
    'kmer_db_path=s'  => \$KMER_DB_PATH,
    'kmer_script=s'  => \$KMER_SCRIPT,
    'kmer_size=i'  => \$KMER_SIZE,
    'recomb_check=s'  => \$RECOMB_CHECK,
    'recomb_flank=i'  => \$RECOMB_FLANK,
    'recomb_evalue=i'  => \$RECOMB_EVALUE,
    'contam_check=s'  => \$CONTAM_CHECK,
    'repeat_file=s'  => \$REPEAT_FILE,
    'min_rep_len=i'  => \$MIN_REP_LEN 
);


##check for proper options specification
$INPUT or die ("\n$usage\n\nERROR: Must provide a file with input filenames with --input.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide a output filename with --output.\n\n");
if ($RECOMB_CHECK or $CONTAM_CHECK){
	can_run ("blastn") or die ("\n$usage\n\nERROR: If the --contam_check or --recomb_check options are specified, blastn must be installed and in your path.\n\n");
}
if ($KMER_DB){
	can_run ($KMER_SCRIPT) or die ("\n$usage\n\nERROR: Could not find $KMER_SCRIPT\. If the --kmer_db option is specified, an external python script will be called. Provide file name and path with --kmer_script option.\n\n");
	substr ($KMER_DB_PATH,-1) eq '/' or $KMER_DB_PATH .= '/';
	my @kmer_list = split (/\,/, $KMER_DB);
	foreach (@kmer_list){
		my $kmer_test = $KMER_DB_PATH . $_;
		-e $kmer_test or die ("\n$usage\n\nERROR: Could not find $kmer_test. Could not find the database specified with the --kmer_db and --kmer_db_path options.\n\n");
	}
	can_run ("python3") or die ("\n$usage\n\nERROR: If the --kmer_db option is specified, python3 must be installed and in your path.\n\n");
}
if ($REPEAT_FILE){
	-e $REPEAT_FILE or die ("\n$usage\n\nERROR: Could not find $REPEAT_FILE\. Could not find the file specified with the --repeat_file option.\n\n");
}

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

##Open temporary output files that will be used to store lines after filtering for read ends and converting repeat-related coordinates. Printer header info to the second of these temp files
my $FHOT = open_output ($OUTPUT . "_TEMP");
my $FHOT2 = open_output($OUTPUT . "_TEMP2");
print $FHOT2 "File\tID\tFamily Size\tReference\tRef Position\tRead Position\tMerge Coverage\tMap Orientation\tRef\tVariant\tRef k-mer\tAlt k-mer\tRead Sequence";
if ($REPEAT_FILE){
	print $FHOT2 "\tRepeat Region\tRepeat Orientation";
}
if ($RECOMB_CHECK){
	print $FHOT2 "\tRecomb Match";
}
if ($CONTAM_CHECK){
	print $FHOT2 "\tContam Match";
}

print $FHOT2 "\n";


#generate summary of repeats if --repeat_file is specified
#define data structure. Hash of hashes of arrays. Key structure: sequence name -> nucleotide position. Values arrays with first element is mapped repeat position and second element is orientation (1 is direct, -1 is inverted).
#note that there is the possibility of multi-copy repeats that are nested or overlapping. So all repeats will be mapped to their earliest position in the genome (including their own position if they are the first occurring repeat copy in the genome.)
#The reported orientation is relative to the copy to which the position is mapped, not necessarily the other repeat copy in the genome. Repeat positions mapped to themseleves will always be reported as forward even if the other copy in the genome is an IR.

my %reps_HoHoA;
if ($REPEAT_FILE){
	my @rep_lines = file_to_array ($REPEAT_FILE);
	foreach (@rep_lines){
		$_ =~ /^\s+$/ and next;
		chomp $_;
		my @sl = split (/\t/, $_);
		$sl[2] - $sl[1] == $sl[4] - $sl[3] or die ("\nERROR: The pair of repeats in the following line from $REPEAT_FILE do not have equal lengths. Repeats in this files should be perfect matches and thus have the same lengths:\n\n$_\n\n");
		$sl[1] < $sl[3] or next;
		$sl[2] - $sl[1] + 1 >= $MIN_REP_LEN or next;
		
		my $strand = $sl[5];
		$strand == 1 or $strand == -1 or die ("\nERROR: Could not interpret orientationinformation in the following line from $REPEAT_FILE:\n\n$_\n\n");
		
		my $pos2;
		if ($strand == 1){
			$pos2 = $sl[3];
		}else{
			$pos2 = $sl[4];
		}
		for (my $pos1 = $sl[1]; $pos1 <= $sl[2]; ++$pos1){

			if (exists ($reps_HoHoA{$sl[0]}->{$pos1})){
				if (min($pos1, $pos2) < $reps_HoHoA{$sl[0]}->{$pos1}[0]){
					$reps_HoHoA{$sl[0]}->{$pos1}[0] = min($pos1, $pos2);
					if ($pos1 < $pos2){
						$reps_HoHoA{$sl[0]}->{$pos1}[1] = 1;
					}else{
						$reps_HoHoA{$sl[0]}->{$pos1}[1] = $strand;
					}
				}
			}else{
				$reps_HoHoA{$sl[0]}->{$pos1}[0] = min($pos1, $pos2);
				if ($pos1 < $pos2){
					$reps_HoHoA{$sl[0]}->{$pos1}[1] = 1;
				}else{
					$reps_HoHoA{$sl[0]}->{$pos1}[1] = $strand;
				}			
			}

			if (exists ($reps_HoHoA{$sl[0]}->{$pos2})){
				if (min($pos1, $pos2) < $reps_HoHoA{$sl[0]}->{$pos2}[0]){
					$reps_HoHoA{$sl[0]}->{$pos2}[0] = min($pos1, $pos2);
					if ($pos1 < $pos2){
						$reps_HoHoA{$sl[0]}->{$pos2}[1] = $strand;
					}else{
						$reps_HoHoA{$sl[0]}->{$pos2}[1] = 1;
					}
				}
			}else{
				$reps_HoHoA{$sl[0]}->{$pos2}[0] = min($pos1, $pos2);
				if ($pos1 < $pos2){
					$reps_HoHoA{$sl[0]}->{$pos2}[1] = $strand;
				}else{
					$reps_HoHoA{$sl[0]}->{$pos2}[1] = 1;
				}			
			}

			
			if ($strand == 1){
				++$pos2;
			}else{
				--$pos2;
			}			
		}	
	}	
}


##define data structure to store recombination check results. Hash of hash of hashes. Key structure: sequence_name-->Position-->Variant. Value is the concatenation of sequence_name, hit_position, percent id, and hit length.
my %recomb_HoHoH;
##define data structure to store contamination check results. Hash of hash of hashes. Key structure: sequence_name-->Position-->Variant. Value is the concatenation of sequence_name and hit_position.
my %contam_HoHoH;




##loop through every provided variant file
foreach (@files){
	my @variant_lines = file_to_array($_);
	shift @variant_lines;
	#loop through every variant line in the provided file
	foreach my $line (@variant_lines){
		chomp $line;
		my @sl = split (/\t/, $line);
		#skip variants that are too close to read ends
		$sl[4] >= $MIN_READ_POS or next;

		print $FHOT "$_\t$line\n";
		
		#check for variants resulting from recombination between short non-identical repeats
		if ($RECOMB_CHECK){
			my %recomb_fasta = fasta2hash ($RECOMB_CHECK);

			my $alt_allele = $sl[8];
			if ($sl[6] eq 'R'){
				$alt_allele = revcom ($alt_allele);
			}
			exists ($recomb_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele}) and next;
			my $left_end = max (0, $sl[3] - $RECOMB_FLANK - 1);
			my $right_end = min (length ($recomb_fasta{$sl[2]}), $sl[3] + $RECOMB_FLANK - 1);
			my $query = substr ($recomb_fasta{$sl[2]}, $left_end, $right_end - $left_end + 1);
			my $temp_recomb_blast_query = $OUTPUT . "_TEMPRECOMBQUERY";
			my $temp_recomb_blast_output = $OUTPUT . "_TEMPRECOMBOUTPUT";
			system ("echo \>query > $temp_recomb_blast_query; echo $query >> $temp_recomb_blast_query");
			system ("blastn -task blastn -evalue $RECOMB_EVALUE -db $RECOMB_CHECK -query $temp_recomb_blast_query -out $temp_recomb_blast_output");
	
			my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => $temp_recomb_blast_output);

			while( my $result_obj = $SearchIO_obj->next_result ) {
				while (my $hit_obj = $result_obj->next_hit ) {
					my $hit_name = $hit_obj->name;
					while (my $hsp_obj = $hit_obj->next_hsp){
						exists ($recomb_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele}) and last;
						my $hit_string = $hsp_obj->hit_string;
						my $query_string = $hsp_obj->query_string;
						my $query_start = $hsp_obj->start('query');
						my $hit_start = $hsp_obj->start('hit');
						my $hit_end = $hsp_obj->end('hit');
						my $strand = $hsp_obj->strand('hit');
						my $target_pos = min ($RECOMB_FLANK + 1, $sl[3]);
						my $ref_coord = $query_start;
						my $hit_coord;
						if ($strand == 1){
							$hit_coord = $hit_start;				
						}else{
							$hit_coord = $hit_end;
						}
						for (my $i = 0; $i < length ($query_string); ++$i){
							substr ($query_string, $i, 1) eq '-' and next;
							if ($ref_coord == $target_pos){
								if (uc(substr($hit_string, $i, 1)) eq $alt_allele){
									$recomb_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele} = $hit_name . ';' . $hit_coord . ';' . $hsp_obj->percent_identity . ';' . $hsp_obj->length('hit');
								}
								last;
							}
							++$ref_coord;
							unless (substr ($hit_string, $i, 1) eq '-'){
								if ($strand == 1){
									++$hit_coord;
								}else{
									--$hit_coord;
								}
							}
						}			
					}
					exists ($recomb_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele}) or $recomb_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele} = "NA";	
				}
			}
			unlink ($temp_recomb_blast_query);
			unlink ($temp_recomb_blast_output);
		}
		
		#check for variants that are exact matches to a contamination database
		if ($CONTAM_CHECK){
			my $alt_allele = $sl[8];
			if ($sl[6] eq 'R'){
				$alt_allele = revcom ($alt_allele);
			}
			exists ($contam_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele}) and next;
			my $temp_contam_blast_query = $OUTPUT . "_TEMPCONTAMQUERY";
			my $temp_contam_blast_output = $OUTPUT . "_TEMPCONTAMOUTPUT";
			system ("echo \>query > $temp_contam_blast_query; echo $sl[-1] >> $temp_contam_blast_query");
			system ("blastn -evalue 1e-6 -db $CONTAM_CHECK -query $temp_contam_blast_query -out $temp_contam_blast_output");

			my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => $temp_contam_blast_output);

			my $result_obj = $SearchIO_obj->next_result;
			if (my $hit_obj = $result_obj->next_hit ) {
				my $hsp_obj = $hit_obj->next_hsp;
				if ($hsp_obj->frac_identical == 1 && $hsp_obj->length('query') == $result_obj->query_length){
					my $hit_start = $hsp_obj->start('hit');
					my $hit_end = $hsp_obj->end('hit');
					my $strand = $hsp_obj->strand('hit');
					my $hit_coord;
					if ($strand == 1){
						$hit_coord = $hit_start + $sl[4] - 1;				
					}else{
						$hit_coord = $hit_end - $sl[4] + 1;
					}
					$contam_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele} = $hit_obj->name . ';' . $hit_coord;					
				}else{
					$contam_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele} = "NA";
				}
			}else{
				$contam_HoHoH{$sl[2]}->{$sl[3]}->{$alt_allele} = "NA";
			}
			unlink ($temp_contam_blast_query);
			unlink ($temp_contam_blast_output);		
		}
	}
}
close $FHOT;


##Process the data stored for repeat mapping, recombination check, and contamination check. Storing in a second temp file.
my @merged_snp_lines = file_to_array($OUTPUT . "_TEMP");
foreach (@merged_snp_lines){

	chomp $_;
	print $FHOT2 $_;
	my @sl = split (/\t/, $_);

	if ($REPEAT_FILE){
		
		if (exists($reps_HoHoA{$sl[3]}->{$sl[4]})){
			print $FHOT2 "\t", $reps_HoHoA{$sl[3]}->{$sl[4]}[0], "\t", $reps_HoHoA{$sl[3]}->{$sl[4]}[1];
		}else{
			print $FHOT2 "\t.\t.";
		}
	
	}
	
	if ($RECOMB_CHECK){
		my $alt_allele = $sl[9];
		if ($sl[7] eq 'R'){
			$alt_allele = revcom ($alt_allele);
		}
		if (exists ($recomb_HoHoH{$sl[3]}->{$sl[4]}->{$alt_allele})){
			print $FHOT2 "\t", $recomb_HoHoH{$sl[3]}->{$sl[4]}->{$alt_allele};
		}else {
			die ("\nERROR: Did not record blast result for $sl[3] $sl[4] $alt_allele\.\n\n");
		}
	}

	if ($CONTAM_CHECK){
		my $alt_allele = $sl[9];
		if ($sl[7] eq 'R'){
			$alt_allele = revcom ($alt_allele);
		}
		if (exists ($contam_HoHoH{$sl[3]}->{$sl[4]}->{$alt_allele})){
			print $FHOT2 "\t", $contam_HoHoH{$sl[3]}->{$sl[4]}->{$alt_allele};
		}else {
			die ("\nERROR: Did not record blast result for $sl[3] $sl[4] $alt_allele\.\n\n");
		}
	}
	print $FHOT2 "\n";
}
close $FHOT2;

if ($KMER_DB){
	my $temp_data = $OUTPUT . "_TEMP2";
	my $total_col_num = 13;
	$REPEAT_FILE and $total_col_num += 2;
	$RECOMB_CHECK and ++$total_col_num;
	$CONTAM_CHECK and ++$total_col_num;
	
	system("$KMER_SCRIPT $temp_data $KMER_DB_PATH $KMER_DB $OUTPUT $total_col_num");
	
}else{
	my $temp_data = $OUTPUT . "_TEMP2";
	system ("cat $temp_data >> $OUTPUT");
}

unlink ($OUTPUT . "_TEMP");
unlink ($OUTPUT . "_TEMP2");

print "\n" . (localtime) . "\nAnalysis Complete\n\n";
