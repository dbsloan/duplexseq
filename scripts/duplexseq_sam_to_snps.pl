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
   SSCS or DCS output from duplex sequencing. It returns a summary of SNP calls.
   
   REQUIRED ARGUMENTS
   
   SAM File
         --sam
         File name of SAM output from bowtie2 mapping of DCS/SSCS sequences
         against a reference genome.
   
   Consensus Sequences Fasta File
         --seqs
         Fasta file containing DCS or SSCS reads from duplex sequencing.

   Reference Genome Fasta
         --ref
         Fasta file containing genome sequence(s) used as reference for read
         mapping.

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
         count. Counts of indentified variants will be reported. More than one 
         database can be specified with a comma-delimited list. This feature may
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
our $REF;
our $OUTPUT;
our $EXCLUDE_REFS;
our $KMER_DB;
our $KMER_SIZE=39;


GetOptions(
    'sam=s'  => \$SAM,
    'seqs=s'  => \$SEQS,
    'ref=s'  => \$REF,
    'output=s'  => \$OUTPUT,
    'exclude_refs=s'  => \$EXCLUDE_REFS,
    'kmer_db=s'  => \$KMER_DB,
    'kmer_size=i'  => \$KMER_SIZE
);


print "\n" . (localtime) . "\nRunning $0 with the following command:\n", qx/ps -o args $$/, "\n";

$SAM or die ("\n$usage\n\nERROR: Must provide a SAM input filename with --sam.\n\n");
$SEQS or die ("\n$usage\n\nERROR: Must provide a DCS or SSCS fasta input filename with --seqs.\n\n");
$REF or die ("\n$usage\n\nERROR: Must provide a reference fasta input filename with --ref.\n\n");
$OUTPUT or die ("\n$usage\n\nERROR: Must provide a base name for output files with --output.\n\n");

#populate hashes of kmer counts if kmer database(s) provided.
my @kmer_AoH;
my @kmer_db_names;
if ($KMER_DB){
	my @split_string = split (/\,/, $KMER_DB);
	my $db_count = 0;
	foreach (@split_string){
		push (@kmer_db_names, $_);
		-e $_ or die ("\n$usage\n\nERROR: Specified k-mer database ($_\) does not exist.\n\n");
		my $FH_KMER = open_file ($_);
		my $kmer_first_line = <$FH_KMER>;
		chomp $kmer_first_line;
		my @sl = split (/\s/, $kmer_first_line);
		my $kmer_test_length = length($sl[0]);
		$kmer_test_length == $KMER_SIZE or die ("\n$usage\n\nERROR: The k-mer size in the specified database ($_\) is $kmer_test_length\, which does not match specified k-mer size ($KMER_SIZE). Change the k-mer size to $kmer_test_length with --kmer_size or regenerate the k-mer database with the desired size.\n\n");
		$kmer_AoH[$db_count]{$sl[0]} = $sl[1];
		while (my $line = <$FH_KMER>){
			chomp $line;
			my @sl2 = split (/\s/, $line);
			$kmer_AoH[$db_count]{$sl2[0]} = $sl2[1];
		}
		++$db_count;
	}
}
my $kmer_left_flank = floor ($KMER_SIZE / 2);
my $kmer_right_flank;
if ($kmer_left_flank == $KMER_SIZE / 2){
	$kmer_right_flank = $kmer_left_flank - 1;
}else{
	$kmer_right_flank = $kmer_left_flank;
}



##read in fasta reference genome and store in hash. replace it with version that removes anything after white space in header (bowtie2 truncates similarly);
my %ref_fasta = fasta2hash($REF);
my %ref_fasta2;
foreach (keys %ref_fasta){
	my @sk = split (/\s/, $_);
	exists ($ref_fasta2{$sk[0]}) and die ("\nERROR: Multiple sequences in reference with the same truncated name: $sk[0]\n\n");
	$ref_fasta2{$sk[0]} = $ref_fasta{$_};
}
%ref_fasta = %ref_fasta2;

##read in all full headers from the input DCS or SSCS fasta file
my ($headersRef, $seqsRef) = get_fasta_names_and_seqs($SEQS);
my @headers = @{$headersRef};

##build a hash of arrays to store info on each read family from fasta headers. Assing values of 0 to merge coordinates if it's an unmerged (i.e., R1 or R2) sequence.
my %readfamily_HoA;
foreach (@headers){
	chomp $_;
	my @sl = split (/\ /, $_);
	$readfamily_HoA{$sl[0]}[0] = $sl[1];
	if ($sl[2]){
		$readfamily_HoA{$sl[0]}[1] = $sl[2];
		$readfamily_HoA{$sl[0]}[2] = $sl[3];	
	}else{
		$readfamily_HoA{$sl[0]}[1] = 0;
		$readfamily_HoA{$sl[0]}[2] = 0;		
	}
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
my $FHO1 = open_output ("$OUTPUT\.cov_stats.txt");
my $FHO2 = open_output ("$OUTPUT\.snps.txt");
print $FHO2 "ID\tFamily Size\tReference\tRef Position\tRead Position\tMerge Coverage\tMap Orientation\tRef\tVariant\tRef k-mer\tAlt k-mer";
if ($KMER_DB){
	foreach (@kmer_db_names){
		print $FHO2 "\tRef k-mer $_\tAlt k-mer $_";
	}
}
print $FHO2 "\tRead Sequence\n";



##loop through all lines in SAM file
my $FH = open_file ($SAM);
#Both of the following ignore reads that have any indels.
#And they include reads with >1 SNP (which are not considered in variant calling).
#And they include single SNPs in match counts.
#So they are approximations, but these differences should be minor (trivial) in most cases.
my %total_cov_hash; #data structure to store raw mapping coverage in bp for each reference seq
#The following is a data structure to store detailed mapping coverage. Key structure is as follows:
#1: Referene sequence name
#2: Distance to flanking sequence
#3: Double coverage
#4: Reference base
my %match_HoHoHoH;

while (<$FH>){

	$_ =~ /^\@/ and next;
	my @sl = split (/\t/, $_);

	##skip lines that do not map
	$sl[2] eq '*' and next;

	##skip lines that map to reference sequences on the exclude list
	exists ($excludes{$sl[2]}) and next;

	##look at lines with gap-free mapping
	if ($sl[5] =~ /^(\d+)M$/){
		
		##Store data for matching
		my $match_len = $1;
		$total_cov_hash{$sl[2]} += $match_len;
		
		
		for (my $i = 0; $i < length($sl[9]); ++$i){	
			my $dist_to_end = min ($i+1, length($sl[9]) - $i);
			my $reference_hit = $sl[2];
			my $double_cov = 0;
			my $match_base;
			if ($sl[1] == 0){
				if ($readfamily_HoA{$sl[0]}[1]){
					if ($i+1 >= $readfamily_HoA{$sl[0]}[1] and $i+1 <= $readfamily_HoA{$sl[0]}[2]){
						$double_cov = 1;
					}
				}	
				#check if this sequence has been reversed during data processing (only applies to SSCS sequences)
				if (substr($sl[0],-2,1) eq '-' and substr($sl[0],-1) eq 'R'){
					$match_base =  revcom(substr($sl[9], $i, 1));
				}else{
					$match_base =  substr($sl[9], $i, 1);				
				}
			}elsif ($sl[1]==16){
				if ($readfamily_HoA{$sl[0]}[1]){
					if (length($sl[9]) - $i >= $readfamily_HoA{$sl[0]}[1] and length($sl[9]) - $i <= $readfamily_HoA{$sl[0]}[2]){
						$double_cov = 1;
					}
				}
				#check if this sequence has been reversed during data processing (only applies to SSCS sequences)
				if (substr($sl[0],-2,1) eq '-' and substr($sl[0],-1) eq 'R'){
					$match_base =  substr($sl[9], $i, 1);
				}else{
					$match_base =  revcom(substr($sl[9], $i, 1));				
				}
			}else{
				die ("\nERROR: bit string not equal to 0 or 16:\n $_\n\n");
			}
			++$match_HoHoHoH{$reference_hit}->{$dist_to_end}->{$double_cov}->{$match_base};
		}		
		##look at reads that show a single substitution relative to reference for variant calling
		if ($_ =~ /XM\:i\:1\s/){
			
			my $read = $sl[9];
			my $read_len = length($read);
			my $start = $sl[3];
			my $refname = $sl[2];
			my $strand = $sl[1];
			my $refseq = substr ($ref_fasta{$refname}, $start - 1, $read_len);
			
			my $kmer = "";
			my $kmer_ref = "";
			
			if ($read_len < $kmer_left_flank + $kmer_right_flank + 1){
				$kmer = "NA";
				$kmer_ref = "NA";
			}
			
			if ($strand == 0){
				for (my $i = 0; $i < $read_len; ++$i){					
					unless (substr($read, $i, 1) eq substr($refseq, $i, 1)){
						my $pos = min ($i+1, $read_len - $i);
						##if these are merged reads, check whether variant is in overlap region
						my $double = 0;
						if ($readfamily_HoA{$sl[0]}[1]){
							if ($i+1 >= $readfamily_HoA{$sl[0]}[1] and $i+1 <= $readfamily_HoA{$sl[0]}[2]){
								$double = 1;
							}
						}
						print $FHO2 "$sl[0]\t", $readfamily_HoA{$sl[0]}[0], "\t$refname\t", $start + $i, "\t$pos\t$double\tF\t";
						#check if this sequence has been reversed during data processing (only applies to SSCS sequences)
						if (substr($sl[0],-2,1) eq '-' and substr($sl[0],-1) eq 'R'){
							print $FHO2 revcom (substr($refseq, $i, 1)), "\t", revcom(substr($read, $i, 1));
						}else{
							print $FHO2 substr($refseq, $i, 1), "\t", substr($read, $i, 1);
						}
						unless ($kmer eq "NA"){
							if ($i >= $kmer_left_flank and $read_len - $i - 1 >= $kmer_right_flank){
								$kmer = substr ($read, $i - $kmer_left_flank, $KMER_SIZE);
								$kmer_ref = substr ($refseq, $i - $kmer_left_flank, $KMER_SIZE);
							}elsif($i < $kmer_left_flank){
								$kmer = substr ($read, 0, $KMER_SIZE);
								$kmer_ref = substr ($refseq, 0, $KMER_SIZE);
							}else{
								$kmer = substr ($read, $read_len - $KMER_SIZE);
								$kmer_ref = substr ($refseq, $read_len - $KMER_SIZE);
							}
						}
					}
				}
			}elsif ($strand==16){
				for (my $i = 0; $i < $read_len; ++$i){
					unless (substr($read, $i, 1) eq substr($refseq, $i, 1)){
						my $pos = min ($i+1, $read_len - $i);
						my $double = 0;
						if ($readfamily_HoA{$sl[0]}[1]){
							if ($read_len - $i >= $readfamily_HoA{$sl[0]}[1] and $read_len - $i <= $readfamily_HoA{$sl[0]}[2]){
								$double = 1;
							}
						}
						print $FHO2 "$sl[0]\t", $readfamily_HoA{$sl[0]}[0], "\t$refname\t", $start + $i, "\t$pos\t$double\tR\t";
						#check if this sequence has been reversed during data processing (only applies to SSCS sequences)
						if (substr($sl[0],-2,1) eq '-' and substr($sl[0],-1) eq 'R'){
							print $FHO2 substr($refseq, $i, 1), "\t", substr($read, $i, 1);
						}else{
							print $FHO2 revcom (substr($refseq, $i, 1)), "\t", revcom (substr($read, $i, 1));
						}
						unless ($kmer eq "NA"){
							if ($i >= $kmer_left_flank and $read_len - $i - 1 >= $kmer_right_flank){
								$kmer = substr ($read, $i - $kmer_left_flank, $KMER_SIZE);
								$kmer_ref = substr ($refseq, $i - $kmer_left_flank, $KMER_SIZE);
							}elsif($i < $kmer_left_flank){
								$kmer = substr ($read, 0, $KMER_SIZE);
								$kmer_ref = substr ($refseq, 0, $KMER_SIZE);
							}else{
								$kmer = substr ($read, $read_len - $KMER_SIZE);
								$kmer_ref = substr ($refseq, $read_len - $KMER_SIZE);
							}
						}
					}
				}	
			}else{
				die ("\nERROR: bit string not equal to 0 or 16:\n $_\n\n");
			}
			if ($kmer eq "NA" or (not defined ($KMER_DB))){
				print $FHO2 "\t$kmer_ref\t$kmer\t$sl[9]\n";
			}else{
				print $FHO2 "\t$kmer_ref\t$kmer";
				foreach my $db (@kmer_AoH){
				
					my %kmer_hash = %{$db};
				
					my $kmer_count = 0;
					my $kmer_ref_count = 0;
					if (exists ($kmer_hash{$kmer})){
						$kmer_count += $kmer_hash{$kmer};
					}
					if (exists ($kmer_hash{revcom($kmer)})){
						$kmer_count += $kmer_hash{revcom($kmer)};
					}
					if (exists ($kmer_hash{$kmer_ref})){
						$kmer_ref_count += $kmer_hash{$kmer_ref};
					}
					if (exists ($kmer_hash{revcom($kmer_ref)})){
						$kmer_ref_count += $kmer_hash{revcom($kmer_ref)};
					}
					print $FHO2 "\t$kmer_ref_count\t$kmer_count";
				}
				print $FHO2 "\t$sl[9]\n";
			}		
		}	
	}	
}

print $FHO1 "Approximate total coverage by reference sequence...\n\n";
print $FHO1 "Sequence\tCoverage (bp)\n";
foreach (sort keys %total_cov_hash){
	print $FHO1 "$_\t$total_cov_hash{$_}\n";
}
print $FHO1 "Detailed matching coverage information (see notes in script about approximations)...\n\n";
print $FHO1 "Sequence\tRead Position\tMerge Coverage\tRef\tCoverage (bp)\n";

foreach my $seq (sort keys %match_HoHoHoH){
	foreach my $pos (sort {$a <=> $b} keys %{$match_HoHoHoH{$seq}}){
		foreach my $double (sort keys %{$match_HoHoHoH{$seq}->{$pos}}){
			foreach my $ref (sort keys %{$match_HoHoHoH{$seq}->{$pos}->{$double}}){
				print $FHO1 "$seq\t$pos\t$double\t$ref\t$match_HoHoHoH{$seq}->{$pos}->{$double}->{$ref}\n";
			}	
		}
	}
}


print "\n" . (localtime) . "\nAnalysis with $0 complete.\n";
