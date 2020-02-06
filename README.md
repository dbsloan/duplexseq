# duplexseq
Scripts for analysis of duplex sequencing data


## Overview

The scripts in this repository are used for processing of raw reads generated by [the duplex sequencing method](https://www.nature.com/articles/nprot.2014.170) and subsequent variant calling and characterization.

## Dependencies

- [BBMerge](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [sloan.pm Perl module](https://github.com/dbsloan/perl_modules)
- Perl
- Python3 (if using a k-mer database to screen variants with `duplexseq_sam_to_indels.pl` or `duplexseq_sam_to_snps.pl`)

## Pipeline Summary

**Raw read processing and filtering of read family consensus sequence**
- `duplexseq_process_reads.pl`
- `duplexseq_family_filter_and_trim.pl`

**Mapping to reference genome**
- `bowtie2`

**Raw variant calling**
- `duplexseq_sam_to_indels.pl`
- `duplexseq_sam_to_snps.pl`

**Aggregating data from multiple libraries, variant filtering, and summary**
- `duplexseq_indels_post_processing.pl`
- `duplexseq_snps_post_processing.pl`
- `duplexseq_coverage_summary.pl`

**Other scripts**
- `duplexseq_batchscripts.pl` (generates a shell script for each input library to automate the first 3 steps above under our standard parameters)
- `duplexseq_indels_kmercount.py` and `duplexseq_snps_kmercount.py` (accessory scripts called during SNP and indel post-processing if using a k-mer database)
- `duplexseq_snp_classification.pl` (functional annotation of filtered SNP calls)

## duplexseq_batchscripts.pl

This script is called with a single input:

`perl duplexseq_batchscripts.pl input_files.txt`

Where `input_files.txt` is a tab-delimited file with one row for each input library and four columns. The first column is a library name that will be applied to all output files. The second and third columns are read 1 and read 2 fastq files, respectively (gzipped is ok). The fourth column is the name of a fasta file with a corresponding bowtie2 database to be used in mapping. See `example_files/input_files.txt`. The output is one shell script for each input library, which are formatted for submission to a standard Slurm-based Linux queuing system and contain the commands to automate read processing, filtering, mapping, and raw variant calling.

## duplexseq_process_reads.pl

This script processes raw fasta inputs to generate single-stranded consensus sequence (SSCS) families and double-stranded consensus sequence (DCS) families. It can be run directly or with the scripts generated by `duplexseq_batchscripts.pl` (see above). Usage instructions for calling it directly are as follows.

`perl duplexseq_process_reads.pl [options/arguments]`

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
                 



