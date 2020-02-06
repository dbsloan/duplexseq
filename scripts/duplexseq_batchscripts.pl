#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 input_file_summary\n\n";

my $file = shift or die ($usage);
my @file_lines = file_to_array ($file);

my $file_count = 0;
foreach (@file_lines){
	$_ =~ /^\s+$/ and next;
	++$file_count;
	chomp $_;
	my @sl = split (/\t/, $_);
	
	my $FH = open_output("$sl[0]\.sh");
	
print $FH "\#!/bin/bash

\#SBATCH --time=10000:00:00
\#SBATCH --job-name=$sl[0]
\#SBATCH --error=$sl[0]\.stderr
\#SBATCH --output=$sl[0]\.stdout

duplexseq_process_reads.pl --r1_fastq=$sl[1] --r2_fastq=$sl[2] --output=$sl[0] --cutadapt_cores=1 --delete_intermediate_fastqs --suppress_sscs --combine_dcs > $sl[0]\.log 2>&1

duplexseq_family_filter_and_trim.pl --input=$sl[0]\.DCS.fas --output=$sl[0]\.DCS.filt.fas --discard_n_reads >> $sl[0]\.log 2>&1

bowtie2 -f -x $sl[3] -U $sl[0]\.DCS.filt.fas -S $sl[0]\.DCS.filt.sam >> $sl[0]\.log 2>&1

duplexseq_sam_to_snps.pl --sam=$sl[0]\.DCS.filt.sam --ref=$sl[3] --seqs=$sl[0]\.DCS.filt.fas --output=$sl[0] >> $sl[0]\.log 2>&1

duplexseq_sam_to_indels.pl --sam=$sl[0]\.DCS.filt.sam --seqs=$sl[0]\.DCS.filt.fas --output=$sl[0] >> $sl[0]\.log 2>&1
";
	close $FH;
	
}



