#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: perl $0 variant_file coord_map\n\n";

my $variant_file = shift or die ($usage);
my $coordmap_file = shift or die ($usage);

my %coordmap_HoH;
my %gene_HoH;

my $FH = open_file($coordmap_file);

my $first_line = <$FH>;

while (<$FH>){
	chomp $_;
	my @sl = split (/\t/, $_);
	$coordmap_HoH{$sl[0]}->{"ref_nuc"} = $sl[1];
	$coordmap_HoH{$sl[0]}->{"type"} = $sl[2];
	$coordmap_HoH{$sl[0]}->{"gene"} = $sl[3];
	$coordmap_HoH{$sl[0]}->{"strand"} = $sl[4];
	$coordmap_HoH{$sl[0]}->{"gene_pos"} = $sl[5];
	$coordmap_HoH{$sl[0]}->{"gene_nuc"} = $sl[6];
	
	if ($sl[2] eq "CDS"){
		$gene_HoH{$sl[3]}->{$sl[5]} = $sl[6];
	}
}

close $FH;

my $FH2 = open_file($variant_file);

$first_line = <$FH2>;
chomp $first_line;

print "$first_line\tSubType\tMutationContext\tRegion\tGene\tGeneStrand\tAA_sub\n";



while (<$FH2>){
	
	chomp $_;
	print $_, "\t";
	my @sl = split (/\t/, $_);
	
	if ($sl[8] eq 'A' or $sl[8] eq 'G'){
		print "$sl[8]\>$sl[9]\t";
		if ($sl[7] eq 'F'){
			print $coordmap_HoH{$sl[4] - 1}->{"ref_nuc"}, 'N', $coordmap_HoH{$sl[4] + 1}->{"ref_nuc"}, "\t";
		}elsif($sl[7] eq 'R'){
			print revcom ($coordmap_HoH{$sl[4] + 1}->{"ref_nuc"}), 'N', revcom ($coordmap_HoH{$sl[4] - 1}->{"ref_nuc"}), "\t";		
		}else{
			die ("\nERROR: Unexpected orientation $sl[7] in following line:\n$_\n\n");
		}
	}elsif ($sl[8] eq 'C' or $sl[8] eq 'T'){
		print revcom($sl[8]), "\>", revcom ($sl[9]), "\t";	
		if ($sl[7] eq 'F'){
			print revcom ($coordmap_HoH{$sl[4] + 1}->{"ref_nuc"}), 'N', revcom ($coordmap_HoH{$sl[4] - 1}->{"ref_nuc"}), "\t";				
		}elsif($sl[7] eq 'R'){
			print $coordmap_HoH{$sl[4] - 1}->{"ref_nuc"}, 'N', $coordmap_HoH{$sl[4] + 1}->{"ref_nuc"}, "\t";		
		}else{
			die ("\nERROR: Unexpected orientation $sl[7] in following line:\n$_\n\n");
		}
	}else{
		die ("\nERROR: Unexpected reference base $sl[8] in following line:\n$_\n\n");
	}
	
	print $coordmap_HoH{$sl[4]}->{"type"}, "\t";
	if ($coordmap_HoH{$sl[4]}->{"type"} eq 'intergenic'){
		print ".\t.\t.";
	}else{
		print $coordmap_HoH{$sl[4]}->{"gene"}, "\t";
		if ($coordmap_HoH{$sl[4]}->{"strand"} == 1){
			print "F\t";
		}elsif ($coordmap_HoH{$sl[4]}->{"strand"} == -1){
			print "R\t";
		}else{
			die ("\nERROR: Unexpected gene orientation in coord map for position $sl[4]\.\n\n");	
		}
	}
	if ($coordmap_HoH{$sl[4]}->{"type"} eq 'CDS'){
		my $ref_codon;
		my $gene = $coordmap_HoH{$sl[4]}->{"gene"};
		my $gene_pos = $coordmap_HoH{$sl[4]}->{"gene_pos"};
		my $gene_strand;
		if ($coordmap_HoH{$sl[4]}->{"strand"} == 1){
			$gene_strand = "F";
		}elsif ($coordmap_HoH{$sl[4]}->{"strand"} == -1){
			$gene_strand = "R";
		}else{
			die ("\nERROR: Unexpected gene orientation in coord map for position $sl[4]\.\n\n");	
		}
		
		my $codon_pos = $coordmap_HoH{$sl[4]}->{"gene_pos"} % 3;
		$codon_pos or $codon_pos = 3;
		if ($codon_pos == 1){
			$ref_codon = $gene_HoH{$gene}->{$gene_pos} . $gene_HoH{$gene}->{$gene_pos + 1} . $gene_HoH{$gene}->{$gene_pos + 2};
		}elsif ($codon_pos == 2){
			$ref_codon = $gene_HoH{$gene}->{$gene_pos - 1} . $gene_HoH{$gene}->{$gene_pos} . $gene_HoH{$gene}->{$gene_pos + 1};		
		}elsif ($codon_pos == 3){
			$ref_codon = $gene_HoH{$gene}->{$gene_pos - 2} . $gene_HoH{$gene}->{$gene_pos - 1} . $gene_HoH{$gene}->{$gene_pos};		
		}else{
			die ("\nERROR: Incorrect codon position parsing for the following line:\n:$_\n\n");	
		}

		my $mut_codon = $ref_codon;
		my $replacement_base;
		if ($gene_strand eq $sl[7]){
			$replacement_base = $sl[9];
			$sl[8] eq substr($ref_codon, $codon_pos-1, 1) or die ("\nERROR: disagreement between $variant_file and $coordmap_file for reference base at position $sl[4]\:\n$_\n\n");
		}else{
			$replacement_base = revcom ($sl[9]);
			revcom ($sl[8]) eq substr($ref_codon, $codon_pos-1, 1) or die ("\nERROR: disagreement between $variant_file and $coordmap_file for reference base at position $sl[4]\:\n$_\n\n");
		}
		substr ($mut_codon, $codon_pos-1, 1, $replacement_base);

		my $ref_aa = codon2aa ($ref_codon);
		my $mut_aa = codon2aa ($mut_codon);
		
		if ($ref_aa eq $mut_aa){
			print "Syn ($ref_codon\>$mut_codon)";
		}else{
			print "$ref_aa\>$mut_aa ($ref_codon\>$mut_codon)";
		}		
	}else{
		unless ($coordmap_HoH{$sl[4]}->{"type"} eq 'intergenic'){
			print ".";
		}
	}
	print "\n";
}


