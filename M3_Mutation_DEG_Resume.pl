#!usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

open (MUTATED, 'file obtained from PercentageVCF_MutatedGene.pl');
my %hash;

while (<MUTATED>) {
	chomp;
		my @data=split(/\t/);   
		my $gen=$data[0];	
		my $percentage=$data[1];  
		my $samples=$data[2];	
		$hash{$gen}{$percentage}=$samples;
}

close MUTATED;

open (GENERAL, "TCGA_DEG_file");
my %hash2;

while (<GENERAL>) {
	chomp;
		my @data=split(/\t/);   
		my $gen=$data[0];	
		my $deg=$data[2];  
		$hash2{$gen}=$deg;
}

close GENERAL;

open (PAIRED, "TCGA_PAIRED_DEG_file");
my %hash3;

while (<PAIRED>) {
	chomp;
		my @data=split(/\t/);   
		my $gen=$data[0];	
		my $deg=$data[2]; 
		$hash3{$gen}=$deg;
}

close PAIRED;

print "GENE\tPERCENTAGE\tDEG in General (LogFC)\tDEG in Paired (LogFC)\tSAMPLE (PATIENT)\n";

foreach my $gen (sort keys %hash) {
	foreach my $percentage (keys %{$hash{$gen}}) {
print "$gen\t$percentage\t" . $hash2{$gen}. "\t" . $hash3{$gen} . "\t" . $hash{$gen}{$percentage} . "\n";
	}
}

