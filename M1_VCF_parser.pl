#!usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;   

my $file;
my $out;    

GetOptions(
	"file|f=s"=>\$file,
	"outfile|o=s"=>\$out   
	);

my $genes; 
my $info;
my $chr;
my $rs;
my $pos;   
my $ref;
my $alt;
my %hash;

if($file){
	open(VCF,'<', $file);	

	my %hash;
	my $id_paciente;
	while (<VCF>) {
		chomp;
		if($_ =~/##SAMPLE=<ID=TUMOR,NAME=(TCGA-..-....)/) {
			$id_paciente=$1;
		}
		if($_ =~/^chr/){		
			my @data=split(/\t/);   
			$info=$data[7];	
			$chr=$data[0];   
			$pos=$data[1];	
			$rs=$data[2];	
			$ref=$data[3];	
			$alt=$data[4];	
		}
		if ($info) {
			my @data2=split(/,/,$info);                                            
				foreach my $infos (@data2){
					my @data3=split(/\|/,$infos);	                                 
					$hash{"$chr\t$pos\t$rs\t$ref\t$alt"}{$data3[1]}{$data3[2]}=$data3[3] if($data3[1]); 
			}
		}
	}
	close VCF;
	
	open(OUT, '>', $out); 
	
	print OUT "PATIENT ID: $id_paciente\n";  
	print OUT "CHROMOSOME\tPOSITION\tID\tREF\tALT\tCONSECUENCE\tIMPACT\tSYMBOL\n";  
	
	foreach my $id (sort keys %hash) {
		foreach my $consecuence (sort keys %{$hash{$id}}) {
			foreach my $impact (keys %{$hash{$id}{$consecuence}}) {
				print OUT "$id\t$consecuence\t$impact\t" . $hash{$id}{$consecuence}{$impact} ."\n";
			}
		}
	}
	close OUT;	
	
}	
else{
	print STDERR "Error:\nUsage: perl $0 -file file.vcf -o out.tsv\n";
}

