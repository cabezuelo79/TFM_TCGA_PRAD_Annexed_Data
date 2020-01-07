#!usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

my @genes;
my $gene_file;
my $dirname;
my $help;
my $file;


GetOptions(
	"dirname|d=s"=>\$dirname,
	"genenames|g=s"=>\$gene_file,
	"help"=>\$help
);

if($dirname and $gene_file){

	my $filecount=0;
	my %hash;

	open(GENES,$gene_file);
		while (<GENES>) {
		chomp;
			my @data=split(/\t/);   
			my $symbol=$data[0];	
			push(@genes,$symbol);
	}
	close (GENES);

	opendir (DIR, $dirname) || die "Error in opening dir $dirname\n";
	my @files = readdir(DIR);
	foreach my $file (@files) {
		next if ($file eq "." or $file eq "..");
	   	open (FILE, "$dirname/$file") or die ("Can't open file\n");	
		$filecount++;		
		while (my $linea=<FILE>) {
			chomp;
			foreach my $gene (@genes){
				next if ($gene =~ /SYMBOL/);
				if($linea =~ $gene){
					$hash{$gene}{$file}++;
					next;
				}
			}
		}
		close FILE;
	}
	closedir(DIR);
	
print "GENE\tPERCENTAGE\tSAMPLE(PATIENT)\n";
	
foreach my $gene (keys %hash){
		print "$gene\t";
		print (scalar((keys %{$hash{$gene}}) / $filecount) * 100);
		print "%\t";
		print join(",", keys %{$hash{$gene}}) . "\n";
	}
}
else{
	print STDERR "USAGE: perl $0 -dirname=vfcs/ -genename=GENERAL.txt\n\n";
}


__END__



