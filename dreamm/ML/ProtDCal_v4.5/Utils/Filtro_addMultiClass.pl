#!/usr/bin/perl -w 
# Syntax
# perl Filtro_addMultiClass.pl <dataset.csv> <comma-separated list of class identifiers>
#
# E.g.:
# perl Filtro_addMultiClass.pl dataset.csv antitumor,antibacterial,antifungus
#
# The code assumes that the first column in the csv file gather the names of every instance, and those names contain the identifiers of the classes (e.g. antitumor,antibacterial,antifungus)  

use strict; 

open (DATATXT,"$ARGV[0]");

my @name = split(/\./,$ARGV[0]);

my $name_csv = $name[0].'_mclass.csv';

open (DATACSV,">$name_csv");

my @classes = split (',',$ARGV[1]);

my $firstline = 0;
while (<DATATXT>) {
	# chomp $_;
	$_ =~ s/,*\n//;
	if ($firstline == 0){
		print DATACSV "$_".",class\n";
		$firstline = 1;
	} else {
		foreach my $c (@classes){
			if ($_ =~ /$c/) {
				print DATACSV "$_".",$c\n";
				last;
			} 
		}
	}
}

close DATACSV;
close DATATXT;
