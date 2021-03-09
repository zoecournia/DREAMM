#!/usr/bin/perl -w 
# Syntax:
# perl Filtro_tab2csv.pl <tab-delimited file> <EN|ES>
#
# ES implies that the input file uses comma ',' for decimal numbers, in which case the separator in the outcome file will be semicolon ';'
# 
# EN implies that the input file uses dot '.' for decimal numbers, in which case the separator in the outcome file will be comma ','

use strict; 

if ($#ARGV != 1) {
	die "Incorrect arguments:\n
	# Syntax:
	# perl Filtro_tab2csv.pl <tab-delimited file> <EN|ES>\n";
}
open (DATATXT,"$ARGV[0]") or die "Couldn't open file $ARGV[0]\n";

my @name = split(/\./,$ARGV[0]);

my $name_csv = $name[0].'.csv';

open (DATACSV,">$name_csv") or die "Couldn't open file $name_csv\n";
if ($ARGV[1] eq "ES") {
	while (<DATATXT>) {
		$_ =~ s/\.//g;
		$_ =~ s/\t+\n/\n/;
		$_ =~ s/\t+/;/g;
		$_ =~ s/,0*;/;/g;
		print DATACSV "$_";
	}
} elsif ($ARGV[1] eq "EN") {
		while (<DATATXT>) {
		$_ =~ s/,//g;
		$_ =~ s/\t+\n/\n/;
		$_ =~ s/\t/,/g;
		$_ =~ s/\.0*,/,/g;
		print DATACSV "$_";
	}
} else {
	die "Please specify a valid format type (ES or EN)\n
	# Syntax:
	# perl Filtro_tab2csv.pl <tab-delimited file> <EN|ES>\n";
}

close DATACSV;
close DATATXT;
