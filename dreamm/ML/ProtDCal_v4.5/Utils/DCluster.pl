#!/usr/bin/perl
##################################################################################################################### 
############ Code developed by Ruiz-Blanco,Y.B., CAMD-BIR Unit. UCLV ################################################
############################## email: yasserblanco@sce.carleton.ca ##################################################
################## Core subroutines were taken from an external source. #############################################
##### The author assumes no responsibility because of the results of using this code ################################
#####################################################################################################################
#
# This version of the code is not practical for datasets with more that 500 descriptors
# The number of instances is not the most rate-limiting factor, because the majority of the processes scaling with the number of instances are 
# properly parallelized
#
# Perl v5.20.1
# This code implements an Spearman-correlation-based filter that removes redundant variables from an input CSV file.
#
# This is an unsupervised filter, the class column MUST be removed from the input file to prevent execution errors.
#
# An initial column with the labels of the instances may be present (this is assumed by default). 
#
# Please remove all thousands separators in your data.
#
# Please check that dots '.' are used for marking decimal numbers 
#
# The correlation coefficient (cc) threshold for initial filtering and clustering are fixable, as well as other parameters described below.
#  
# The code implements a single-linkage clustering algorithm using the Spearman cc as the metric for clustering. 
# The sum of the cc for every member of a cluster is used to identify the closest element to the centroid,
# which is thus selected as representative of all the attributes in the cluster in the finally reduced data set.
#
# Syntax: 
# perl Clustering_CSV_Spearman_Correlation_v1.6.pl -h
#
# Minimum syntax to access default parameters values:
# perl Clustering_CSV_Spearman_Correlation_v1.6.pl -d <input data>.csv
#
# Optional arguments:
# -t <threshold value of corr. coef. for clustering> # Default: 0.95 
#
# -d <data name>
#
# -n <0 or 1> # Default: '1' (true). It tells to the code about the presence of a column with the name of the instances at the 
# beginning of the data file. 1 for true and 0 for false. NOTE: if default value '1' is used without the presence of  
# the names column, the first-column's values are treated as labels of each instance.
#
# -c <0 or 1> # Default: '1' (true). It determines weather to save a file with the centroid of every formed cluster.
#
# -m <-1, 0, 1 or matrix_file_name> # Default: '1' (print true). It determines weather to save a file with the superior half of the 
# correlation matrix without the main diagonal elements. Alternatively to values '0' or '1' (false or true). Besides the value '-1' 
# instructs the code to die after saving the matrix file. Furthermore, it allows receiving directly the name of a file containing 
# a correlation matrix, in which case the matrix is directly read before clustering (NOte: This option allows avoiding the highly
# time-consuming step of building the corr. matrix if different threshold values will be assessed in separate calculations).  
# 
# -k <0 or 1> # Default '1' (true). It determines weather to check the input data.
#
# -x <positive integer> # Default: '2'. It determines the number of threads to use for computing the Corr. Mat. We recommend 
# a number larger than the amount of available cores for the calculation.
#
# -r <file name> # Default: 'null'. Should contain the name of a previously computed rank file

use strict;
use warnings;
use threads;
use threads::shared;

use List::Util qw(sum);

PrintReadme();
my $help = grep (/-h/,@ARGV);
if ($help > 0) {
	print "	Minimum syntax to access default parameters values:\n
	perl Clustering_CSV_Spearman_Correlation_v1.6.pl -d <input data>.csv
	\n";
	die "	Syntax and intructions were written in the file README.txt\n";
}

print "#### Setling Parameters: ####\n";
my %opt = @ARGV;
my $t = ((my $findt = grep (/-t/,@ARGV)) > 0)? $opt{'-t'}:0.95;
my $n = ((my $findn = grep (/-n/,@ARGV)) > 0)? $opt{'-n'}:1;
my $m = ((my $findm = grep (/-m/,@ARGV)) > 0)? $opt{'-m'}:1;
my $c = ((my $findc = grep (/-c/,@ARGV)) > 0)? $opt{'-c'}:1;
my $k = ((my $findch = grep (/-k/,@ARGV)) > 0)? $opt{'-k'}:1;
my $x = ((my $findx = grep (/-x/,@ARGV)) > 0)? $opt{'-x'}:2;
my $r = ((my $findr = grep (/-r/,@ARGV)) > 0)? $opt{'-r'}:"null";
my $data;
my $findd = grep (/-d/,@ARGV);
if ($findd > 0){
	$data = $opt{'-d'};
} else {
	print "Key in the name of your input file: ";
	$data = <STDIN>;
	chomp $data;
}

if (scalar(@ARGV) == 0){
print "Do you intend to change default settings: '-t 0.95' '-n 1' -m '1' -c '1' -k '1' -x '2' ? (y/n):";
my $answer = <STDIN>;
my $ok = index($answer,"y");
	if ($ok == 0){
		print "Threshold value (default: -t 0.95): ";
		$t = <STDIN>;
		print "There is a names column (default: -n 1 [true]): ";
		$n = <STDIN>;
		print "Print Corr. Matrix (default: -m 1 [print true]): ";
		$m = <STDIN>;
		print "Print Centroids List (default: -c 1 [true]): ";
		$c = <STDIN>;
		print "Check input file (default: -k 1 [true]): ";
		$k = <STDIN>;
		print "Number of threads to use for computing Corr. Mat. (default: -x 2): ";
		$x = <STDIN>;
		$x = int($x);
	}
}

##### Warnings ######
if ($n != 0 && $n != 1){
die "# Wrong value in option: -n. It must be 0 or 1 \n";
}

if ($c != 0 && $c != 1){
die "# Wrong value in option: -c. It must be 0 or 1 \n";
}

if ($k != 0 && $k != 1){
die "# Wrong value in option: -k. It must be 0 or 1 \n";
}

if ($m =~ /^-{0,1}\d$/ && $m != 0 && $m != 1 && $m != -1){
die "# Wrong numeric value in option: -m. It must be -1, 0 or 1 \n";
}

if ($t > 1 || $t < 0){
die "# Wrong value in option: -t. It must be a value within the range [0:1) \n";
}

if ($x < 1){
die "# Wrong value in option: -x. It must be a positive integer \n";
}

#### END Warnings #####

if ($k == 1){ 
	print "#### Checking input file #####\n";
	DoDataChecking ($data,$n);
}
open (ALL,"$data");
my @all = readline ALL;
close ALL;

my $data_changed; 
if ($r eq "null" && $m =~ /^-{0,1}\d$/) {
	$data_changed = PrintRanks (\@all,$data);
} else {
	$data_changed = $r;
}

##### Creating Corr. Matrix Array
my @matrix_array = ();
if ($m !~ /^-{0,1}\d$/) {
	@matrix_array = ReadCorrMatrix_mod($m);
} else {
	@matrix_array = ThreadedComputeCorrMat2($data_changed,$m);
	if ($m == -1){
	  die "### Execution ended after saving the corr. matrix\n";
	}
}

my @ClustersList = DoClustering (\@all,\@matrix_array);

my @centroids = BuildingCentroidList($data,\@ClustersList,\@matrix_array);

ReducingData(\@centroids,\@all,$data,$n);


############### Subroutine definitions ##################

################## From external source #################

sub convert_values_to_ranks {
	my $values = shift;

	# code below is slightly unintuitive, but we have to compute
	# average rank in cases where there are ties.
	my $idx = 0;
	my %sorted; # $sorted{$val} = [ $rank1, $rank2, $rank3 ];
	foreach my $val ( sort { $a <=> $b } @$values ) {
		push @{$sorted{$val}}, $idx;
		$idx++
		}

	# compute the average rank for a given value
	my %average_ranks =
	map { $_ => List::Util::sum(@{$sorted{$_}}) / scalar(@{$sorted{$_}}) }
	keys %sorted;

	# encode the values using average rank of the value in the list
	my @ranks = map { $average_ranks{$_} } @$values;

	if ( wantarray ) { return @ranks; }
	return \@ranks;
}

sub calculate_spearman_correlation_mod2 {
	my $ranked_n1 = shift;
	my $ranked_n2 = shift;

	if ( scalar(@{$ranked_n1}) != scalar(@{$ranked_n2}) ) {
	die "Error: spearman correlation given two lists of unequal size!\n";
	}

	my $sum_diff_squared = 0;
	foreach my $idx ( 0 .. scalar(@$ranked_n1)-1 ) {
		my $diff = $ranked_n1->[$idx] - $ranked_n2->[$idx];
		$sum_diff_squared += $diff * $diff;
		}

	my $N = scalar(@$ranked_n1);
	my $rho = 1 - ( 6 * $sum_diff_squared / ($N * ($N*$N-1)) );
	return $rho;
}

################## Own subroutines #################

sub ReducingData {
	print "######### Filtering Initial Dataset #######\n";
	my $selv = $_[0];
	
	# Getting variable's names:
	my @data_array = @{$_[1]};
	chomp @data_array;
	my $var_names = shift(@data_array);
	my @var_names_array = split(/['\t'|',']/,$var_names);

	# Opening output file for reduced data:
	my @name = split(/\./,$_[2]);
	my $name_csv = $name[0].'_red.csv';
	open (DATAREDCSV,">$name_csv");

	# Creating an array with the positions of the selected variables
	my $selected = join (',',@$selv);
	$selected .= ','; # See explanation in line 218
	my @pos_list = ();
	if ($_[3] == 1){
		push (@pos_list,0);
	}
	my $count = 0;
	foreach my $var (@var_names_array){
		my $var1 = $var.','; # this dues to preventing the case in which the name of a variable is contained in the name of a second one 
		my $match = index($selected,$var1);
		if ($match > -1){
			push (@pos_list,$count);
		}	
		$count++
	}

	# Copying selected item from each line to reduced data file
	my $num_selected = scalar(@pos_list);
	my @red_line = ();
	# printing first line with variable names
	for (my $i = 0; $i < $num_selected; $i++){
		push (@red_line,$var_names_array[$pos_list[$i]]);
	}
	my $new_line = join (',',@red_line);
	print DATAREDCSV "$new_line\n";

	# printing the rest of the lines
	foreach my $raw_line (@data_array){
		chomp $raw_line;
		my @line = split(/[\t|,|\s]/,$raw_line);
		@red_line = ();
		for (my $i = 0; $i < $num_selected;$i++){
			push (@red_line,$line[$pos_list[$i]]);
		}
		$new_line = join (',',@red_line);
		print DATAREDCSV "$new_line\n";
	}
	close DATAREDCSV;
}

sub BuildingCentroidList {
	# Opening list of centroid variables for writing if requested:
	my @name = split(/\./,$_[0]);
	my $cent = $name[0].'_list.csv';
	if ($c == 1){
		open (VARCSV,">$cent") || die "# Couldn't create $cent\n";
		print VARCSV "Centroid attributes,Size of Cluster\n";
		close VARCSV;
	}

	print "##### Building list of the centriods of each cluster ####\n";
	my $clusters = $_[1];
	my @matrix = @{$_[2]};
	my @centroids = ();
	if ($c == 1) {
		open (VARCSV,">>$cent") || die "# Couldn't find $cent\n";
	}
	
	for (my $cl = 0; $cl < ($#$clusters + 1); $cl++){
		my @clust = @{@$clusters[$cl]};
		my $csize = $#clust + 1;
		my $max_sum = 0;
		my $min_var = $clust[0];
		for(my $i = 0; $i < $csize; $i++){
			my $var1 = $clust[$i];
			my $sum = 0;
			for(my $j = 0; $j < $csize; $j++){
				my $var2 = $clust[$j];
				if ($var1 ne $var2){
					$sum += abs(GetCorr_mod($var1,$var2,\@matrix));
				}
			}
			if ($sum > $max_sum) {
			$max_sum = $sum;
			$min_var = $var1;
			} 
		}
		if ($c == 1){
			print VARCSV "$min_var,$csize\n";
		}
		push (@centroids,$min_var);
	}
	if ($c == 1){
		print VARCSV "#### END ####\n";
		close VARCSV;
	}
	return @centroids;
}

sub ReadCorrMatrix_mod {

	print "#### Reading Corr Matrix from file ####\n";
	open (CORRMAT,"$_[0]") || die "# Couldn't find $_[0]\n";
	my @corr_mat;
	while (<CORRMAT>){
	  my @temp = split (/,/,$_);
	  push (@corr_mat, @temp);
	}
	shift(@corr_mat);
	shift(@corr_mat);
	pop(@corr_mat);
	chomp(@corr_mat);
	return @corr_mat;
	close CORRMAT;
}

sub GetCorr_mod {
  my %mat = @{$_[2]};
  my $key1 = "$_[0]".":"."$_[1]";
  my $key2 = "$_[1]".":"."$_[0]";
  if ($mat{$key1}){
    return $mat{$key1};
  } else {
    return $mat{$key2};
  }
}

sub DoClustering {

	print "#### Building Clusters (This step may take some time) ####\n";

	# Getting variable's names:
	my @data_array = @{$_[0]};
	my @matrix = @{$_[1]};
	close DATATXT;
	chomp @data_array;
	my $var_names = shift(@data_array);
	my @copy_var_names_array = split(/["\t"|',']/,$var_names);

	my @clusters_list = ();
	if ($n == 1) {
		shift(@copy_var_names_array);
	}

	# Clustering
	for (my $i = 0; $i < scalar(@copy_var_names_array);){
		my $var1 = shift(@copy_var_names_array);
		my @cluster = ($var1);
		my $size = scalar(@cluster);
		for (my $k = 0; $k < $size; $k++) {
			my $var = $cluster[$k];
			my $size2 = scalar(@copy_var_names_array);
			for (my $j = 0; $j < $size2; $j++){
				my $var2 = $copy_var_names_array[$j];
				my $cc = abs(GetCorr_mod($var,$var2,\@matrix));
				if ($cc > $t){
					$var2 = splice(@copy_var_names_array,$j,1);
					$j--;
					push (@cluster,$var2);
					$size2 = scalar(@copy_var_names_array);
				}
			}
			$size = scalar(@cluster);
		}
		push(@clusters_list,[@cluster]);
	}
	return @clusters_list;
}

sub DoDataChecking {
	open (DATATXT,$_[0]) || die "# Couldn't open $_[0]\n";
	my @data_array = readline DATATXT;
	close DATATXT;
	chomp @data_array;
	
	# Deleting line-change character and field separators at the end of a line
	foreach my $file_line (@data_array) {
		$file_line =~ s/^(.+)[\t,\s]*$/$1/;
	}
	
	my $var_names = shift(@data_array);
	my @copy_var_names_array = split(/["\t"|',']/,$var_names);
	
	my @name = split('\.',$_[0]);
	open (LOG,">$name[0].log");
	
	my $total = scalar(@copy_var_names_array); 
	
	for (my $line = 0; $line < scalar(@data_array); $line++){
		my @line_array = split(/["\t"|',']/,$data_array[$line]);
		my $length = scalar(@line_array);
		if ($length != $total){
			print LOG "The size of line ".($line + 2)." differs from first line\n";
		}
		for(my $i = $_[1]; $i < $length; $i++){
			if (($line_array[$i] !~ /\d/) || ($line_array[$i] =~ /[a-df-zA-DF-Z]/)) {
				print LOG "Found problem in < $line_array[$i] > at: line ".($line + 2)." column ".($i + 1)."\n";
				close LOG;
				open (LOG,">>$name[0].log");
			}
		}
	}
	close LOG;
}

sub PrintRanks {
	print "####### Building Ranks #######\n";
	my @ranks_file_array : shared;
	# Getting variable's names:
	my @data_array = @{$_[0]};
	# Deleting line-change character and field separators at the end of a line
	foreach my $file_line (@data_array) {
		$file_line =~ s/^(.+)[\t,\s]*$/$1/;
	}
	
	my $var_names = shift(@data_array);
	my @var_names_array = split(/[,\t]/,$var_names);
	my $num_var = scalar(@var_names_array);
	my $v = 0;
	my @threads;	
	for (my $itr1 = $n; $itr1 < $num_var; $itr1++){
		
		# Obtaining the values of the first variable:
		
		
		if ($v < $x){
			push (@threads,threads->create(\&Ranks_core,\@var_names_array,$itr1,\@data_array,\@ranks_file_array));
			$v++;
		} else {
			foreach my $thr (@threads){
				$thr->join();
			}
			$v = 0;
			@threads = ();
			$itr1--;
		}
		
	}
	foreach my $thr (@threads){
		$thr->join();
	}
	
	my @name = split (/\./,$_[1]);
	my $new_data = $name[0]."_ranks.csv";
	open (RANKS, ">$new_data");
	foreach my $k (@ranks_file_array){
	  print RANKS "$k\n";
	}
	
	close RANKS;
	return $new_data;

}

sub Ranks_core {
	my @var1 = ();
	my $ranks_line = "";
	foreach my $file_line (@{$_[2]}) {
		my @line = split(/[\t,]/,$file_line); 
		push @var1,"$line[$_[1]]";
	}
	my @var_names_array1 = @{$_[0]};
	my $copyvar1 = convert_values_to_ranks(\@var1);
	$ranks_line .= "$var_names_array1[$_[1]]";
	foreach my $k (@$copyvar1){
	  $ranks_line .= ",$k";
	}
	
	{
	lock($_[3]);
	push (@{$_[3]},($ranks_line));
	}
}

sub ThreadedComputeCorrMat2 {
	
	# Opening raw data file for reading
	open (DATATXT,"$_[0]") || die "# Couldn't open $_[0]\n";
	my @mod_data : shared; 
	@mod_data = readline DATATXT;
	my $num_var = scalar(@mod_data);

	# defining output matrix array
	my @matrix_array :shared;
	my $m_array	:shared = \@matrix_array;
	my @name = split(/\./,$_[0]);
	my $corr_mat = $name[0].'_matrix.csv';
	my $Z;
	if ($_[1] != 0) {
		print "##### Building Correlation matrix #####\n";
		open (MATRIX,">$corr_mat") || die "# Couldn't create $corr_mat\n";
		print MATRIX "FirstAttribute:SecondAttribute,SpearmCorrelation\n";
		close MATRIX;
		open (MATRIX,">>$corr_mat") || die "# Couldn't find $corr_mat\n";
	} else {
		print "##### Building Correlation matrix #####\n";
	}
	
	my @threads = ();
	
	for (my $itr1 = 0; $itr1 < $num_var - 1;){
		
		
		#******* Exploring the array of variables forward using threads:
		# *********** threaded algorithms *********
		my $itr2 = $itr1 + 1;
		if ($itr1 > $x - 1) {
			while ($itr1 < $num_var - 1) {
				for (my $t = 0; $t < scalar(@threads);) {
					my $thr = $threads[$t];
					if ($thr->is_joinable() && $itr1 < $num_var - 1){
						$thr->join();
						splice (@threads,$t,1);
						my @var1 = split(/[\t,]/,$mod_data[$itr1]);
						
						my $name1 = shift(@var1);
						my $copyvar1 = \@var1;
						
						push (@threads,threads->create(\&PrintMat_mod3,$copyvar1,$itr2,$name1,$corr_mat,$m_array,\@mod_data,$num_var));
						$itr2++;
						$itr1++;
						
					} elsif ($itr1 >= $num_var - 1) {
						last;
					} else {
						$t++;
					}
				}
			}
		} else {
			# Obtaining the values of the first variable:
			my @var1 = split(/[\t,]/,$mod_data[$itr1]);
			
			my $name1 = shift(@var1);
			my $copyvar1 = \@var1; 
		
			push (@threads,threads->create(\&PrintMat_mod3,$copyvar1,$itr2,$name1,$corr_mat,$m_array,\@mod_data,$num_var));
			$itr1++;
		}
		
		# *********** END of threaded algorithms *********
	}
		
	foreach my $thr (@threads){
		$thr->join();
	}
	
	if ($_[1] != 0) {
		my $z = scalar(@matrix_array);
		
		for (my $t = 0; $t < $z - 1; $t += 2){
			my $o = $t + 1;
			print MATRIX "$matrix_array[$t],$matrix_array[$o]\n";
		}
		print MATRIX "#### END ####\n";
		close MATRIX;
	}
	return @matrix_array;
}

sub PrintMat_mod3 {
	my $i2 = $_[1];
	my @v1 = @{$_[0]};
	my @mod_data1 = @{$_[5]};
	my $n1 = $_[2];
	while ($i2 < $_[6]){
		my @var2 = split(/[\t,]/,$mod_data1[$i2]);
		
		my $name2 = shift(@var2);
		my $copyvar2 = \@var2;
		PrintMat_mod2(\@v1,$i2,$n1,$name2,$_[3],$_[4],$copyvar2);
		$i2++;
	}
}

sub PrintMat_mod2 {
	my $corr_coef = calculate_spearman_correlation_mod2($_[6],$_[0]);
	my $temp = "$_[2]".":"."$_[3]";
	# Saving Correlation Matrix 
	{
	lock($_[5]);
	push (@{$_[5]},($temp,$corr_coef));
	}
}

sub PrintReadme {
	open (HELP,">README.txt");
	
	print HELP 
"##################################################################################################################### 
############ Code developed by Ruiz-Blanco,Y.B., CAMD-BIR Unit. UCLV ################################################
######################## email: yasserblanco[at]sce.carleton.ca #####################################################
################## Core subroutines were taken from an external source. #############################################
##### The author assumes no responsibility because of the results of using this code ################################
#####################################################################################################################
#
# This version of the code is not practical for datasets with more that 500 descriptors
# The number of instances is not the most rate-limiting factor, because the majority of the processes scaling with the number of instances are 
# properly parallelized
#
# Perl v5.20.1
# This code implements an Spearman-correlation-based filter that removes redundant variables from an input CSV file.
#
# This is an unsupervised filter, the class column MUST be removed from the input file to prevent execution errors.
#
# An initial column with the labels of the instances may be present (this is assumed by default). 
#
# Please remove all thousands separator in your data.
#
# Please check that dots '.' are used for marking decimal numbers 
#
# The correlation coefficient (cc) threshold for initial filtering and clustering are fixable, as well as other parameters described below.
#  
# The code implements a single-linkage clustering algorithm using the Spearman cc as the metric for clustering. 
# The sum of the cc for every member of a cluster is used to identify the closest element to the centroid,
# which is thus selected as representative of all the attributes in the cluster in the finally reduced data set.
#
# Syntax: 
# perl Clustering_CSV_Spearman_Correlation_v1.6.pl
#
# Minimum syntax to access default parameters values:
# perl Clustering_CSV_Spearman_Correlation_v1.6.pl -d <input data>.csv
#
# Optional arguments:
# -t <threshold value of corr. coef. for clustering> # Default: 0.95 
#
# -d <data name>
#
# -n <0 or 1> # Default: '1' (true). It tells to the code about the presence of a column with the name of the instances at the 
# beginning of the data file. 1 for true and 0 for false. NOTE: if default value '1' is used without the presence of  
# the names column, the first-column's values are treated as labels of each instance.
#
# -c <0 or 1> # Default: '1' (true). It determines weather to save a file with the centroid of every formed cluster.
#
# -m <-1, 0, 1 or matrix_file_name> # Default: '1' (print true). It determines weather to save a file with the superior half of the 
# correlation matrix without the main diagonal elements. Alternatively to values '0' or '1' (false or true). Besides the value '-1' 
# instructs the code to die after saving the matrix file. Furthermore, it allows receiving directly the name of a file containing 
# a correlation matrix, in which case the matrix is directly read before clustering (NOte: This option allows avoiding the highly
# time-consuming step of building the corr. matrix if different threshold values will be assessed in separate calculations).  
# 
# -k <0 or 1> # Default '1' (true). It determines weather to check the input data.
#
# -x <positive integer> # Default: '2'. It determines the number of threads to use for computing the Corr. Mat. We recommend 
# a number larger than the amount of available cores for the calculation.
#
# -r <file name> # Default: 'null'. Should contain the name of a previously computed rank file\n";
close HELP;
}
