#!/usr/bin/perl
#####################################################################################################################  
# 
# Minimum syntax to access default parameters values (using FASTA files and projects files of type '*.proj'):
#
# perl ProtDCal-parallel.pl -type <proj|idl|ppi|ppi_idl> -file <fasta|pdb>
# 
# The code assumes that Java is settled within the enviroment variables
#
# Optional arguments:
# -file <fasta|pdb> # Default 'fasta' It determines the type of input data (currently the code only accepts 'fasta')
#
# -ppi <path to directory containing the *.ppi file> # Default: 'Datasets/' 
#
# -v <0|1> # Default: 0 (not to include) It determines whether to include or not the name of teh project file in the features tag 
#
# -threads <number of threads to launch> # Default: '4' This parameter determines the number of used threads to distribute the run, as well as
#  the number of packages the data will be divided into.
#
# -in <Path to input directory> # Default: 'Datasets/Fasta_Protein_Format'
#
# -out <Path to output directory> # Default: 'Outputs'  
# 
# -proj <Path to projects directory> # Default: 'Projects' 
#
# -idl <Path to IDL-projects directory> # Default: 'Projects'
#
# -jar <Path to .jar file of ProtDcal> # Default: 'ProtDCal.jar'
#
# -Xmx <maximum heap space for Java VM> # Default: '1000m'  Please not that MUST be an space between the option and its value e.g. '-Xmx 2000m'
#
# -type <proj|idl|ppi|ppi_idl> # It determines the type of calculation
#
# -gz <true|false> # Default: "true" It determines whether to zip the results file.
#

use strict;
use Cwd;
use warnings;
use threads;
use threads::shared;
use File::Spec::Functions;
use File::Path;
use File::Copy;
use File::Find;
use File::Basename;
no warnings 'File::Find';
use List::Util qw(sum);
use IO::Compress::Gzip qw(gzip $GzipError);

PrintReadme();

my $help = grep (/-h/,@ARGV);
if ($help > 0) {
	print "	Minimum syntax to access default parameters values:\n
	perl ProtDCal-parallel.pl -type <idl|proj|ppi|ppi_idl> -file <fasta|pdb>
	\n";
	die "	Syntax and intructions were written in the file README-parallel.txt\n";
}

print "#### Setting Parameters: ####\n";
my %opt = @ARGV;
my $t = ((my $findt = grep (/-ppi/,@ARGV)) > 0)? $opt{'-ppi'}:"Datasets";
my $v = ((my $findv = grep (/-v/,@ARGV)) > 0)? $opt{'-v'}:0;
my $type = ((my $findtype = grep (/-type/,@ARGV)) > 0)? $opt{'-type'}:"null";
my $f = ((my $findf = grep (/-file/,@ARGV)) > 0)? $opt{'-file'}:"null";
if ($f eq "null") {
	print "Please enter a valid data type for the calculation <fasta|pdb>: ";
	$f = <STDIN>;
	chomp $f;
}
my $default_i;
if ($f eq "pdb") {
	$default_i = catdir('Datasets','Pdb_Protein_Format');
} elsif ($f eq "fasta") {
	$default_i = catdir('Datasets','Fasta_Protein_Format');
} else {
	die "Please check that -f must be either 'fasta' or 'pdb'\n";
}
my $i = ((my $findi = grep (/-in/,@ARGV)) > 0)? $opt{'-in'}:$default_i;
my $o = ((my $findo = grep (/-out/,@ARGV)) > 0)? $opt{'-out'}:"Outputs";
my $p = ((my $findp = grep (/-proj/,@ARGV)) > 0)? $opt{'-proj'}:"Projects";
my $x = ((my $findx = grep (/-idl/,@ARGV)) > 0)? $opt{'-idl'}:"Projects";
my $pdc = "ProtDCal.jar";
my $jar = ((my $findjar = grep (/-jar/,@ARGV)) > 0)? $opt{'-jar'}:$pdc;
my $n = ((my $findn = grep (/-threads/,@ARGV)) > 0)? $opt{'-threads'}:4;
my $xmx = ((my $findxmx = grep (/-Xmx/,@ARGV)) > 0)? $opt{'-Xmx'}:"1000m";
my $gz = ((my $findgz = grep (/-gz/,@ARGV)) > 0)? $opt{'-gz'}:"true";
my $ppi_name = "pairs";

#### Warnings #####

if (($type ne "ppi") && ($type ne "idl") && ($type ne "ppi_idl") && ($type ne "proj")) {
	print "Please enter a valid type for the calculation <proj|idl|ppi|ppi_idl>: ";
	$type = <STDIN>;
	chomp $type;
}

if (($f ne "pdb") && ($f ne "fasta")) {
	print "Please enter a valid data type for the calculation <fasta|pdb>: ";
	$f = <STDIN>;
	chomp $f;
}

if ($f ne "fasta" && ($type eq "ppi" || $type eq "ppi_idl")) {
	die "The calculations type $type are only compatible with data type 'fasta' (-f fasta). Please check your input parameters\n";
}

print "Parameters Summary:
-ppi = $t (Only used for -type <ppi|ppi_idl>)
-v = $v
-proj = $p (Only used for -type <ppi|proj>)
-out = $o
-in = $i (Only used for -type <idl|ppi_idl>)
-idl = $x (Only used for -type <idl|ppi_idl>)
-file = $f (Only used for -type <idl|ppi_idl>)
-jar = $jar
-threads = $n
-Xmx = $xmx
-type = $type
-gz = $gz\n";

###### END Warnings #####

###### Splitting Dataset #####
print "#### Distributining Instances ####\n";
	
if ($f eq "fasta" && ($type eq "proj" || $type eq "idl")) {
	opendir (my $dir, $i);
	my @fasta_files = grep (/\.fasta/,readdir($dir));
	unless (scalar(@fasta_files) >= 1) {
	  	die "No file with 'fasta' extension was found in $i\n";
	}
	my $path = catdir($i,'SPLIT');
	mkdir $path,0775;
	foreach my $file (@fasta_files){
		SplitFasta ($i,$file,$n);
	}
} elsif ($f eq "pdb") {
	opendir (my $dir, $i);
	my @pdb_files = grep (/\.pdb/,readdir($dir));
	push (@pdb_files,grep (/\.ent/,readdir($dir)));
	unless (scalar(@pdb_files) >= 1) {
		die "Please check that $i is the valid directory for the input dataset. 
Otherwise declare a valid path using the '-i' option.
If '-type proj' is being used, the same directory path must appear in the projects files (*.proj)\n";
	}
	my $currdir = getcwd();
	chdir $i;
	mkdir 'SPLIT',0775;
	chdir "SPLIT";
	for (my $k = 0; $k < $n; $k++) {
		mkdir $k,0775; 
	}
	
	chdir (updir());
	my $block = int (scalar(@pdb_files)/$n);
	my $c = 0;
	for (my $k = 0; $k < scalar(@pdb_files); $k++) {
		if ($k < $block * ($c + 1)) {
			my $destination = catfile ('SPLIT',$c,$pdb_files[$k]);
			copy($pdb_files[$k],$destination); 
		} elsif ($k >= $n * int (scalar(@pdb_files / $n))) {
			my $destination2 = catfile ('SPLIT',$c,$pdb_files[$k]);
			copy($pdb_files[$k],$destination2); 
		} else {
			$c++;
			$k--;
		}
	}
	chdir $currdir;
} 


##### END Splitting Dataset ####

#### Splitting project files ####
if ($type eq "proj") {
 	my $currdir = getcwd();
 	opendir (my $pdir, $p);
	my @proj = grep (/\.proj/,readdir($pdir));
	chdir $p;
	for (my $d = 0; $d < $n; $d++) {
		mkdir $d,0775;
		foreach my $k (@proj){
			open (PROJ,"$k");
      
			$k =~ /(.+)\.proj/;
			my $new_p = $1."_".$d.".proj";
			chdir $d;      
			open (NEW,">$new_p");
			while (<PROJ>){
				if ($. == 2){
					my $datadir = catdir($i,'SPLIT',$d);
					$_ =~ s/^.+$/$datadir/;
				}
				print NEW "$_";
			}
			close NEW;
      
			chdir (updir());
			close PROJ;
		}
	}
	chdir $currdir;
}

### Splitting .ppi file ####

if ($type eq "ppi" || $type eq "ppi_idl") {
	my $currdir = getcwd();
	opendir (my $tdir, $t);
	my @ppi = grep (/\.ppi/,readdir($tdir));
	if (scalar(@ppi) > 1) {
		die "More than 1 *.ppi files were found. Please make sure than the input directory contains only 1 *.ppi file\n";
	}	
	my $start = 0;
	my $end;
	my $path = catdir($t,'SPLIT');
	mkdir $path,0775;
	chdir $path;
	
	for (my $d = 0; $d < $n; $d++) {
		mkdir $d,0775;
		
		foreach my $k (@ppi){
			chdir (updir());
			open (PPI,"$k");
			chdir "SPLIT"; 
			$k =~ /(.+)\.ppi/;
			$ppi_name = $1;
			my $new_p = $1."_".$d.".ppi";
			chdir $d;
			my @ppi_lines = readline PPI;
			my $lines_num = scalar (@ppi_lines);
			if (($n - $d) > 1) {
				$end = int ($lines_num / $n) * ($d + 1); 
			} else {
				$end = $lines_num;
			}
			$start = int ($lines_num / $n) * $d;     
			open (NEW,">$new_p");
			for (my $c = 0; $c <= scalar(@ppi_lines); $c++){
				if ($c >= $start && $c < $end){
					print NEW "$ppi_lines[$c]";
				} elsif ($c == $end){
					last;
				}
			}
			close NEW;
			chdir (updir());
			close PPI;
		}
	}
	chdir $currdir;
}

#### Launching threads ####

print "#### Initializing PDC jobs ####\n";
my $count = 0;
my $fname = "PDC.".$count.".out";
if (-e "PDC.0.out"){
	my $currdir = getcwd();
	opendir (my $outdir, $currdir);
	my @outs = grep (/\.out/,readdir($outdir));
	  
	while (-e $fname){
		$count++;
		$fname = "PDC.".$count.".out";
	}
}
my @threads = ();
for (my $d = 0; $d < $n; $d++) {
	my $inputdir = catdir($i,'SPLIT',$d);
	push (@threads,threads->create(\&PDCjob,$inputdir,$d,$fname,$ppi_name));
}
 
foreach my $thr (@threads){
	$thr->join();
}

print "#### Processing split outcomes ####\n";

#Deleting split input files

if ($type ne "ppi" && $type ne "ppi_idl") {
	rmtree(catdir($i,'SPLIT')); 	 
} else {
	my $ppi_dir = catdir($t,'SPLIT');	
	if (-d $ppi_dir) {
		rmtree($ppi_dir); 	
	}
}

#Moving split output files to directory $o
my @outfiles;
for (my $d = 0; $d < $n; $d++) {
  if ($type ne "proj") {
  	my $outcome = catdir($o,$d);
	@outfiles = ();
	find(\&wanted,$outcome);
	foreach my $file (@outfiles) {
		my @parse = fileparse($file,qr/\.[^.]*/);
		my $name = "$parse[0]"."_$d"."$parse[2]";
		my $filedest = catfile($o,$name);
		move($file,$filedest);
	}
	rmtree($outcome);#Deleting empty directories
  } else {
	  my $projdir = catdir($p,$d); 
	  rmtree($projdir);  # Removing modified projects

	  my $outcome = catdir($o,$d);
	  @outfiles = ();
	  find(\&wanted,$outcome);
	  foreach my $file (@outfiles) {
		my @parse = fileparse($file,qr/\.[^.]*/);
		my $name = "$parse[0]"."$parse[2]";
		my $filedest = catfile($o,$name);
		move($file,$filedest);
	  }
	  rmtree($outcome);#Deleting empty directories
  }
}

#Combining output files and removing parts
if ($type eq "proj") {
	my $currdir = getcwd();
  # Reading projects' names
	opendir (my $pdir, $p);
	my @proj = grep (/\.proj/,readdir($pdir));
	chdir $o;
  # Combining results files
	foreach my $k (@proj){
    # Creating an array with all the output files names
		my @outnames;
		$k =~ /(.+)\.proj/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_p = $nam."_".$d."_Prot.txt";
		push (@outnames,$new_p); 
		}
	
	# Copying the content of all the partial results files into the first one
		foreach my $m (@outnames){
			open (PROJ,$m);
			my $name0 = $nam."_"."0_Prot.txt";
			open (PROJ0,">>$name0");
			my @content = readline (PROJ);
			close PROJ;
			shift (@content);
			foreach my $l (@content){
				print PROJ0 $l;
			}
			close PROJ0;
			unlink($m);
		}  
	}
  # Combining log files
	foreach my $k (@proj){
    # Creating an array with all the log files names
		my @outnames;
		$k =~ /(.+)\.proj/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_p = $nam."_".$d.".log";
			push (@outnames,$new_p); 
		}
	# Copying the content of all the log files into the first one
		foreach my $m (@outnames){
			open (PROJ,$m);
			my $name0 = $nam."_"."0.log";
			open (PROJ0,">>$name0");
			my @content = readline (PROJ);
			close PROJ;
			foreach my $l (@content){
				print PROJ0 $l;
			}
 			close PROJ0;
			unlink($m);
		}  
	}
	chdir $currdir;
} elsif ($type eq "idl") {
	my $currdir = getcwd();
	# Reading idl-projects' names
  	opendir (my $xdir, $x);
  	my @idl = grep (/\.idl/,readdir($xdir));
  	chdir $o;
	
	# Combining results files
	foreach my $k (@idl){
		# Creating an array with all the output files names
		my @outnames;
		$k =~ /(.+)\.idl/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_x = $nam."_Prot"."_$d.txt";
			push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial results files into the first one
		foreach my $m (@outnames){
			open (IDL,$m);
			my $name0 = $nam."_Prot_0.txt";
			open (IDL0,">>$name0");
			my @content = readline (IDL);
			close IDL;
			shift (@content);
			foreach my $l (@content){
				print IDL0 $l;
			}
			close IDL0;
			unlink($m);
		}  
	}
	# Combining log files
	foreach my $k (@idl){
		# Creating an array with all the log files names
		my @outnames;
		$k =~ /(.+)\.idl/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_x = $nam."_$d".".log";
			push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial log files into the first one
		foreach my $m (@outnames){
			open (IDL,$m);
			my $name0 = $nam."_0.log";
			open (IDL0,">>$name0");
			my @content = readline (IDL);
			close IDL;
			foreach my $l (@content){
				print IDL0 $l;
			}
			close IDL0;
			unlink($m);
		}  
	}
	chdir $currdir;
} elsif ($type eq "ppi_idl") {
	my $currdir = getcwd();
	# Reading idl-projects' names
  	opendir (my $xdir, $x);
  	my @idl = grep (/\.idl/,readdir($xdir));
  	chdir $o;
	
	# Combining results files
	foreach my $k (@idl){
		# Creating an array with all the output files names
		my @outnames;
		$k =~ /(.+)\.idl/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_x = $nam."_PPI"."_$d.txt";
			push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial results files into the first one
		foreach my $m (@outnames){
			open (IDL,$m);
			my $name0 = $nam."_PPI_0.txt";
			open (IDL0,">>$name0");
			my @content = readline (IDL);
			close IDL;
			shift (@content);
			foreach my $l (@content){
				print IDL0 $l;
			}
			close IDL0;
			unlink($m);
		}  
	}
	# Combining log files
	foreach my $k (@idl){
		# Creating an array with all the log files names
		my @outnames;
		$k =~ /(.+)\.idl/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
			my $new_x = $nam."_$d".".log";
			push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial log files into the first one
		foreach my $m (@outnames){
			open (IDL,$m);
			my $name0 = $nam."_0.log";
			open (IDL0,">>$name0");
			my @content = readline (IDL);
			close IDL;
			foreach my $l (@content){
				print IDL0 $l;
			}
			close IDL0;
			unlink($m);
		}  
	}
	chdir $currdir;
} else {
	my $currdir = getcwd();
	# Reading proj-projects' names
  	opendir (my $pdir, $p);
  	my @proj = grep (/\.proj/,readdir($pdir));
  	chdir $o;
	
	# Combining results files
	foreach my $k (@proj){
		# Creating an array with all the output files names
		my @outnames;
		$k =~ /(.+)\.proj/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
		  my $new_x = $nam."_PPI"."_$d.txt";
		  push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial results files into the first one
		foreach my $m (@outnames){
		  open (PROJ,$m);
		  my $name0 = $nam."_PPI_0.txt";
		  open (PROJ0,">>$name0");
		  my @content = readline (PROJ);
		  close PROJ;
		  shift (@content);
		  foreach my $l (@content){
			print PROJ0 $l;
		  }
		  close PROJ0;
		  unlink($m);
		}  
	}
	# Combining log files
	foreach my $k (@proj){
		# Creating an array with all the log files names
		my @outnames;
		$k =~ /(.+)\.proj/;
		my $nam = $1;
		for (my $d = 1; $d < $n; $d++) {
		  my $new_x = $nam."_$d".".log";
		  push (@outnames,$new_x); 
		} 

		# Copying the content of all the partial log files into the first one
		foreach my $m (@outnames){
		  open (PROJ,$m);
		  my $name0 = $nam."_0.log";
		  open (PROJ0,">>$name0");
		  my @content = readline (PROJ);
		  close PROJ;
		  foreach my $l (@content){
			print PROJ0 $l;
		  }
		  close PROJ0;
		  unlink($m);
		}  
	}
	chdir $currdir;
}

# Gzipping output files

if ($gz eq "true") {
  print "#### Compacting output files (the original files are deleted) ####\n";
  my $currdir = getcwd();
  opendir (my $odir, $o);
  my @txt = grep (/\.txt/,readdir($odir));
  chdir $o;
  foreach my $z (@txt) {
    my @zname = split (/\./,$z);
	gzip $z => "$zname[0].gz" or die "gzip failed: $GzipError";
    unlink ($z);
  }
  chdir $currdir;
}

# Deleting empty log files
my $currdir = getcwd();
opendir (my $logdir, $currdir);
my @logs = grep (/\.log/,readdir($logdir));
foreach my $l (@logs) {
	unless (-s $l) {
		unlink ($l);
	}
}

opendir (my $ldir, $o);
my @ls = grep (/\.log/,readdir($ldir));

chdir $o;
foreach my $l (@ls) {
	unless (-s $l) {
		unlink ($l);
	}
}
chdir $currdir;

############### Subroutine definitions ##################

sub wanted {
  if (-f $_ ) {
    push (@outfiles,$File::Find::name);
  }
}

sub SplitFasta {
	my $fasta = catfile($_[0],$_[1]);
	open (FASTA,"$fasta") || die "# Couldn't find file $_[1] in $fasta\n";
	my @names_array;
	my $currentdir = getcwd();
	my $go = catdir($_[0],'SPLIT');
	chdir $go;
	for (my $j = 0; $j < $_[2]; $j++){
	  $_[1] =~ /(.+)\.fasta/;
	  my $name = $1."_$j.fasta";
	  if (-e $j){
	    chdir $j;
	  } else {
	    mkdir $j,0775;
	    chdir $j;
	  } 
	  
	  my $new = catfile ($j,$name);
	  open (NEW,">$name") || die "# Couldn't create file $name in $j\n";
	  push (@names_array,$new);
	  close NEW;
	  chdir (updir());
	}
	my $j = -1;
	while (<FASTA>){
	  
	if ($_ =~ /^>/) {
	    if ($j < ($_[2] - 1)) {
	      $j++;
	      open (my $mi, ">>$names_array[$j]");
	      print $mi "$_";
	      close $mi;
	    } else {
	      $j = 0;
	      open (my $mi, ">>$names_array[$j]");
	      print $mi "$_";
	      close $mi;
	    }
	  } elsif ($_ !~ /^$/) {
	    open (my $mi, ">>$names_array[$j]");
	    print $mi "$_";
	    close $mi;
	  }
	}
	chdir $currentdir;
}

sub PDCjob {
	my $output_dir = catdir($o,$_[1]);
	mkdir $output_dir,0775;
	my $proj_dir;

	if ($type eq "ppi" || $type eq "ppi_idl") { 
  		my $new_ppi = $_[3]."_".$_[1].".ppi";		
  		my $ppi_file = catfile($t,'SPLIT',$_[1],$new_ppi);
    
  		if ($type eq "ppi") {
  			system ("java -Xmx$xmx -jar $jar -o $output_dir -p $p -ppi $ppi_file >> $_[2]");
		} elsif ($type eq "ppi_idl") {
			system ("java -Xmx$xmx -jar $jar -f $f -i $i -o $output_dir -x $x -ppi $ppi_file >> $_[2]");
		}
	} elsif ($type eq "proj"){
		$proj_dir = catdir($p,$_[1]);
		system ("java -Xmx$xmx -jar $jar -o $output_dir -p $proj_dir -v $v >> $_[2]");
	} elsif ($type eq "idl") {
		system ("java -Xmx$xmx -jar $jar -f $f -i $_[0] -o $output_dir -x $x >> $_[2]");
	} else {
		die "Please enter a valid job type: 'proj', 'idl', 'ppi' or 'ppi_idl'\n";
	}
}


sub PrintReadme {
	open (HELP,">README-parallel.txt");
	
	print HELP 
"#####################################################################################################################  
# 
# Minimum syntax to access default parameters values (using FASTA files and projects files of type '*.proj'):
#
# perl ProtDCal-parallel.pl -type <proj|idl|ppi|ppi_idl> -file <fasta|pdb>
# 
# The code assumes that Java is settled within the enviroment variables
#
# Optional arguments:
# -file <fasta|pdb> # Default 'fasta' It determines the type of input data (currently the code only accepts 'fasta')
#
# -ppi <path to directory containing the *.ppi file> # Default: 'Datasets/' 
#
# -v <0|1> # Default: 0 (not to include) It determines whether to include or not the name of teh project file in the features tag 
#
# -threads <number of threads to launch> # Default: '4' This parameter determines the number of used threads to distribute the run, as well as
#  the number of packages the data will be divided into.
#
# -in <Path to input directory> # Default: 'Datasets/Fasta_Protein_Format'
#
# -out <Path to output directory> # Default: 'Outputs'  
# 
# -proj <Path to projects directory> # Default: 'Projects' 
#
# -idl <Path to IDL-projects directory> # Default: 'Projects'
#
# -jar <Path to .jar file of ProtDcal> # Default: 'ProtDCal.jar'
#
# -Xmx <maximum heap space for Java VM> # Default: '1000m'  Please not that MUST be an space between the option and its value e.g. '-Xmx 2000m'
#
# -type <proj|idl|ppi|ppi_idl> # It determines the type of calculation
#
# -gz <true|false> # Default: 'true' It determines whether to zip the results file.
#.\n";
close HELP;
}
