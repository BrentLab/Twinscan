#!/usr/bin/perl -w
use strict;

my $usage = "$0  <output directory > <old parameter file to split> <new parameter file name> <program directory>

Currently the EpaCons.bntree model is copied and used for the CNC.bntree model

This requires an old zhmm parameter file. The new phylogenetic parameters will replace just the old 
phylogenetic parameters to create a new nscan parameter file.

This wrapper acutally runs 3 other programs

	fake_CNC.pl
	combine_params.pl
	combine_zhmm.pl

which finishes construction of the new parameter file.

\n";


die $usage unless @ARGV == 4; 

my $param_dir = shift @ARGV;
my $split_zhmm = shift @ARGV;
my $nscan_name = shift @ARGV;
my $prog_dir = shift @ARGV;

# copy EpaCons.bntree to CNC.bntree to fake a CNC state
printf STDERR "making CNC state from EpaCons\n";
#system ("/home/bl/rhb/ssgross/scripts/pipeline/fake_CNC.pl $param_dir");
system ("$prog_dir/fake_CNC.pl $param_dir");

# combine phylogenetic trees for each feature into a single phylogenetic paramter file
printf STDERR "combining bntree files into a phylogenetic parameter file\n";
#system ("/home/bl/rhb/ssgross/scripts/pipeline/combine_params.pl $param_dir > $param_dir/phylo_params");
system ("$prog_dir/combine_params.pl $param_dir > $param_dir/phylo_params");

# split an old zhmm nscan parameter file into pieces so that the new phylogenetic parameters can replace the old ones
printf STDERR "splitting the old nscan parameter file into blocks\n";
system ("mkdir $param_dir/split_old_zhmm");
#my $command = "/home/bl/rhb/ssgross/scripts/pipeline/split_zhmm.pl  $split_zhmm $param_dir/split_old_zhmm";
my $command = "$prog_dir/split_zhmm.pl  $split_zhmm $param_dir/split_old_zhmm";
printf STDERR "COMMAND\n";
printf STDERR "$command\n";
system ("$command");

# combine the new phylogenetic parameters and the rest of the old zhmm file into a new nscan parameter file
printf STDERR "combing the old blocks from the nscan parameter file with the new phylogenetic parameters\n";
#system ("/home/bl/rhb/ssgross/scripts/pipeline/combine_zhmm.pl $param_dir/split_old_zhmm $param_dir/phylo_params $param_dir/$nscan_name");
system ("$prog_dir/combine_zhmm.pl $param_dir/split_old_zhmm $param_dir/phylo_params $param_dir/$nscan_name");
