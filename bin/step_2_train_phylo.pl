#!/usr/bin/perl -w
use strict;


my $usage = "$0 <train command> <tree> <model type> <output directory>

The train command refers to the C-code executable that trains the phylogenetic parameters. 
This is part of the download. 

The tree is the phylogenetic tree rerooted at the target. For example, for human as target with 
mouse, rat, and chicken as informants, 

human		hg17
mouse		mm5
rat		rn3
chicken		galGal2

the input tree looks like

	[[galGal2:4,[mm5:2,rn3:3]],none] 

where the number after the colon indicates the line where that specie  appears in the alignment file counting the fasta header
line as 0 and the target line as 1. 
 
The standard model type is R1.

\n";

@ARGV == 4 || die $usage;

my ($train_command, $tree, $model_type, $ss_dir) = @ARGV;
my $SCRIPT_DIR = $ss_dir;

opendir (SS_DIR, $ss_dir);
my @ss_files = grep(/\.ss$/, readdir(SS_DIR));
closedir (SS_DIR);


my $command_list = "";

for (my $i=0; $i<@ss_files; $i++) {
    $ss_files[$i] =~ /(\S+)\.ss$/;
    my $model_name = $1;

    my $ss_file = "$ss_dir/$ss_files[$i]";
    my $command = "$train_command $tree $model_type $ss_dir/$ss_files[$i] $model_name $ss_dir/$model_name.bntree";
    $command_list .= "$command\n";
}

open (SCRIPT_FILE, ">$SCRIPT_DIR/step_train_nscan");
print SCRIPT_FILE ("$command_list\n");
close (SCRIPT_FILE);

#####################################################################################
#
# The following command is used on the Brent Lab cluster at
# Washington University to submit batch jobs
# to estimate the phylogenetic tree parameters. It should be changed to the 
# correct form to submit batch jobs on your system.
# The 4000 indicates that it needs at most 4000 Mb of memory to run.

my $command = "enqueue $SCRIPT_DIR/step_train_nscan -P medium -mem 4000";

####################################################################################

system ($command);
