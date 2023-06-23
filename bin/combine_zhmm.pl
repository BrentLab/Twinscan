#!/usr/bin/perl -w
use strict;

# program builds a zhmm file

my ($i,$j,$k,$m,$n);

my $usage = "usage: $0 <old zhmm file dir> <full-path phylogenetic parameter file> <full-path new parameter file, use nscan.name.zhmm >

This program uses and an old nscan parameter file and builds a new one that is identical to the old one except it has
new phylogenetic parameters
 \n";
die $usage unless @ARGV == 3;

my $zhmm_file = $ARGV[0];
my $phylo_file = $ARGV[1];
my $new_param_file = $ARGV[2];

system ( "cat $zhmm_file/header.zhmm > $new_param_file" );
system ( "cat $zhmm_file/STATES.zhmm >> $new_param_file" );
system ( "cat $zhmm_file/STATE_TRANSITIONS.zhmm >> $new_param_file" );
system ( "cat $zhmm_file/STATE_DURATIONS.zhmm >> $new_param_file" );
system ( "cat $zhmm_file/SEQUENCE_MODELS.zhmm >> $new_param_file" );
    
system ( "cat $phylo_file >> $new_param_file" );

system ( "cat $zhmm_file/GTF_CONVERSION.zhmm >> $new_param_file" );
