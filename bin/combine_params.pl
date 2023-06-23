#!/usr/bin/perl -w
use strict;

my $NAME = 0;
my $TYPE = 1;
my $LENGTH = 2;
my $FOCUS = 3;

my $order = 1;

# IntronCons model had its name changed to NCCons before Feb 10, 2005
my @models = (["NCCons", "BNTREE", 1, -1],
	      ["CNC", "BNTREE", 1, -1],
	      ["CDSCons", "BNTREE_CDS", 1, -1],
	      ["StartCodonCons", "BNTREE_ARRAY", 12, 6],
	      ["StopCodonCons", "BNTREE_ARRAY", 6, 0],
	      ["DonorCons", "BNTREE_ARRAY", 9, 3],
	      ["AcceptorCons", "BNTREE_ARRAY", 43, 39],
	      ["EpCons", "BNTREE", 1, -1],
	      ["EncCons", "BNTREE", 1, -1],
	      ["EaCons", "BNTREE", 1, -1],
	      ["EpaCons", "BNTREE", 1, -1],
	      ["UtrDonorCons", "BNTREE_ARRAY", 9, 3],
	      ["UtrAcceptorCons", "BNTREE_ARRAY", 43, 39]);

my $usage = "combine_params.pl <output directory>\n";

@ARGV == 1 || die $usage;

my ($param_dir) = @ARGV;

print "<PHYLOGENETIC_MODELS>\n\n";

for (my $i=0; $i<@models; $i++) {
    if ($models[$i][$TYPE] eq "BNTREE") {
	print STDERR "$models[$i][$NAME]\n";
	#just print the .bntree file
#       open (BNTREE_FILE, "$param_dir/$models[$i][$NAME].bntree");
	open (BNTREE_FILE, "$param_dir/$models[$i][$NAME].bntree") || die "$param_dir/$models[$i][$NAME].bntree IS MISSING!\n";
	while (<BNTREE_FILE>) {
	    print $_;
	}
	close (BNTREE_FILE);
    }
    elsif ($models[$i][$TYPE] eq "BNTREE_ARRAY") {
	print STDERR "$models[$i][$NAME]\n";
	#print header
	print "$models[$i][$NAME]\tBNTREE_ARRAY\t$order\t$models[$i][$LENGTH]\t$models[$i][$FOCUS]\n";
	for (my $j=0; $j<$models[$i][$LENGTH]; $j++) {
	    print STDERR "$models[$i][$NAME]: $j\n";
#           open (BNTREE_FILE, "$param_dir/$models[$i][$NAME]_$j.bntree");
	    open (BNTREE_FILE, "$param_dir/$models[$i][$NAME]_$j.bntree") || die "$param_dir/$models[$i][$NAME]_$j.bntree IS MISSING!\n";
	    while (<BNTREE_FILE>) {
		print $_;
	    }
	    close (BNTREE_FILE);
	}
    }
    elsif ($models[$i][$TYPE] eq "BNTREE_CDS") {
	print STDERR "$models[$i][$NAME]\n";
	#print header
	print "$models[$i][$NAME]\tBNTREE_CDS\t$order\n";
	for (my $j=0; $j<3; $j++) {
	    print STDERR "$models[$i][$NAME]: $j\n";
#           open (BNTREE_FILE, "$param_dir/$models[$i][$NAME]_$j.bntree");
	    open (BNTREE_FILE, "$param_dir/$models[$i][$NAME]_$j.bntree") || die "$param_dir/$models[$i][$NAME]_$j.bntree IS MISSING!\n";
	    while (<BNTREE_FILE>) {
		print $_;
	    }
	    close (BNTREE_FILE);
	}
    }
    print "\n";
}
