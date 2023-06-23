#!/usr/bin/perl -w
use strict;

my $usage = "$0  <output directory >

Currently the EpaCons.bntree model is copied and used for the CNC.bntree model

\n";


die $usage unless @ARGV == 1; 

my $output_dir = shift @ARGV;



&copy_model;



#############################################################

sub copy_model {

    my $start = 1;
    my $line = "";
    my $header;
    my @list;
    if (open (I, "< $output_dir/EpaCons.bntree") ) {
        open (CNC, "> $output_dir/CNC.bntree");
        while ($line = <I>) {
            chomp $line;
            if ($start) {
                $start = 0;
                @list = split (/\s+/, $line);
                die "ERROR on first line of EpaCons.bntree\n" unless ($list[0] eq "EpaCons");
                $header = "CNC" . substr ($line, 7, length ($line) - 7);
                printf CNC "%s\n", $header;
            } else { 
                printf CNC "%s\n", $line;
            }
        }
    }

                
    printf STDERR "Faking CNC state  with EpaCons.bntree\n";
}

