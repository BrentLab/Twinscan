#!/usr/bin/perl -w
use strict;

# taking input 8way mm5, rn3, hg17, canFam1 galGal2 to output 2way mm5, hg17
my $number_of_genomes = 8;
my @select_genomes = qw (1 5 );

my $usage = "$0 <input align file> <output align file>  \n";

@ARGV == 2 || die $usage;

my ($in_align_file) = $ARGV[0];
my ($out_align_file) = $ARGV[1];
die "input and output align files cannot be the same!\n" unless ($in_align_file ne $out_align_file); 


my ($i,$j,$k,$m,$n);

my $line;
my $char;
my $header =  "";

my @line_list;


open (IN_ALIGN, "< $in_align_file") || die "Couldn't open $in_align_file\n";
open (OUT_ALIGN, "> $out_align_file");



for ($i = 0; $i <= $number_of_genomes; $i++) {
 
    # skip rest of input if the last genome needed has been read in
    next if ($i > $select_genomes [-1]);

    $line = <IN_ALIGN>;
    chomp $line;
    $char = substr ($line, 0, 1) ;
    if ($char ne ">") {
        for ($j = 0; $j < @select_genomes; $j++) {
            if ( $i == $select_genomes [$j]) {
                printf OUT_ALIGN "%s\n", $line;
            }
        }
    } else {
        @line_list = split (/\s+/, $line);
        $header = $header . $line_list[0];
        for ($j = 1; $j < @select_genomes; $j++) {
            $header = $header . " " . $line_list[ $select_genomes[$j] - 1];
        }
        printf OUT_ALIGN "$header\n";
    }       
}
 
