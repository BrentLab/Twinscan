#!/usr/bin/perl -w
use strict;

# program splits a zhmm file

my ($i,$j,$k,$m,$n);

my $usage = "usage: $0 <zhmm file> <output dir> \n";
die $usage unless @ARGV == 2;

my $zhmm_file = $ARGV[0];
my $output_dir = $ARGV[1];

    

my ($line);
#my $title = "header." . $zhmm_file;
my $title = "header.zhmm" ;
my $header;
my $first = 1;
my $body = "";
my $first_char;
my $last_char;

my @line_list;

open (Z, "< $zhmm_file") || die "cannot open $zhmm_file\n";

while ( defined ($line = <Z>)) {
    chomp $line;
    @line_list = split (/\s+/, $line);
    if (defined ($line_list[0])) {
#       printf "line defined\n";
        $first_char = substr ($line, 0, 1);
        $last_char = substr ($line, length($line) - 1, 1);
        if ($first_char eq "<") {
            &print_section;
            if ($last_char eq ">") {
#               $title = substr ($line, 1, length($line) -2) . "." . "$zhmm_file";
                $title = substr ($line, 1, length($line) -2) . ".zhmm" ;
            } else {
                die "header error $line\n";
            }
        }
        $body = $body . $line . "\n";
    } else {
#       printf "line not defined\n";
        if ($first) {
            $header = $header  . "\n";
        }
        $body = $body  . "\n";
    }
}
&print_section;

########################################################################################

sub print_section {

    printf STDERR "$title\n";
    open (P, "> $output_dir/$title");
    printf P "$body";
    $body = "";
    $first = 0;
}
