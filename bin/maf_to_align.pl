#!/usr/bin/perl -w
use strict;

my $TRUE = 1;
my $FALSE = 0;

my $START = 0;
my $STOP = 1;
my $SEQ = 2;

# 8 way alignment for human
#
#my @seq_names = qw (hg17 mm5 rn3 galGal2 canFam1 fr1 danRer1 panTro1);
#my %name_map = ("hg17" => 0,
#	"mm5"  => 1,
#	"rn3"  => 2,
#	"galGal2" => 3,
#	"canFam1" => 4,
#	"fr1" => 5,
#	"danRer1" => 6,
#	"panTro1" => 7);
#my $target_name = "hg17";

#printf STDERR "current alignment is 3 way for C. elegans\n";
# 3 way alignment C. elegans
#
#my @seq_names = qw (ce2 cb1 caeRei0 );
#my %name_map = ("ce2" => 0,
#       "cb1"  => 1,
#       "caeRei0"  => 2);
#my $target_name = "ce2";

printf STDERR "current alignment is 5 way for M. musculus\n";
# 5 way alignment M. musculus

my @seq_names = qw (mm5 rn3 hg17 canFam1 galGal2 );
my %name_map = ("mm5" => 0,
       "rn3"  => 1,
       "hg17"  => 2,
       "canFam1"  => 3,
       "galGal2"  => 4);
my $target_name = "mm5";

my $unaligned_char = ".";
my $gap_char = "_";

my $usage = "maf_to_align.pl <target sequence> <MAF file>\n";

@ARGV == 2 || die $usage;

my ($seq_file, $maf_file) = @ARGV;

print STDERR "Reading $seq_file\n";

#read in the target sequence
my $target_seq = "";
open (SEQ_FILE, "$seq_file");
<SEQ_FILE>;  #skip header
my $line_num = 0;
while (<SEQ_FILE>) {
  my $seq_line = $_;
  chomp ($seq_line);
  $target_seq .= uc($seq_line);
  $line_num++;
  if ($line_num % 100000 == 0) {
    print STDERR "Read $line_num lines\n";
  }
}
close (SEQ_FILE);
my $seq_length = length($target_seq);
print STDERR "Sequence $seq_file is $seq_length bp\n";

my @blocks;
my $block_num = 0;

open (MAF_FILE, "$maf_file");
my @maf_lines = <MAF_FILE>;
close (MAF_FILE);

my $j = 0;
while ($j<@maf_lines) {
  if ($maf_lines[$j] =~ /^a/) {  #alignment block
    my $start_coor;
    my $block_length;
    my @this_block;
    $j++;
    while ($maf_lines[$j] =~ /^s\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
		my $seq_name = $1;
		my $start = $2;
		my $size = $3;
		my $strand = $4;
		my $source_size = $5;
		my $text = $6;

		$seq_name =~ /(\S+)\.\S+/;
		$seq_name = $1;

		if (! exists($name_map{$seq_name})) {
		    die "Alignment block starting at $start: unrecognized sequence name $seq_name\n";
		}

		my @text_arr = split(//, uc($text));
		$this_block[$name_map{$seq_name}] = \@text_arr;

		if ($seq_name eq $target_name) {
		    $block_length = length($text);	
		    $start_coor = $start;
		    $blocks[$block_num][$START] = $start_coor;
		}
		$j++;
	    }
	   
	    my @target_arr = split(//, substr($target_seq, $start_coor, $block_length));
	    	    
	    for (my $u=0; $u<@seq_names; $u++) {
		$blocks[$block_num][$SEQ][$u] = "";
	    }

	    my $seq_coor = 0;
	    for (my $v=0; $v<$block_length; $v++) {
		if ($this_block[0][$v] ne "-") { #ignoring gaps in target sequence
		    #first check to make sure target sequence matches
		    #this block
		    if ($target_arr[$seq_coor] ne $this_block[0][$v] &&
			$target_arr[$seq_coor] ne "N") {
			print STDERR "Mismatch between sequence and alignment at block $start_coor, offset $v\n";
		    }
		    for (my $u=0; $u<@seq_names; $u++) {
			if (defined $this_block[$u]) {
			    if ($this_block[$u][$v] eq "-") {
				$blocks[$block_num][$SEQ][$u] .= $gap_char;
			    }
			    else {
				if ($u == 0) {
				    #include N's from target sequence in output file
				    $blocks[$block_num][$SEQ][$u] .= $target_arr[$seq_coor];
				}
				else {
				    $blocks[$block_num][$SEQ][$u] .= $this_block[$u][$v];
				}
			    }
			}
			else {  #no alignment for sequence $u in this block
			    $blocks[$block_num][$SEQ][$u] .= $unaligned_char;
			}
		    }
		    $seq_coor++;
		}
            }
            $blocks[$block_num][$STOP] = $blocks[$block_num][$START] + 
		length($blocks[$block_num][$SEQ][0]) - 1;
	    $block_num++;
	    if ($block_num % 5000 == 0) {
		print STDERR "Processed $block_num blocks\n";
	    }
	}
	$j++;
    }
    print STDERR "Processed $block_num alignment blocks in $j lines\n";

    #sort blocks by ascending stop coordinate
    my @unsorted = @blocks;
    @blocks = sort by_stop_coor @unsorted; 
    
    #now write out the data
    #first print the header
    print ">";
    for (my $i=0; $i<@seq_names; $i++) {
	print "$seq_names[$i] ";
    }
    print "\n";

    #print the target sequence
    print "$target_seq\n";
    my $target_length = length($target_seq);
    print STDERR "Printed target sequence, length $target_length\n";

    #print the aligned sequences
    for (my $i=1; $i<@seq_names; $i++) {
	my $length = 0;
	print STDERR "Printing sequence for $seq_names[$i]\n";
	my $current_coor = 0;
	for (my $block_num = 0; $block_num < @blocks; $block_num++) {
	    #start writing this block from the beginning or the end of the previous
	    #block, whichever is greater
	    my $offset;
	    if ($block_num == 0 || $blocks[$block_num][$START] > $blocks[$block_num-1][$STOP]) {
		$offset = 0;
	    }
	    else {
		$offset = $blocks[$block_num-1][$STOP] - $blocks[$block_num][$START] + 1;
	    }
	    if ($offset == 0) {
		my $unaligned_length = $blocks[$block_num][$START] - $current_coor;
		my $unaligned_str = $unaligned_char x $unaligned_length;
		print $unaligned_str;
		$length += length($unaligned_str);
	    }
	    else {
                #print STDERR "Overlap detected, $offset\n";
            }
	    if ($offset == 0) { #skip overlapping blocks for now
	      my $block_string = substr($blocks[$block_num][$SEQ][$i], 
		     $offset, length($blocks[$block_num][$SEQ][$i]) - $offset);
	      print $block_string;
	      $length += length($block_string);
	      $current_coor = $blocks[$block_num][$STOP] + 1;
            }
	}
	#Finally, the last string of unaligned characters
	my $unaligned_length = $seq_length - $current_coor;
	my $unaligned_str = $unaligned_char x $unaligned_length;
	print $unaligned_str;
	$length += length($unaligned_str);
	print "\n";
	print STDERR "Length $length, ending coor $current_coor\n";
    }

sub by_stop_coor {
    my @a_arr = @$a;
    my @b_arr = @$b;
    return $a_arr[$STOP] <=> $b_arr[$STOP];
}
