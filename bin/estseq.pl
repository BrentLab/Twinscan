#!/usr/bin/perl 
use strict;
use Getopt::Std;
use vars qw($opt_a $opt_E $opt_l $opt_g $opt_t $opt_x $opt_s 
            $opt_p $opt_S $opt_P $opt_v $opt_u $opt_w);
getopts('aEfugtxsSplP:vw');
my $usage = "
usage: estseq.pl [options] <fasta file> <blat file> <...blat file>

generate est seq from blat alignments (psl format)


  -s                    output sequence(s) of symbols instead of sequences 
                        of numbers.
  -P <match_percent>    The minimum match percentage for an EST alignment to 
                        be included.
                        default value is 0.95 (95%).

  -S                    use Spliced EST alignments only.
  -u                    output 4-symble EST alignment sequence. 

  -a                    input file is UCSC all field file. 
                        default is blat output with ESTs as query and 
                        the target genomic sequence as database.

  -l                    <blat_file> is a list of blat alignment file. 
                        default is a blat alignment file.
";

die $usage unless @ARGV >= 2;
my ($fasta_file, @blat_file) = @ARGV;

my $mpercentage = 0.95;
if ($opt_P) {
    $mpercentage = $opt_P;
}


my %Transition = (
	A => {'G'=>1, 'R'=>1},
	C => {'T'=>1, 'Y'=>1},
	G => {'A'=>1, 'R'=>1},
	T => {'C'=>1, 'Y'=>1},
	R => {'A'=>1, 'G'=>1},
	Y => {'C'=>1, 'T'=>1}
);

# note that this order may get changed later
my %Order = (
	'N' => 1, # unaligned region 
	'E' => 2, # transcription
	'I' => 3, # unaligned region between two aligned regions, i.e. intron region
	'U'=> 3.5,# uncertain
	'/' => 4, # transition
	'-' => 5, # 
);

#my %seq = get_seqs($fasta_file);
my %seq_length = get_seq_lengths($fasta_file);

#----------------------------------------------------------------
# collect all the HSPs - going to sort by best HSP eventually
# organise the HSPs by subject
#----------------------------------------------------------------

my @ests;
my %ests;
my $f_n = 0;
my @blat_files;
if (!$opt_l) {
    @blat_files = @blat_file;
}
else {
    foreach my $blat_file (@blat_file) {
	open(DATA, $blat_file) or die "can't open file $blat_file:$!\n";
	 @blat_files=(@blat_files, <DATA>); 
	close(DATA);
    }
    chomp(@blat_files);
}


foreach my $blat_file (@blat_files) {
    $f_n++;
    print STDERR "Blat_file:$blat_file\n";
    my $blat ;
    if ($opt_a) {
	$blat =  blat_align($blat_file, 1);
    }
    else {
	$blat =  blat_align($blat_file);
    }


    my %regs  = %$blat;

    foreach my $sbj(keys %regs ) {
	my @als = @{$regs{$sbj}};

	@als = sort{$$b{matched_size} <=> $$a{matched_size}} @als;
	for (my $i = 0; $i < @als; $i ++) {
	    my %al = %{$als[$i]};
	    my $dbname = $al{dbname};

	    my ($strand) = $al{strand};
	    my $match_size = $al{matched_size};
	    my $t_size = $al{t_size};
		my $block_size  = $al{b_size};
	    my @block_sizes = split(/,/, $block_size);
	    my $b_n = @block_sizes;

	    my $q_begs=$al{q_start};
	    if ($q_begs =~ /^\s*(\S+),\s*$/) {
		$q_begs= $1;
	    }    
	    print STDERR "\n$i: $dbname, $strand, matched $match_size, block_size:$t_size $block_size start: ", $al{t_start}, , " qeury size:". $al{q_size}."\n" if $opt_v;

	    if ($al{matched_size}/$al{q_size} >= $mpercentage) {
		    if ($opt_S and $b_n < 2) {
		next;
	    }
		else {
		    #print STDERR "est: $als[$i]\n";
		    push(@{$ests{$dbname}}, $als[$i]);
		}
	    }
	}
    }
}



#------------------------------------------------------------
#  construct the est_seq from the hsps organized by subjects
#  
#  the region without EST alignment is N-> 0;
#  the region only one EST alignment is E-> 1;
#  if two or more ESTs aligned to the same region, the region not covered by all of alignments is uncertain, U-> 3, or if the window centered in this position with size 5 has mismatch 
#  a region between two E region is set to I->2 if it is not Uncertain.
#
#------------------------------------------------------------


my $i = 0;
foreach my $seq(keys %seq_length) {
    print STDERR "seq: $seq, length: ", $seq_length{$seq}{seq_length}, "\n";
    if (defined $ests{$seq}) {
	@ests = @{$ests{$seq}};
    }
    else {
	@ests = ();
    }
    generate_estseq(\@ests, $seq_length{$seq}{seq_length}, $seq);
}







sub generate_estgtf{ #output gtf files

    my ($ests, $seq_length, $seq_name) = @_;
    my @ests = @$ests;
    #my @G; # genomic index - contains the conservation symbols
    my $G; 
    my $jj = 0; 
    my $ests_n = @ests;
    print STDERR "total est numbers: $ests_n\n";
    my $test = 0; 

    print ">$seq_name\n";
    my $gn_id = 0; 

    foreach my $est(@ests) {
	 $jj ++;
	 print STDERR "$jj\t " if $jj%1000 == 1;
	 my %al = %$est;
	 my $block_size  = $al{b_size};
	 my @block_sizes = split(/,/, $block_size);
	 my $b_n = @block_sizes;

	 my $q_begs = $al{t_start};

         if ($q_begs =~ /^\s*(\S+),\s*$/) {
	     $q_begs= $1;
	 }
	 my @q_begs = split(/,/, $q_begs);
	 my $r_n = @q_begs;
	 if ($r_n != $b_n ) {
	     die  "error: block size $b_n and the start point number $r_n  are not equal.\n";
	 }

	 for (my $j = 0;  $j < $r_n; $j ++) {
	     my $b = $q_begs[$j];
	     my $l = $block_sizes[$j];

	     print "$seq_name\tEST\tCDS\t$b+1\t", $b+$l, "\t.\t", $al{strand},"\t.\tgene_id \"$seq_name\_$gn_id\"; transcript_id \"$seq_name\_$gn_id\";\n";
	 }

     }


 }




sub generate_estseq{
    my ($ests, $seq_length, $seq_name) = @_;
    my @ests = @$ests;
    #my @G; # genomic index - contains the conservation symbols
    my $G; 
    my $jj = 0; 
    my $ests_n = @ests;
    print STDERR "total est numbers: $ests_n\n";
    my $test = 0; 
    for (my $i = 0; $i < $seq_length; $i ++) {
	substr($G, $i, 1) = 'N';
    }
    foreach my $est(@ests) {
	 $jj ++;
	 print STDERR "$jj\t " if $jj%1000 == 1;
	 my %al = %$est;
	 my $block_size  = $al{b_size};
	 my @block_sizes = split(/,/, $block_size);
	 my $b_n = @block_sizes;

	 my $q_begs = $al{t_start};

         if ($q_begs =~ /^\s*(\S+),\s*$/) {
	     $q_begs= $1;
	 }
	 my @q_begs = split(/,/, $q_begs);
	 my $r_n = @q_begs;
	 if ($r_n != $b_n ) {
	     die  "error: block size $b_n and the start point number $r_n \n";
	 }

	 for (my $j = 0;  $j < $r_n; $j ++) {
	     my $b = $q_begs[$j];
	     my $l = $block_sizes[$j];

	     for (my $p = $b; $p < $b + $l; $p ++) { # for exons
		 if (substr($G, $p, 1) eq 'N') {
		     substr($G, $p, 1) =  'E';
		 }
		 else {
		     if (substr($G, $p, 1) ne 'E' ) {  # $G[$i] eq  'I' or 'U'
			 substr($G, $p, 1) = 'U';
		     }
		 }
	     }

	     if ($j > 0 & $j <= $r_n) { # for introns
		 my ($i_b, $i_e) = ($q_begs[$j - 1] + $block_sizes[$j - 1], 
				    $q_begs[$j] - 1);
		 if ($i_e - $i_b > 10) {
		     for (my $p = $i_b; $p <= $i_e; $p ++ ) {
			 if (substr($G, $p, 1) eq 'N') {
			     substr($G, $p, 1) = 'I';
			 }
			 else {
			     if (substr($G, $p, 1) ne 'I') { # $G[$i] eq 'E' or 'U'
				 substr($G, $p, 1) = 'U';
			     }
			 }
		     }
		 }
	     }
	 }

     }

     # change undefined values to the unaligned symbol
     for(my $i=0;$i<$seq_length;$i++) {
	 substr($G, $i, 1) = 'N' unless defined substr($G, $i, 1); # unaligned symbol is N
     }

     if ($seq_length != length($G)) {
	 print $seq_length, "\n", length($G), "\n";
	 die "conseq length difference! $seq_length != ", length($G);
     }

     #------------------------------------------------------------
     # output
     #------------------------------------------------------------
     #----------------  ???  -------------------------------------------
     if (!$opt_t) {$G =~ tr/\//:/} # change transition to mismatch
     if ($opt_x) {$G =~ tr/\-/E/}  # change gaps	   to match
     if ($opt_x)  {$G =~ tr/:/U/}  # chance mismatch   to match
     if (!$opt_g) {$G =~ tr/\-/U/} # change gaps       to mismatch
     if (!$opt_u) {$G =~ tr/U/N/} # change uncertain  to N
     

     my %symbol=(
		 'U' => 0, # mismatch
		 'E' => 0, # match
		 'N' => 0, #unaligned
		 'I' =>0   #intron
		);
     if ($opt_u){$symbol{'U'}=0;}
     if ($opt_g){$symbol{'-'}=0;}
     if ($opt_t){$symbol{'/'}=0;}
     
     my $n = 0;
     foreach my $char (sort symbolic keys %symbol) {
	 $symbol{$char} = $n;
	 $n++;
     }


     if (!$opt_s) {
#	print STDERR "translating.";
	my @key   = sort symbolic    keys %symbol;
	my @value = sort {$a <=> $b} values %symbol;
	my $symbols  = join("", @key);
	my $numbers  = join("", @value);

	 my $code = "\$G =~ tr[$symbols][$numbers]";
	 eval "$code";
     }


    print ">$seq_name\n$G\n";
 }



#------------------  subroutines ---- --------------




sub symbol_seq{
    my ($symb, $length) = @_;
    my $seq = "";
    for (my $i = 0; $i < $length; $i ++) {
	$seq.= "$symb";
    }
    return $seq;
}



sub symbolic {
	return $Order{$a} <=> $Order{$b}
}




sub overlap {
    my ($bb, $ee, $regs) = @_;


    if (defined @$regs) {
	for (my $i = 0; $i < @$regs; $i ++) {
	    #my ($r_b, $r_e) = ($$regs[$i]->qb, $$regs[$i]->qe);
	    my ($r_b, $r_e) = ($$regs[$i]->sb, $$regs[$i]->se);

	    if (($bb <= $r_b && $ee >= $r_b) or
		($bb <= $r_e && $ee >= $r_e ) or
		($r_b <= $bb and $ee<=$r_e)) {
		return 1;
	    }
	}

    }

    return -1;

}

sub get_seq_lengths{
    my ($file) = @_;
    my %seqs;
    open(FASTA, $file) or die "can't open file $file:$!\n";
    my $seq_length = 0;
    my $start = 0;
    my $seq_name;

    while (<FASTA>) {
	chomp;
	if (/^\s*$/) {
	    next;
	}
	elsif (/^>(\S+)/) {
	    if ($start > 0) {
		$seqs{$seq_name}{seq_length} = $seq_length;
		$seq_length = 0;
	    }
	    else {
		$start = 1;
	    }
	    $seq_name = $1;
	    $seqs{$seq_name}{def} = $_;

	}
	else {
	    $seq_length += length($_);
	}
    }
    if (defined $seq_name) {
	$seqs{$seq_name}{seq_length} = $seq_length;
    }
    close(FASTA);

    return %seqs;
}


sub get_seqs{
    my ($file) = @_;
    my %seqs;
    open(FASTA, $file) or die "can't open file $file:$!\n";
    my $seq = "";
    my $start = 0;
    my $seq_name;
    while (<FASTA>) {
	chomp;
	if (/^\s*$/) {
	    next;
	}
	elsif (/^>(\S+)/) {
	    if ($start > 0) {
		$seqs{$seq_name}{seq} = $seq;
	    }
	    else {
		$start = 1;
	    }
	    $seq_name = $1;
	    $seqs{$seq_name}{def} = $_;
	}
	else {
	    $seq .= $_;
	}
    }
    if (defined $seq_name) {
	$seqs{$seq_name}{seq} = $seq;
    }
    close(FASTA);

    return %seqs;
}




sub blat_align{ # get the aligned region from the blat result file, 
                    # return hash with subject as keys
    my ($f, $format) = @_;  # $format defined is for all field file downloaded from ucsc table browser

    my $start = 0;
    if ($f =~ /\.gz\s*$/) {
	open(DATA, "zcat $f |") or die "-can't open file $f:$!\n";
    }
    else {
	open(DATA, $f) or die "-can't open file $f:$!\n";
    } 
    my @a_regions;
    my %hit_regions;
    while (<DATA>) {
	chomp;
	my $line = $_;
	if (/^\s*$/) {
	    next;
	}
	elsif (/^---------/) {
	    $start = 1;
	}
	else {
	    #print STDERR "line: $line\n";
	    if (defined $format) {
		($line) = $line =~/^\S+\s+(.*)$/;
	    }
	    #print STDERR "#line: $line\n";

	    if ($line =~ /^(\d+)\s+(\S+\s+){8}(\S+)\s+(\S+)\s+(\S+\s+){3}(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/ ) {
		#      1       2          3       4          5       6       7       8      9       10     11
		my ($match_size, $strand, $query, $q_size, $dbname, $t_size, $t_beg,
		$t_end, $b_count, $b_size, $q_start, $t_start ) =
		( $1,        $2,      $3,      $4,      $5,     $6, $7,
		  $8,        $9,      $10,       $11,    $12   );

		($strand) = $strand =~ /^(\S+)\s*$/;

		#if ( $query =~ /^gi\|([^\|]+)/) {
		if ( $query =~ /^(\S+)/) {
		    $query = $1;
		}
		($dbname) = $dbname =~ /^(\S+)\s*$/;
		if ($b_size =~ /^(\S+),\s*$/) {
		    $b_size = $1;
		}

		my %a_region = (dbname=>$dbname, matched_size=>$match_size,
				strand=>$strand, beg=>$t_beg, end=>$t_end, 
				query => $query, b_size=>$b_size, t_size =>$t_size, q_size => $q_size,
				q_start=>$q_start, t_start => $t_start);
		push(@{$hit_regions{$dbname}}, \%a_region); # the only line different with blat_align()
		push(@a_regions, \%a_region);
	    }
	}
    }
    close(DATA);
    #return \@a_regions;
    return \%hit_regions;
}
