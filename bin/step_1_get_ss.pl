#!/usr/bin/perl -w
use strict;
use GTF_Parser_UTR_SU;
use MultiRandom;

my $NC_SAMPLING_RATE = 0.01;

my $order = 1;


# human parameters for 4way alignment
##########################################################################################
my @human_chr_names = qw (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y); # target chromosomes
##########################################################################################



my $MC = 0;
my $WAM = 0;

my $CDS0 = 0;
my $CDS1 = 1;
my $CDS2 = 2;
my $NC_MC = 3;
my $EP_MC = 4;
my $EA_MC = 5;
my $ENC_MC = 6;
my $EPA_MC = 7;
my @mc_name = qw (CDSCons_0 CDSCons_1 CDSCons_2 NCCons EpCons EaCons EncCons EpaCons);

my $DONOR = 0;
my $ACCEPTOR = 1;
my $START_CODON = 2;
my $STOP_CODON = 3;
my $UTR_DONOR = 4;
my $UTR_ACCEPTOR = 5; 
my @wam_name = qw (DonorCons AcceptorCons StartCodonCons StopCodonCons  UtrDonorCons  UtrAcceptorCons);
my @wam_length = (     9,          43,            12,            6,          9,              43);


my $usage = "$0 
                <alignment directory> 
                <output directory> 
                < GTF directory>
                <number of genomes including target>
                <target and informant names> such as hg17,mm5,rn3,galGal2  in order, separated by commas, without any spaces

The chromosome names are hard-coded for human. If another specie is used as the target, then
the chromosome names must be changed in the code.
                \n";

@ARGV == 5 || die $usage;

my $align_dir = shift @ARGV;
my $out_dir = shift @ARGV;
my $gtf_dir = shift @ARGV;
my $n = shift @ARGV;
my $names_str = shift @ARGV;

$names_str =~ /(^[A-Za-z0-9]+)/;
my $genome = $1;
printf "genome $genome\n";
my @chr_names;
if ($genome eq "hg17") {
    @chr_names = @human_chr_names;
} else {
    die "ERROR: chromosome list must be hard-coded for any genome other than hg17\n";
}


my @mc;
my @wam;

my $mr;

for (my $arg=0; $arg<@ARGV; $arg++) {
  my $gtf_dir = $ARGV[$arg];
  for (my $chr_index=0; $chr_index<@chr_names; $chr_index++) {
      my $chr_name = "chr$chr_names[$chr_index]";
      my $gtf_file = "$gtf_dir/$chr_name.gtf";
      my $align_file = "$align_dir/$chr_name.align";
      $mr = MultiRandom->new($align_file);
      print STDERR "Processing file $gtf_file\n";
      
      my @unsorted_genes = @{GTF_Parser_UTR_SU::parse_gtf($gtf_file)};
      #sort genes by start coord
      my @genes = sort by_start @unsorted_genes;

      for (my $i=0; $i<@genes; $i++) {
	  my @this_gene = @{$genes[$i]};

	  if ($i != @genes-1 && rand() < $NC_SAMPLING_RATE) {
	      #get intergenic counts
	      my $this_gene_stop = $this_gene[@this_gene - 1][$STOP];
	      my $next_gene_start = $genes[$i+1][0][$START];
	      countIntergenic ($this_gene_stop+7, $next_gene_start-7);
	  }

	  for (my $j=0; $j<@this_gene; $j++) {
	      my $type = $this_gene[$j][$TYPE];
	      my $strand = $this_gene[$j][$STRAND];
	      my $start = $this_gene[$j][$START];
	      my $stop = $this_gene[$j][$STOP];    
	      my $frame = $this_gene[$j][$FRAME];

	      if (! (defined $genes[$i][$j][$TYPE])) {
		  print STDERR "Type undefined for $i, $j\n";
	      }

	      if ($type == $EINITU || $type == $EINITS) {
		  if ($strand eq "+") {
		      countWAM ($START_CODON, $start-6, $start+5, $strand);
		      countCDS ($frame, $start+6, $stop-3, $strand);
		      countWAM ($DONOR, $stop-2, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($DONOR, $start-6, $start+2, $strand);
		      countCDS ($frame, $start+3, $stop-6, $strand);
		      countWAM ($START_CODON, $stop-5, $stop+6, $strand);		      
		  }
	      }
	      elsif ($type == $EXON) {
		  if ($strand eq "+") {
		      countWAM ($ACCEPTOR, $start-40, $start+2, $strand);

##################### BUG FIX #############################################################################
#	RHB June 6, 2005
#                     countCDS ($frame, $start-3, $stop-3, $strand);
		      countCDS ($frame, $start+3, $stop-3, $strand);
##########################################################################################################

		      countWAM ($DONOR, $stop-2, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($DONOR, $start-6, $start+2, $strand);
		      countCDS ($frame, $start+3, $stop-3, $strand);
		      countWAM ($ACCEPTOR, $stop-2, $stop+40, $strand);
		  }
	      }
	      elsif ($type == $ETERM) {
		  if ($strand eq "+") {
		      countWAM ($ACCEPTOR, $start-40, $start+2, $strand);
		      countCDS ($frame, $start+3, $stop, $strand);		   
		      countWAM ($STOP_CODON, $stop+1, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($STOP_CODON, $start-6, $start-1, $strand);
		      countCDS ($frame, $start, $stop-3, $strand);
		      countWAM ($ACCEPTOR, $stop-2, $stop+40, $strand);
		  }
	      }
	      elsif ($type == $ESNGL) {
		  if ($strand eq "+") {
		      countWAM ($START_CODON, $start-6, $start+5, $strand);
		      countCDS ($frame, $start+6, $stop, $strand);
		      countWAM ($STOP_CODON, $stop+1, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($STOP_CODON, $start-6, $start-1, $strand);
		      countCDS ($frame, $start, $stop-6, $strand);
		      countWAM ($START_CODON, $stop-5, $stop+6, $strand);
		  }
	      }
	      elsif ($type == $INTRON && rand() < $NC_SAMPLING_RATE) {
		  if ($strand eq "+") {
		      countMC ($NC_MC, $start+6, $stop-39, $strand);
		  }
		  elsif ($strand eq "-") {
		      countMC ($NC_MC, $start+39, $stop-6, $strand);
		  }
	      }
	      elsif ($type == $EP) {
		   if ($strand eq "+") {
		      countMC ($EP_MC, $start, $stop-3, $strand);
		      countWAM ($UTR_DONOR, $stop-2, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($UTR_DONOR, $start-6, $start+2, $strand);
		      countMC ($EP_MC, $start+3, $stop, $strand);
		  } 
	       }
	      elsif ($type == $EA) {
		  if ($strand eq "+") {
		      countWAM ($UTR_ACCEPTOR, $start-40, $start+2, $strand);
		      countMC ($EA_MC, $start+3, $stop-6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countMC ($EA_MC, $start+6, $stop-3, $strand);
		      countWAM ($UTR_ACCEPTOR, $stop-2, $stop+40, $strand);
		  }
	      }
	      elsif ($type == $ENC) {
		  if ($strand eq "+") {
		      countWAM ($UTR_ACCEPTOR, $start-40, $start+2, $strand);
		      countMC ($ENC_MC, $start+3, $stop-3, $strand);
		      countWAM ($UTR_DONOR, $stop-2, $stop+6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countWAM ($UTR_DONOR, $start-6, $start+2, $strand);
		      countMC ($ENC_MC, $start+3, $stop-3, $strand);
		      countWAM ($UTR_ACCEPTOR, $stop-2, $stop+40, $strand);
		  }
	      }
	      elsif ($type == $EPA) {
		  if ($strand eq "+") {
		      countMC ($EPA_MC, $start, $stop-6, $strand);
		  }
		  elsif ($strand eq "-") {
		      countMC ($EPA_MC, $start+6, $stop, $strand);
		  }
	      }
	  }
      }
  }
}

#output sufficient statistics files

#Markov chains
for (my $i=0; $i<@mc_name; $i++) {
    open (MC_OUT, ">$out_dir/$mc_name[$i].ss");

    #count length and number of tuples
    my $length = 0;
    my $num_tuples = 0;
    foreach my $key (keys %{$mc[$i]}) {
	if ($key =~ /^[ACGT\._ ]*$/) {
	    $num_tuples++;
	    $length += $mc[$i]{$key};
	}
    }

    #print header
    print MC_OUT "NSEQS = $n\n";
    print MC_OUT "LENGTH = $length\n";
    print MC_OUT "TUPLE_SIZE = ".($order+1)."\n";
    print MC_OUT "NTUPLES = $num_tuples\n";
    print MC_OUT "NAMES = $names_str\n";
    print MC_OUT "ALPHABET = ACGT_.\n";
    print MC_OUT "NCATS = -1\n\n";

    my $j = 0;
    foreach my $key (keys %{$mc[$i]}) {
	if ($key =~ /^[ACGT\._ ]*$/) {
	    print MC_OUT "$j\t$key\t$mc[$i]{$key}\n";
	    $j++;
	}
    }
    close MC_OUT;
}


#WAMs
for (my $i=0; $i<@wam_name; $i++) {
    for (my $j=0; $j<$wam_length[$i]; $j++) {
	open (WAM_OUT, ">$out_dir/$wam_name[$i]_$j.ss");
	
	#count length and number of tuples
	my $length = 0;
	my $num_tuples = 0;
	foreach my $key (keys %{$wam[$i][$j]}) {
	    if ($key =~ /^[ACGT\._ ]*$/) {
		$num_tuples++;
		$length += $wam[$i][$j]{$key};
	    }
	}

	#print header
	print WAM_OUT "NSEQS = $n\n";
	print WAM_OUT "LENGTH = $length\n";
	print WAM_OUT "TUPLE_SIZE = ".($order+1)."\n";
	print WAM_OUT "NTUPLES = $num_tuples\n";
	print WAM_OUT "NAMES = $names_str\n";
	print WAM_OUT "ALPHABET = ACGT_.\n";
	print WAM_OUT "NCATS = -1\n\n";

	my $k = 0;
	foreach my $key (keys %{$wam[$i][$j]}) {
	    if ($key =~ /^[ACGT_\. ]*$/) {
		print WAM_OUT "$k\t$key\t$wam[$i][$j]{$key}\n";
		$k++;
	    }
	}
	close WAM_OUT;
    }
}

#########################################################################################

sub countWAM {
    my ($type, $start, $stop, $strand) = @_;

    if ($stop < $start) {
	return;
    }

    my @seq_arr;

    if ($strand eq "+") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getSeq($i, $start-$order, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    elsif ($strand eq "-") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getRevCompSeq($i, $start, $stop+$order);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }

    for (my $i=$order; $i<@{$seq_arr[0]}; $i++) {
	my $obs_str = "";
	for (my $k = $order; $k>=0; $k--) {
	    for (my $j=0; $j<$n; $j++) {
		$obs_str .= $seq_arr[$j][$i - $k];
	    }
	    if ($k != 0) {
		$obs_str .= " ";
	    }
	}
	if (exists $wam[$type][$i-$order]{$obs_str}) {
	    $wam[$type][$i-$order]{$obs_str}++;
	}
	else {
	    $wam[$type][$i-$order]{$obs_str} = 1;
	}
    }
}

#########################################################################################

sub countMC {
    my ($type, $start, $stop, $strand) = @_;

    if ($stop < $start) {
	return;
    }

    my @seq_arr;

    if ($strand eq "+") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getSeq($i, $start-$order, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    elsif ($strand eq "-") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getRevCompSeq($i, $start, $stop+$order);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    
    for (my $i=$order; $i<@{$seq_arr[0]}; $i++) {
	my $obs_str = "";
	for (my $k=$order; $k>=0; $k--) {
	    for (my $j=0; $j<$n; $j++) {
		$obs_str .= $seq_arr[$j][$i - $k];
	    }
	    if ($k != 0) {
		$obs_str .= " ";
	    }
	}
	if (exists $mc[$type]{$obs_str}) {
	    $mc[$type]{$obs_str}++;
	}
	else {
	    $mc[$type]{$obs_str} = 1;
	}
    }
}

#########################################################################################

sub countIntergenic {
    my ($start, $stop) = @_;

    if ($stop < $start) {
	return;
    }

    #intergenic sequence gets 1/2 counts on plus strand,
    #1/2 counts on minus strand
    my @seq;

    # plus strand
    for (my $i=0; $i<$n; $i++) {
	my $seq_ref = $mr->getSeq($i, $start-$order, $stop);
	$seq[$i] = $$seq_ref;
    }
    for (my $i=$order; $i<length($seq[0]); $i++) {
	my $obs_str = "";
	for (my $k=$order; $k>=0; $k--) {
	    for (my $j=0; $j<$n; $j++) {
		$obs_str .= substr ($seq[$j], $i - $k, 1);
	    }
	    if ($k != 0) {
		$obs_str .= " ";
	    }
	}
	if (exists $mc[$NC_MC]{$obs_str}) {
	    $mc[$NC_MC]{$obs_str} += 0.5;
	}
	else {
	    $mc[$NC_MC]{$obs_str} = 0.5;
	}
    }

    #minus strand
    for (my $i=0; $i<$n; $i++) {
	my $seq_ref = $mr->getRevCompSeq($i, $start, $stop+$order);
	$seq[$i] = $$seq_ref;
    }
    for (my $i=$order; $i<length($seq[0]); $i++) {
	my $obs_str = "";
	for (my $k=$order; $k>=0; $k--) {
	    for (my $j=0; $j<$n; $j++) {
		$obs_str .= substr ($seq[$j], $i - $k, 1);
	    }
	    if ($k != 0) {
		$obs_str .= " ";
	    }
	}
	if (exists $mc[$NC_MC]{$obs_str}) {
	    $mc[$NC_MC]{$obs_str} += 0.5;
	}
	else {
	    $mc[$NC_MC]{$obs_str} = 0.5;
	}
    }
}


#########################################################################################

sub countCDS {
    my ($frame, $start, $stop, $strand) = @_;

    if ($stop < $start) {
	return;
    }

    my @seq_arr;
    
    if ($strand eq "+") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getSeq($i, $start-$order, $stop);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }
    elsif ($strand eq "-") {
	for (my $i=0; $i<$n; $i++) {
	    my $seq_ref = $mr->getRevCompSeq($i, $start, $stop+$order);
	    my @seq_split = split (//, $$seq_ref);
	    $seq_arr[$i] = \@seq_split;
	}
    }

    for (my $i=$order; $i<@{$seq_arr[0]}; $i++) {
	my $type;
	if ($frame == 0) {
	    $type = $CDS0;
	}
	elsif ($frame == 1) {
	    $type = $CDS1;
	}
	elsif ($frame == 2) {
	    $type = $CDS2;
	}

	my $obs_str = "";
	for (my $k=$order; $k>=0; $k--) {	
	    for (my $j=0; $j<$n; $j++) {
		$obs_str .= $seq_arr[$j][$i - $k];
	    }
	    if ($k != 0) {
		$obs_str .= " ";
	    }
	}
	if (exists $mc[$type]{$obs_str}) {
	    $mc[$type]{$obs_str}++;
	}
	else {
	    $mc[$type]{$obs_str} = 1;
	}
	$frame = ($frame + 1) % 3;
    }
}

##################################################################################################

sub by_start {    
    my @gene_a = @$a;
    my @gene_b = @$b;

    return $gene_a[0][$START] <=> $gene_b[0][$START];
}
