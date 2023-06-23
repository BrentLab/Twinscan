#!/usr/bin/perl -w
use strict;
use lib "$ENV{TWINSCAN}/lib";
use BPlite;
use Getopt::Std;
use vars qw($opt_u $opt_g $opt_t $opt_x $opt_s $opt_p);
getopts('ugtxsp');
my $usage = "
usage: conseq.pl [options] <fasta file> <blast file> <...blast file>
  -g include gaps 
  -u include unaligned (0 = mismatch, 1 = match, 2 = unaligned) 
  -t include transitions and transversions
  -x use whole HSP (0 = unaligned, 1 = hsp) 
  -s use symbols instead of numbers
  -p sort HSPs by percent identity (default is score)
";

die $usage unless @ARGV >= 2;
my ($fasta_file, @blast_file) = @ARGV;

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
	':' => 1, # mismatch
	'|' => 2, # match
	'.' => 3, # unaligned
	'/' => 4, # transition
	'-' => 5, # gap
);


open(FASTA, $fasta_file) or die;
my ($def, @seq) = <FASTA>;
chomp @seq;
close FASTA;
my $seq = join("", @seq);

# collect all the HSPs - going to sort by best HSP eventually
#print STDERR "Reading";
my @hsp;
my @database;
foreach my $blast_file (@blast_file) {
#	print STDERR ".";
    my $rezip;
    if(-e $blast_file){
	$rezip = 0;
    }
    elsif(-e "$blast_file.gz"){
	unless(system("gunzip $blast_file.gz") == 0){
	    die "Couldn't unzip $blast_file.gz\n";
	}
	$rezip = 1;
    }
    else{
	die "Couldn't find $blast_file\n";
    }   
    open(BLAST, $blast_file) or die;
    my $blast = new BPlite::Multi(\*BLAST);
    while(my $report = $blast->nextReport) {
		push @database, $report->database;    
		while(my $sbjct = $report->nextSbjct) {
		    while(my $hsp = $sbjct->nextHSP) {
			push @hsp, $hsp;
		    }
		}
    }    
    close BLAST;


    if($rezip){
	system("gzip $blast_file");
    }
}
# sort by score
#print STDERR "sorting";
if ($opt_p) {
	@hsp = sort {$b->percent <=> $a->percent or $b->score <=> $a->score} @hsp;
}
else {
	@hsp = sort {$b->score <=> $a->score or $b->percent <=> $a->percent} @hsp;
}
#print STDERR ".";

#print STDERR "conseq";
my @G; # genomic index - contains the conservation symbols
foreach my $hsp (@hsp) {
	my ($begin, $end) = ($hsp->qb, $hsp->qe);
	my ($qs, $hs, $as) = ($hsp->qa, $hsp->sa, $hsp->as);
	if ($begin > $end) {
		($begin, $end) = ($end, $begin);
		$qs = reverse $qs;
		$hs = reverse $hs;
		$as = reverse $as;
	}

	my $conseq = "";
	for (my $i=0;$i<length($qs);$i++) {
		my $qb = substr($qs, $i, 1);
		my $hb = substr($hs, $i, 1);
		my $ab = substr($as, $i, 1);
		if    ($qb eq '-') {next}           # skip query gaps, no indexing
		elsif ($hb eq '-') {$conseq .= '-'} # sbjct gap
		elsif ($ab eq ' ') {
			if (exists $Transition{$qb}{$hb}) {$conseq .= '/'}
			else                              {$conseq .= ':'}
		} # mismatch
		else {$conseq .= '|'} # match
		
	}
	
	for (my $i=0;$i<length($conseq);$i++) {
		my $index = $i + $begin -1;
		my $prev  = $G[$index];
		my $char = substr($conseq, $i, 1);
		
		# precedence rule - if defined by a better HSP, leave it alone
		if (not defined $prev) {$G[$index] = $char}
	}
}
#print STDERR ".";
	

# change undefined values to the unaligned symbol
for(my $i=0;$i<length($seq);$i++) {
	$G[$i] = '.' unless defined $G[$i]; # unaligned symbol is .
}
if (length($seq) != scalar(@G)) {
	print length($seq), "\n", scalar(@G), "\n";
	die "conseq length difference! ", length($seq) , " != ", scalar(@G);
}


# output
my $conseq = join("", @G);
if (!$opt_t) {$conseq =~ tr/\//:/} # change transition to mismatch
if ($opt_x) {$conseq =~ tr/\-/|/}  # change gaps	   to match
if ($opt_x)  {$conseq =~ tr/:/|/}  # chance mismatch   to match
if (!$opt_g) {$conseq =~ tr/\-/:/} # change gaps       to mismatch
if (!$opt_u) {$conseq =~ tr/\./:/} # change unaligned  to mismatch

my %symbol=(
	':' => 0, # mismatch
	'|' => 0, # match
);
if ($opt_u){$symbol{'.'}=0;}
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

	my $code = "\$conseq =~ tr[$symbols][$numbers]";
	eval "$code";
}

print ">Informant database(s): ";
for (my $i = 0; $i <= $#database; $i++) {
	print "$database[$i]\t";
}
print "\n";
print $conseq, "\n";




#print STDERR "\n";


sub symbolic {
	return $Order{$a} <=> $Order{$b}
}
