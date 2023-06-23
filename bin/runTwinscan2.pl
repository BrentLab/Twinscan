#!/usr/bin/perl -w
##################################################################
#   Copyright (c) 2000-2004 Washington University, St. Louis, MO    
#   All Rights Reserved    
#   Send all comments to twinscan@cse.wustl.edu
#   Driver for Twinscan 2.03
#   Based on the driver for Twinscan 2.0
##################################################################
$| = 1; # don't buffer STDOUT (otherwise output is hidden during long ops)

use strict;
use File::Basename;
use Getopt::Long;

my ($opt_d, $opt_blat, $opt_B, $opt_RM, $opt_r, $opt_x, $opt_v);
GetOptions(
	   "d=s" => \$opt_d,
           "B=s" => \$opt_B,
	   "RM=s" => \$opt_RM,
           "blat=s" => \$opt_blat,
	   "r=s" => \$opt_r, 
	   "x" => \$opt_x,
	   "v" => \$opt_v);
my $usage = "
runTwinscan2.pl - An example script for the Twinscan 2.03 process

usage: $0 [options] <target> <informant> [<ESTdb>]

Target Options:
  -r  <file>                Twinscan zoe HMM Parameter file
  -B  <file>                BLAST parameter file
  -b  <blat>                BLAT parameter file

REPEATMASK OPTION:
-RM RepeatBlaster's option, please place in double quotes, such as \"-ar\" 
    for Arabidopsis, \"-mus\" for mouse, (default human)

    If RepeatBlaster is not available, will switch automatically 
    to RepeatMasker

Informant Options:
  -x  format informant with xdformat  [default assumes blast db in path]

Output Options:
  -d  Specify Directory for output file [default is current directory]

EST Options: 
 <ESTdb> is a list of EST files in FASTA format. 

Mandatory Files (will be created if they don't exist):
  <target>.{masked, blast, conseq, twinscan, gff}
";

##########################################################
# $TWINSCAN_HOME should ideally be set with the TWINSCAN #
# environment variable.                                  #
##########################################################

my $TWINSCAN_HOME = "";
if ( $ENV{TWINSCAN} && -e $ENV{TWINSCAN}){
    $TWINSCAN_HOME = $ENV{TWINSCAN};
}
else{
    my $TWINSCAN_BIN_PATH="./";
    $TWINSCAN_HOME="..";
    my $tmp_command= $0;
    if ($tmp_command =~ /(^\S+\/)(.+)$/){
	$TWINSCAN_BIN_PATH = $1;
	if ($TWINSCAN_BIN_PATH =~ /(^\S+)\/bin/)
        {
	    $TWINSCAN_HOME =$1;
	}
	elsif($TWINSCAN_BIN_PATH =~ /^bin/){
	    $TWINSCAN_HOME =".";
	}
	elsif(($TWINSCAN_BIN_PATH eq "./")
	      ||($TWINSCAN_BIN_PATH eq ".")
	      ||($TWINSCAN_BIN_PATH eq "")){
	    $TWINSCAN_HOME ="..";
	}
	elsif($TWINSCAN_BIN_PATH eq "../"){
	    $TWINSCAN_HOME ="../..";
	}
    }
    $ENV{TWINSCAN} = "$TWINSCAN_HOME";
}

#  Arguements
die $usage unless @ARGV == 2 or @ARGV == 3;
my ($TARGET, $INFORMANT, $ESTdb) = @ARGV;

#  Paths Assigned
my $TWINSCAN_PATH 	= $ENV{TWINSCAN};
my $OUTPUT_PATH     	= $opt_d ? $opt_d : "$TWINSCAN_PATH/output";
my $REPEATBLASTER_EXE	= "RepeatBlaster";	# Format for local environment
my $REPEATMASKER_EXE	= "RepeatMasker";	# Format for local environment
my $BLASTN		= "blastn";		# Format for local environment
my $BLAT                = "blat";               # Format for local environment
my $XDFORMAT		= "xdformat";		# Format for local environment
my $PRESSDB		= "pressdb";		# Format for local environment
my $CONSEQ    		= "$TWINSCAN_PATH/bin/conseq.pl";
my $ESTSEQ    		= "$TWINSCAN_PATH/bin/estseq.pl";
my $TWINSCAN_EXE 	= "$TWINSCAN_PATH/bin/iscan";
my $TWINSCAN_GTF 	= "$TWINSCAN_PATH/bin/zoe2gtf";


#  Parameters
my $REPEATPARA = $opt_RM ? $opt_RM : "";
my $RUN_XDFORMAT     = $opt_x;
my $TWINSCAN_PARAM   = $opt_r ? $opt_r : "$TWINSCAN_PATH/parameters/worm_iscan_est-3889-genes-11-3-2004.zhmm";
my $BLAST_PARAM      = $opt_B ? $opt_B : "$TWINSCAN_PATH/parameters/Celegans.blast.param";
my $BLAT_PARAM      = $opt_blat ? $opt_blat : "$TWINSCAN_PATH/parameters/blat.param";
my $VERBOSE	 = $opt_v;


if (!(-e $OUTPUT_PATH)){ system "mkdir $OUTPUT_PATH";};

my $TARGET_PATH = "";
my $TARGET_FILE = $TARGET;
while($TARGET_FILE =~ /(^\S+\/)(.+)$/){
	$TARGET_PATH .= $1;
	$TARGET_FILE = $2;
}
 
open(PROGRESS,">$TARGET.twinscan.progress") or die "could not open $TARGET.twinscan.progress\n";

open(VERIFY, "<$TARGET") or die "could not open $TARGET\n";
my %v;
<VERIFY>;
while(<VERIFY>){
    if(/a/i){$v{a} = 1}
    if(/c/i){$v{c} = 1}
    if(/g/i){$v{g} = 1}
    if(/t/i){$v{t} = 1}
}
close(VERIFY);
my @acgt = keys %v;
exit_with_blank_gtf() unless @acgt == 4;

#######################
# input sequence check#
#######################
open(SEQ, "<$TARGET") || die "couldn't open $TARGET file\n";
my $dna="";
my $msg="";
my $description=">$TARGET";
my $count=0;
while (my $line = <SEQ>) {
    chomp $line;
    if ($line =~ /^>/) {
	$description = $line;
	$count++;
    }else{
	$dna .=$line; 
    }   
}

close SEQ;
if ($count >1){
    print "Error: $count sequences in the file. Twinscan allows only one sequence.\n";
    print PROGRESS "Error: $count sequences in the file. Twinscan allows only one sequence,\n";
    exit -1;
}

my $iupac =  ($dna =~ s/bdhkmrsuvwy/bdhkmrsuvwy/gi);

if ($iupac  > 0 ) {
   print "Non [A,C,G,T] IUPAC letters in sequence will be changed to N.  Other symbols will be deleted.";
   print PROGRESS "Non [A,C,G,T] IUPAC letters in sequence will be changed to N.  Other symbols will be deleted.";  
}

my $copy = $dna;
for ($dna){s/([^acgtn\s])//gi;}

my $prelen=length($copy);
for ($copy){ tr/acgtnACGTN//d;}
my $cutlen = length ($copy);
my $afterlen = $prelen -$cutlen;
if ($cutlen != 0){
   print "$cutlen characters deleted\n";
   print PROGRESS "$cutlen characters deleted\n";
}
if (($cutlen != 0)||($count==0)){   
   open(OUT, ">$OUTPUT_PATH/$TARGET_FILE.mod") || die "could not open $OUTPUT_PATH/$TARGET_FILE.mod  file\n";
   print OUT $description, "\n";
   print OUT $dna, "\n";
   close OUT;
   $TARGET_FILE= "$TARGET_FILE.mod";
   $TARGET = "$OUTPUT_PATH/$TARGET_FILE";
   $TARGET_PATH = $OUTPUT_PATH;
}



####################
# run RepeatMasker #
####################

if (-e "$OUTPUT_PATH/$TARGET_FILE.masked") {
	print PROGRESS "Did not mask, masked file already exists.\n";
} 

else {
	my $x = system("$REPEATBLASTER_EXE -s $REPEATPARA $TARGET");
	if ($x == -1) {
	    system("$REPEATMASKER_EXE -w -s $REPEATPARA $TARGET")==0 or die "repeat mask failed\n";
	   
	}
	if($OUTPUT_PATH ne $TARGET_PATH){	    
	    system("mv $TARGET.cat $OUTPUT_PATH/.");
	    system("mv $TARGET.masked $OUTPUT_PATH/.");
	    system("mv $TARGET.out $OUTPUT_PATH/.");
	    system("mv $TARGET.stderr $OUTPUT_PATH/.");
	    system("mv $TARGET.tbl $OUTPUT_PATH/.");	
	}
	print PROGRESS "Masked file created.\n";
}

############
# xdformat #
############
my $BLAST_DB;

if (-e "$INFORMANT.xnt") {
		$BLAST_DB = "$INFORMANT";
	}
	   
if (not defined $BLAST_DB) {
	my $x = system("$XDFORMAT");
	if ($x == -1) {$x = system("$PRESSDB");}
	else {
 		system("$XDFORMAT -n $INFORMANT") == 0 
			or die "ERROR: xdformat failed\n";
	}	
	if ($x == -1) { die "Error: neither xdformat nor pressdb found\n";}
	else { 
		system("$PRESSDB $INFORMANT > /dev/null") == 0 
				or die "PressDB Failed\n";
	} 
}

#############
# run BLAST #
#############
my $empty_blast = 0;

my $blast_param = get_blast($BLAST_PARAM);

if (-e "$OUTPUT_PATH/$TARGET_FILE.blast"){
	print PROGRESS "Did not create blast, file already exists\n";
}
    
else{
	my $out = "$OUTPUT_PATH/$TARGET_FILE";
	
	my $rc = system("$BLASTN $INFORMANT $out.masked $blast_param > $out.blast");
	if    ($rc == 256*23) {exit_with_blank_gtf()}
	elsif ($rc != 0)      {die "ERROR: blastn failed\n"}
	print PROGRESS "Created blast report.\n";
}


################################
# create conservation sequence #
################################
if(-e "$OUTPUT_PATH/$TARGET_FILE.conseq") {
	print PROGRESS "Did not conseq, file exists.\n";
}
else{
	system("$CONSEQ -u $TARGET $OUTPUT_PATH/$TARGET_FILE.blast > $OUTPUT_PATH/$TARGET_FILE.conseq ") 
		== 0 or die "ERROR: conseq failed\n";
	print PROGRESS "Conseq created from blast report.\n";
}


#############
# run BLAT #
#############



if (defined $ESTdb ) {
    warn("BLAT should be installed in your system in order to generate ESTseqs.\n");
    my $empty_blat = 0;
    my $blat_param = get_blast($BLAT_PARAM);

    if (-e "$OUTPUT_PATH/$TARGET_FILE.psl"){
	print PROGRESS "Did not create blat, file already exists\n";
    }
    else{

        my @ESTdbs = get_ESTdbs($ESTdb);

	for (my $i = 0; $i < @ESTdbs; $i ++) {
	    my $out = "$OUTPUT_PATH/$TARGET_FILE";
	    if ($i < 1) {
		my $rc = system("$BLAT  $out.masked ".$ESTdbs[$i]. " $blat_param $out.psl");
		if    ($rc == 256*23) {exit_with_blank_gtf()}
		elsif ($rc != 0)      {die "ERROR: blat failed\n"}
	    }
	    else {
		my $tmp="$out.psl.$$.$i";

		my $rc = system("$BLAT  $out.masked " . $ESTdbs[$i]. " -noHead $blat_param $tmp > /dev/null");
		if    ($rc == 256*23) {exit_with_blank_gtf()}
		elsif ($rc != 0)      {die "ERROR: blat failed\n"}
		system("cat $tmp >> $out.psl");
		system("rm -f $tmp");
	    }
	}
	print PROGRESS "Created blat report.\n";
    }
}

################################
# create ESTseqs               #
################################

if (defined $ESTdb ) {

    if(-e "$OUTPUT_PATH/$TARGET_FILE.estseq") {
	print PROGRESS "Did not estseq, file exists.\n";
    }
    else{
	system("$ESTSEQ  $TARGET $OUTPUT_PATH/$TARGET_FILE.psl > $OUTPUT_PATH/$TARGET_FILE.estseq")
	    == 0 or die "ERROR: estseq failed\n";
	print STDERR "$ESTSEQ  $TARGET $OUTPUT_PATH/$TARGET_FILE.psl > $OUTPUT_PATH/$TARGET_FILE.estseq\n"; 
	print PROGRESS "ESTseq created from blat report.\n";
    }

}


####################
# run Twinscan 2.0 #
####################
if (-e "$OUTPUT_PATH/$TARGET_FILE.twinscan") {
    print PROGRESS "Twinscan not run, file already exists.\n";
}
else{
    my $command = "$TWINSCAN_EXE -i $TWINSCAN_PARAM" ;
    if (defined $ESTdb) {
	$command = "$TWINSCAN_EXE -pe -i $TWINSCAN_PARAM" ;
    }
    $command .=" $OUTPUT_PATH/$TARGET_FILE.masked -c=$OUTPUT_PATH/$TARGET_FILE.conseq";

    if (defined $ESTdb) {
	$command .= " -e=$OUTPUT_PATH/$TARGET_FILE.estseq ";
    }

    system("$command > $OUTPUT_PATH/$TARGET_FILE.twinscan") == 0
	or die "ERROR: twinscan failed\n";
    print PROGRESS "Twinscan run.\n";
}


#run twinscan2gtf
system("$TWINSCAN_GTF $TWINSCAN_PARAM $OUTPUT_PATH/$TARGET_FILE.twinscan -s=$TARGET_FILE > $OUTPUT_PATH/$TARGET_FILE.gff") == 0
    or die "ERROR: twinscan2gtf failed\n";
print PROGRESS "Done.\n";
close(PROGRESS);



sub exit_with_blank_gtf{
    open(OUT, ">$OUTPUT_PATH/$TARGET_FILE.gff") or die;
    close(OUT);
    print "Exit with blank\n";
    print PROGRESS "All N sequence, twinscan not run.\n";
    print PROGRESS "Done.\n";
    close(PROGRESS);
    exit;
}

# Get BLAST parameters from file
sub get_blast{
    my($filename) = @_;
    open(PARA, "<$filename") || die "couldn't open $filename file\n";
    my $para;
    while (my $line = <PARA>) {
        chomp $line;
        $para .=" ". $line;
    }
    close PARA;
    return $para;
}


sub  get_ESTdbs {

    my ($f) = @_;
    my @ests; 
    open(DATA, $f) or die "can't open $f:$!\n";
    while (<DATA>) {
	chomp;
	if (/^\#/ or /^\s*$/) {
	    next;
	}
	else {
	    push(@ests, $_);
	}
    }
    close(DATA);
    return @ests;
}
