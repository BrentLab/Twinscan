#!/usr/bin/perl -w

##################################################################
#   Copyright (c) 2006 Washington University, St. Louis, MO      #
#   All Rights Reserved                                          #
#   Send all comments to Mani - mozhiyan@cse.wustl.edu           #
#   Driver script for Pairagon v0.95 and up                      #
#   Driver version July 2006                                     #
##################################################################

package pairagon;
use strict;
use BPlite;
use FAlite;
use File::Basename;
use Getopt::Long;

############################################
#          GLOBAL VARIABLES
#
# All the global variables are defined here. 
# This should be kept to a minimum. (Command
# line option variables excluded.)
############################################

my ($opt_a, $opt_b, $opt_d, $opt_e, $opt_g, $opt_n, $opt_r, $opt_s);
my ($opt_L, $opt_P);
my ($opt_lowQ, $opt_debug, $opt_unopt, $opt_unseeded, $opt_nocleanup);

# $opt_lowQ is defined for future compatibility with qPairagon.
# It is not used anywhere.

# Only THREE non-command-line-option global variables for now

my $SCRIPT_NAME = "runPairagon.pl";
my $PROGRESS_FH;
my $DEBUG = 0;

############################################
# Help message
############################################

my $usage = "
$SCRIPT_NAME - A driver script for the Pairagon program

usage: $SCRIPT_NAME [options] <cdna> [<genomic>] [<quality_seq_file>]

Program Options:
  -a       <mode>           direction of cDNA sequence to align to the genomic sequence (forward|reverse|both, default:both)
  -b       <blastdb>        blastdb of all genomic sequences (only with -g)
  -e       <directory>      directory where executables are located (default: \$HOME/bin)
  -g       <directory>      directory where genomic files are located (only with -b)
  -r       <file>           Pairagon's zoe HMM Parameter file (default: <exe_directory>/pairagon.zhmm)
  -s       <mode>           direction of sense strand with respect to the genomic sequence
                            (forward = +, reverse = -, cdna = same strand as the cDNA, default: cdna)
  -L       <length>         extra length on each side of blast aligned region to 
			    be extracted for fine alignment (default: 10000)
 --unopt                    do not use the optimized viterbi (default: false)
 --unseeded                 do not use seeded alignment (Global optimal alignment takes a VERY LONG time to run)
                            - CANNOT be used with -g and -b (default: false)


Output Options:
  -d       <string>         Specify Directory for output file (default: <pwd>/output)
  -n       <string>         Instance name (Alignment will be in <output_dir>/<instance_name>.estgen)
 --nocleanup                do not clean up temporary files at the end. 
                            Normally only the <instance_name>.{pair,estgen,progress} files remain.
  -P                        only generate a seed alignment file, then stop. 

Debug Options:
 --debug   <number>         prints out messages to debug
";

############################################
# Begin Execution
############################################

GetOptions(
	"a=s" => \$opt_a,
	"b=s" => \$opt_b,
	"d=s" => \$opt_d,
	"e=s" => \$opt_e,
	"g=s" => \$opt_g,
	"n=s" => \$opt_n,
	"r=s" => \$opt_r,
	"s=s" => \$opt_s,
	"L=s" => \$opt_L,
	"P"   => \$opt_P,
	"-debug=n"    => \$opt_debug,
	"-unopt"      => \$opt_unopt,
	"-unseeded"   => \$opt_unseeded,
	"-nocleanup"  => \$opt_nocleanup
);
if ((!$opt_b && @ARGV < 2) ||
    ( $opt_b && @ARGV < 1)) {
	die $usage;
}

run_me(@ARGV);
exit(0);

############################################
# End Execution
############################################


sub run_me {

	############################################
	#          run_me() VARIABLES
	#
	# All the run_me() variables are defined here. There
	# should not be any definition inside the
	# main function, except for iterators and
	# fully local variables (scope is within 10
	# lines). All run_me() variables are UPPER
	# case and local variables are lower case
	############################################

	# Arguments
	my ($CDNA, $GENOMIC, $QUALITY);

	# File names
	my ($QUALITY_FILE);
	my $HMM_FILE;
	my $SEED_FILE;
	my $PROGRESS_FILE;
	my $BLAST_DB;
	my $OUTPUT_DIR;
	my $EXE_DIR;
	my $GENOMIC_DIR;

	# Instance name
	my $NAME;

	# Seed alignment
	my $SEED_PROGRAM;

	# Program options
	my $EXTRA_GENOMIC;
	my $ALIGNMENT_MODE;
	my $SPLICE_MODE;
	my $OUTPUT_MODE;
	my $OPTIMIZED_MODE;
	my $SEEDED_MODE;
	my $CLEANUP;
	my $MAX_INTRON_LENGTH = 5000000;

	# Seeds
	my %SEEDS     = ();
	my @CDNA_DEFS = ();

	############################################
	# Grab arguments
	############################################

	($CDNA, $GENOMIC, $QUALITY) = @_;
	
	############################################
	# Set file names from command line options
	############################################


	$SEED_PROGRAM   = "WU-BLAST";   # Change this to whatever program you want 
					# and implement a function to generate seed
					# alignments from its output

	$NAME           = basename($CDNA);
	$EXE_DIR        = $ENV{"HOME"}."/bin/";
	$OUTPUT_DIR     = "output";
	$BLAST_DB       = undef;
	$GENOMIC_DIR    = undef;
	$EXTRA_GENOMIC  = 10000;
	$ALIGNMENT_MODE = undef;
	$SPLICE_MODE    = "cdna";
	$OPTIMIZED_MODE = "-o";
	$OUTPUT_MODE    = "-i";
	$SEEDED_MODE    = 1;
	$CLEANUP        = 1;

	if ($opt_a) {$ALIGNMENT_MODE = $opt_a;}
	if ($opt_b) {$BLAST_DB       = $opt_b;}
	if ($opt_d) {$OUTPUT_DIR     = $opt_d;}
	if ($opt_e) {$EXE_DIR        = $opt_e;}
	if ($opt_g) {$GENOMIC_DIR    = $opt_g;}
	if ($opt_n) {$NAME           = $opt_n;}
	if ($opt_s) {$SPLICE_MODE    = $opt_s;}
	if ($opt_L) {$EXTRA_GENOMIC  = $opt_L;}

	$HMM_FILE       = "$EXE_DIR/pairagon.zhmm";
	if ($opt_r) {$HMM_FILE       = $opt_r;}

	# Long command line options

	if ($opt_debug)    {$DEBUG          = $opt_debug;}
	if ($opt_unopt)    {$OPTIMIZED_MODE = "";}
	if ($opt_unseeded) {$SEEDED_MODE    = 0;}
	if ($opt_nocleanup){$CLEANUP    = 0;}

	$SEED_FILE      = "$OUTPUT_DIR/$NAME.seed";
	$PROGRESS_FILE  = "$OUTPUT_DIR/$NAME.progress";

	############################################
	# Check integrity of the command line
	############################################

	if (($opt_b && !$opt_g) || (!$opt_b && $opt_g)) {
		die "-b and -g options must be used together\n";
	}

	if ($opt_b) {
		if (defined($GENOMIC)) {
			print "# Ignoring $GENOMIC since -b and -g are set\n";
		}
		$GENOMIC = undef;
	}

	if ($SEEDED_MODE == 0 && ($opt_b || $opt_g)) {
		die "Unseeded alignment should be run on a single (and reasonably short) genomic sequence. Do not use -b and -g\n";
	}

	############################################
	# Check integrity of the input files
	############################################

	my $fasta_status;
	($fasta_status, @CDNA_DEFS) = get_valid_fasta_entries($CDNA);
	if ($fasta_status > 1 && $opt_b) {
		die "Cannot use multiple cDNA sequences in multiple genomic sequence mode";
	} elsif ($fasta_status == -1) {
		print "# Non [A,C,G,T,N] letters in input file $CDNA\n";
		print "# These will be converted to N\n";
		create_clean_fasta_entries($CDNA, "$OUTPUT_DIR/$NAME.cdna.fa.modified");
		$CDNA = "$OUTPUT_DIR/$NAME.cdna.fa.modified";
		($fasta_status, @CDNA_DEFS) = get_valid_fasta_entries($CDNA);
	}

	############################################
	# Don't do this for genomic sequences since
	# you will be reading and writing huge
	# files. Pairagon can convert non-ACGTN to N
	# internally now.
	############################################

	$fasta_status = undef; # End my scope

	open(PROGRESS, ">>$PROGRESS_FILE") || die "Cannot open progress file $PROGRESS_FILE: $!";
	$PROGRESS_FH = \*PROGRESS;
	print $PROGRESS_FH "Started  at ".get_today()." on ".get_hostname()."\n";

	####################################################
	# I should make this object oriented at some point!
	####################################################

	if ($SEEDED_MODE != 0) {

		############################################
		# Run seed alignment program and
		# generate seed alignment for Pairagon
		############################################

		check_seed_alignment_program($SEED_PROGRAM);
		%SEEDS = generate_seed_alignment($SEED_PROGRAM, $CDNA, $GENOMIC, $OUTPUT_DIR, $NAME, $EXTRA_GENOMIC, $BLAST_DB);
		open(FH, ">$SEED_FILE") || die "Cannot open $SEED_FILE: $!";
		foreach my $cdna (@CDNA_DEFS) { # Assuming even empty alignments make it to the hash
			my $seed = $SEEDS{$cdna};
			if (defined(%$seed)) {
				process_seed_alignment($MAX_INTRON_LENGTH, $seed);
				print_seed(\*FH, $cdna, $EXTRA_GENOMIC, %$seed);
			} else {
				die "Seed alignment cannot be found for $cdna. Please check the input sequences\n";
			}
		}
		close(FH);

		if (scalar(keys %SEEDS) == 0) {
			print "# Cannot find seed alignments for $CDNA\n";
			print $PROGRESS_FH "# Cannot find seed alignments for $CDNA\n";
			print $PROGRESS_FH "Finished at ".get_today()." on ".get_hostname()."\n";
			exit(-1);
		}

		if ($opt_P) {
			print "# Seed alignment is in $SEED_FILE\n";
			print $PROGRESS_FH "Finished at ".get_today()." on ".get_hostname()."\n";
			exit(0);
		}


		###############################################
		# Get the right genomic file from BLAST report
		###############################################

		if (defined($BLAST_DB)) {
			my ($cdna)   = sort keys %SEEDS;
			my $seed_ref = $SEEDS{$cdna};
			my %seed     = %$seed_ref;
			my $genomic  = $seed{'subject'};
			$GENOMIC     = "$GENOMIC_DIR/$genomic.fa";
			my $fasta_status = get_valid_fasta_entries($GENOMIC);
			if ($fasta_status > 1) {
				die "Cannot use multiple fasta file for genomic sequence. Try using -b instead";
			} elsif ($fasta_status == -1) {
				print "# Non [A,C,G,T,N] letters in input file $GENOMIC\n";
				print "# These will be converted to N\n";
				create_clean_fasta_entries($GENOMIC, "$OUTPUT_DIR/$NAME.genomic.fa.modified");
				$GENOMIC = "$OUTPUT_DIR/$NAME.genomic.fa.modified";
			}
			print "# Best alignment found for $genomic. Using $GENOMIC\n";
		}
	}

	############################################
	# Run Pairagon
	############################################

	my $pairagon = pairagon::instance->new(
			'PROGRAM_DIR'    => $EXE_DIR,
			'HMM_FILE'       => $HMM_FILE,
			'CDNA'           => $CDNA,
			'GENOMIC'        => $GENOMIC,
			'OUTPUT_DIR'     => $OUTPUT_DIR,
			'OPTIMIZED_MODE' => $OPTIMIZED_MODE,
			'OUTPUT_MODE'    => $OUTPUT_MODE,
			'ALIGNMENT_MODE' => $ALIGNMENT_MODE,
			'SPLICE_MODE'    => $SPLICE_MODE,
			'INSTANCE_NAME'  => $NAME
			);
	if ($SEEDED_MODE == 0) {
		# $pairagon->{SEED} will not be set!
	}  elsif (-f $SEED_FILE) {
		$pairagon->{SEED} = $SEED_FILE;
	} else {
		die "Seed alignment could not be generated. Check $PROGRESS_FILE for details";
	}
	$pairagon->check_executables;
	$pairagon->run;

	print "# Alignment is in ".$pairagon->estgen_output_file."\n";

	if ($CLEANUP == 1) {
		cleanup_after_seeding($SEED_PROGRAM, $OUTPUT_DIR, $NAME, $GENOMIC);
	}

	print $PROGRESS_FH "Finished at ".get_today()." on ".get_hostname()."\n";
	close(PROGRESS);
}

############################################
############################################
##       Subroutines
############################################
############################################

############################################
#        Utility Functions
############################################

# Uses the global variable $DEBUG
sub print_debug {
	my ($level, $subroutine, $message) = @_;
	if ($DEBUG >= $level) {
		printf "%50s:\t%s\n", $subroutine, $message;
	}
}

sub get_today {
	return scalar localtime();
}

# For lack of a faster and better way, it calls UNIX hostname.
# Feel free to fix it.

sub get_hostname {
	my $string = `hostname`;
	chomp($string);
	return $string;
}

#####################################################################
# Check validity of fasta file
#####################################################################

sub get_valid_fasta_entries {
	my ($file) = @_;
	my @defs   = ();
	my $count  = 0;
	print_debug(5, "get_valid_fasta_entries", "Trying to open $file");
	open(SEQ, "<$file") || die "Cannot open $file: $!";
	my $fasta = new FAlite(\*SEQ);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		$def =~ s/ .*//g;
		$def =~ s/>//;
		$count++;
		push(@defs, $def);
		print_debug(5, "get_valid_fasta_entries", "Checking entry $count: ".$entry->def);
		if ($entry->seq =~ /[^acgtnACGTN]/) {
			return -1;
		}
	}
	if ($count == 0) {
		die "Invalid fasta file: $file";
	}
	close(SEQ);
	return ($count, @defs);
}

sub create_clean_fasta_entries {
	my ($old_file, $new_file) = @_;
	open (ORIGINAL, "<$old_file") || die "Cannot open $old_file";
	open (MODIFIED, ">$new_file") || die "Cannot open $new_file for modified cDNA";
	my $fasta = new FAlite(\*ORIGINAL);
	while (my $entry = $fasta->nextEntry) {
		my $def = $entry->def;
		my $seq = $entry->seq;
		if ($seq =~ /[^acgtnACGTN]/) {
			$seq =~ s/[^acgtnACGTN]/N/g;
		}
		print MODIFIED "$def\n$seq\n";
	}
	close(ORIGINAL);
	close(MODIFIED);
}

############################################
#  Seed alignment hash handler
############################################

sub print_seed {
	my ($file_handler, $cdna, $offset, %seeds) = (@_);
	my ($subject, $sub_len, $score, $strand, $gb_start, $gb_end, $hsp_ref, @hsps);

	$subject  = $seeds{"subject"};
	$sub_len  = $seeds{"sublen"};
	$score    = $seeds{"score"};
	$strand   = $seeds{"strand"};
	$gb_start = $seeds{"gb_start"};
	$gb_end   = $seeds{"gb_end"};
	$hsp_ref  = $seeds{"hsps"};
	@hsps     = @$hsp_ref;

	if (defined($offset)) {
		if ($gb_start - $offset > 0) {
			$gb_start -= $offset;
		} else {
			$gb_start = 1;
		}
		if ($gb_end + $offset < $sub_len) {
			$gb_end += $offset;
		} else {
			$gb_end = $sub_len;
		}
	}

	print $file_handler ">$cdna\n";
	print $file_handler "genomic_boundary_start=$gb_start genomic_boundary_end=$gb_end strand=$strand subject=$subject\n";
	print $file_handler "count=".(scalar @hsps)."\n";
	foreach my $hsp (@hsps) {
		print $file_handler "$hsp\n";
	}
}

##########################################################################################################
# Makes sure that there isn't a VERY LONG gap in the seed alignment that is too long to be an intron.
# If there is one or more such gaps, this chooses the largest contiguous set of HSPs without such a gap.
# Remember that the reference to the hash is passed here so that you can change the seed alignment here.
##########################################################################################################
sub process_seed_alignment {
	my ($intron_limit, $seeds_ref) = @_;
	my $hsp_ref  = $seeds_ref->{"hsps"};
	my @hsps     = @$hsp_ref;
	my @long_breaks = ();
	my @new_hsps = ();

	#########################################################
	# Check for long gaps, and store the index of such a gap
	#########################################################

	push(@long_breaks, 0);
	for (my $i = 1; $i <= $#hsps; $i++) {
		if ($hsps[$i]->sb - $hsps[$i-1]->se > $intron_limit) {
			push(@long_breaks, $i);
		}
	}
	push(@long_breaks, scalar(@hsps));

	if (scalar(@long_breaks) == 2) {
		# Fine, do nothing
	} else {
		##########################################################################################################
		# There is at least one long intron in here. Need to choose a section that does not include a long intron. 
		# Choose the set with the most HSPs without long gaps
		##########################################################################################################

		print STDERR "# One or more gaps longer than $intron_limit found in the seed alignment.\n";
		print STDERR "# Using the largest subset of contiguous HSPs without such a gap.\n";
		print $PROGRESS_FH "# One or more gaps longer than $intron_limit found in the seed alignment.\n";
		print $PROGRESS_FH "# Using the largest subset of contiguous HSPs without such a gap.\n";

		my $best_break = 0;
		my $best_size  = -1;
		for (my $i = 1; $i <= $#long_breaks; $i++) {
			if ($long_breaks[$i] - $long_breaks[$i-1] > $best_size) {
				$best_size = $long_breaks[$i] - $long_breaks[$i-1];
				$best_break = $i;
			}
		}

		##########################################################################################################
		#chop from $long_breaks[$best_break-1] to $long_breaks[$best_break]-1
		# Get that new set
		##########################################################################################################

		for (my $i = $long_breaks[$best_break-1]; $i < $long_breaks[$best_break]; $i++) {
			push(@new_hsps, $hsps[$i]);
		}
		$seeds_ref->{"hsps"} = \@new_hsps;
		$seeds_ref->{"gb_start"}= $new_hsps[0]->sb;
		$seeds_ref->{"gb_end"}  = $new_hsps[$#new_hsps]->se;
	}
}

############################################
#  HSP object handler
############################################

#####################################################################
# Checks if an hsp is valid. Validity is defined as:
#  >96% identity for high quality sequence 
#      or >93% for low quality sequences, hsp->length >= 100
#      or >95% for low quality sequences, hsp->length < 100
#  and <50% A and T bases in the hsp
#  and 100% if hsp->length < 30
#####################################################################

sub hsp_validity {
	my ($hsp, $cdna_length) = @_;

	# 96% for high quality sequence only. If you suspect 
	# low quality sequence, lower the threshold

	if ((!$opt_lowQ && $hsp->percent < 96) ||
	    ($hsp->length < 100 && $hsp->percent < 95) ||
	    ($hsp->length >= 100 && $hsp->percent < 93)) {
		return -1;
	}

	#Ignore polyA and polyT stretches for pin generation

	my $polyA_p = 0.5;
	my $query = $hsp->queryAlignment;
	$query =~ s/A//ig;
	if ((length $query) < (1 - $polyA_p) * (length $hsp->queryAlignment)) {
		return -1;
	}

	$query = $hsp->queryAlignment;
	$query =~ s/T//ig;
	if ((length $query) < (1 - $polyA_p) * (length $hsp->queryAlignment)) {
		return -1;
	}

	# Short exons - should not have any mismatches. If they do, skip them
	if ($hsp->length < 30 && $hsp->percent != 100) {
		return -1;
	}

	# Short terminal exons should actually be terminal
	if ($hsp->length < 20 && ($hsp->qb > 0.2*$cdna_length || $hsp->qe < 0.8*$cdna_length)) {
		return -1;
	}

	return 1;
}

#####################################################################
# If the coordinates of hsp begin position are monotonic
#####################################################################

sub compatible_hsp {
	return (sort_genomic(@_) eq 
	        sort_cdna(@_));
}

sub sort_genomic {
	my (@list) = @_;
	@list = sort {($a->sb) <=> ($b->sb)} @list;
	return @list;
}

sub sort_cdna {
	my (@list) = @_;
	@list = sort {$a->qb <=> $b->qb} @list;
	return @list;
}

############################################
#  Generic seed alignment program handler
############################################

####################################################################################################################
# To add a new seed alignment program PROGX:
# Write the following subroutines: 
#     check_progX()                         to check the presence of the executables of PROGX
#     generate_seed_alignment_with_progX()  to generate the seed alignment using PROGX
#     cleanup_after_progX()                 to clean up intermediate files (if applicable)
# and call them from the generic seed alignment handler routines.
#
# generate_seed_alignment_with_progX() must return a hash with the following structure:
#	foreach my $cdna (@names_of_cdnas_in_cdna_file) {
#		$seed{$cdna}{"subject"} = $subject;
#		$seed{$cdna}{"sublen"}  = $sub_len;
#		$seed{$cdna}{"score"}   = $score;
#		$seed{$cdna}{"strand"}  = $strand;
#		$seed{$cdna}{"gb_start"}= 16010501;
#		$seed{$cdna}{"gb_end"}  = 16040600;
#		my @hsps                = ();
#               push(@hsps, pairagon::hsp::new(16020501, 16020600) (1, 100));
#               push(@hsps, pairagon::hsp::new(16030501, 16030600) (101, 200));
#		$seed{$cdna}{"hsps"}    = \@hsps;
#	}
#
# The routine that calles generate_seed_alignment() will then process the seed alignments and create the
# seed alignment file for Pairagon.
####################################################################################################################

sub check_seed_alignment_program {
	print "# Checking seed alignment program\n";
	print $PROGRESS_FH "# Checking for seed alignment program\n";
	my ($program) = @_;
	if ($program eq "WU-BLAST") {
		check_wu_blast();
	} else {
		die "$SCRIPT_NAME cannot generate seed alignment using $program. Please see the documentation to add this functionality";
	}
}

sub generate_seed_alignment {
	print "# Generating seed alignment\n";
	print $PROGRESS_FH "# Generating seed alignment\n";
	my ($program) = @_;
	if ($program eq "WU-BLAST") {
		return generate_seed_alignment_with_wu_blast(@_);
	} else {
		die "$SCRIPT_NAME cannot generate seed alignment using $program. Please see the documentation to add this functionality";
	}
}

sub cleanup_after_seeding {
	print "# Cleaning up temporary files\n";
	print $PROGRESS_FH "# Cleaning up temporary files\n";
	my ($program) = @_;
	if ($program eq "WU-BLAST") {
		cleanup_after_wu_blast(@_);
	} else {
		die "$SCRIPT_NAME cannot generate seed alignment using $program. Please see the documentation to add this functionality";
	}
}

############################################
#  WU-BLAST program handler
############################################

sub get_wu_blast_executables {
	return ("xdformat", "blastn");
}

sub get_wu_blast_intermediate_filenames {
	my ($name) = @_;
	return ("$name.fwd.blast", "$name.rev.blast", "$name.seed");
}

sub get_wu_blast_blastdb_filenames {
	my ($name) = @_;
	return ("$name.xnd", "$name.xnt", "$name.xns", "$name.xni");
}

sub cleanup_after_wu_blast {
	my (undef, $output_dir, $instance_name, $genomic_filename) = @_;
	my @seed_files = get_wu_blast_intermediate_filenames($instance_name);
	foreach my $file (@seed_files) {
		if (-f "$output_dir/$file") {
			system("rm $output_dir/$file");
		}
	}
	if (defined($genomic_filename)) {
		my $genomic_name = basename($genomic_filename);
		my @blastdb_files = get_wu_blast_blastdb_filenames($genomic_name);
		foreach my $file (@blastdb_files) {
			if (-f "$output_dir/$file") {
				system("rm $output_dir/$file");
			}
		}
	}
}

sub check_wu_blast {
	my ($xdformat, $blastn) = get_wu_blast_executables();
	if (!system("$xdformat 1>/dev/null 2>/dev/null")) {
		die "$xdformat program in WU-BLAST package not found in path";
	}
	if (!system("$blastn 1>/dev/null 2>/dev/null")) {
		die "$blastn program in WU-BLAST package not found in path";
	}
}

# BLAST parameters are stored in variable $blast_parameters.
sub generate_seed_alignment_with_wu_blast {
	my ($program, $cdna_file, $genomic_file, $output_dir, $name, $extra_genomic, $db_name) = @_;

	my ($xdformat, $blastn) = get_wu_blast_executables();
	my $blast_db            = "";
	my $command             = "";

	my $blast_parameters;
	my %global_seed;
	my %fwd_seeds;
	my %rev_seeds;
	my $genomic_name;

	my $alignable = 0;

	# Set the blast parameters

	$blast_parameters  = "B=10000 V=100 ";
	$blast_parameters .= "-gi -cpus=1 -warnings -nogap -links -topComboN=1 -span1 ";

	###################################################################################
	# Consider using the following for sequences with long gaps or multiple hits
	# $blast_parameters .= "-gi -cpus=1 -warnings -nogap -links -topComboN=1 -lcmask ";
	###################################################################################

	if ($opt_lowQ) {
		$blast_parameters .= "M=1 N=-3"; # Low quality. If there is a mismatch I only need 3 bases to compensate for that.
	} else {
		$blast_parameters .= "M=1 N=-20 E=1e-10"; # Need high quality. If there is a mismatch I need 20 bases to compensate for that.
	}

	# Create the blastdb if necessary

	if (defined($db_name) && (chomp($db_name) ne "")) {
		$blast_db = $db_name;
		#check for blastdb consistency
	} else {
		$genomic_name = basename($genomic_file);
		$blast_db     = "$output_dir/$genomic_name";
		$command      = "$xdformat -n -o $blast_db $genomic_file";
		print $PROGRESS_FH "$command\n";
		if (system("$command 1>/dev/null 2>/dev/null") != 0) {
			die "$xdformat failed";
		}
	}

	# Run WU-BLAST

	my ($subject, $score, $hsp_ref, @hsps);
	my ($gb_start, $gb_end, $strand);
	my $blast_file;

	# Search forward strand of cDNA
	
	$blast_file = "$output_dir/$name.fwd.blast";
	$command    = "$blastn $blast_db $cdna_file $blast_parameters -top > $blast_file";
	print $PROGRESS_FH "$command\n";
	if (system("$command") != 0) {
		die "$blastn failed";
	}
	%fwd_seeds  = process_blast_output_file($blast_file, "+");

	# Search reverse strand of cDNA

	$blast_file = "$output_dir/$name.rev.blast";
	$command    = "$blastn $blast_db $cdna_file $blast_parameters -bottom > $blast_file";
	print $PROGRESS_FH "$command\n";
	if (system("$command") != 0) {
		die "$blastn failed";
	}
	%rev_seeds  = process_blast_output_file($blast_file, "-");

	# Choose between forward and reverse for each cDNA

	foreach my $cdna (sort keys %fwd_seeds) { # Assuming even empty alignments make it to the hash
		my %seed = %{$fwd_seeds{$cdna}};
		if ($rev_seeds{$cdna}{"score"} > $fwd_seeds{$cdna}{"score"}) {
			%seed = %{$rev_seeds{$cdna}};
			print_debug(5, "generate_seed_alignment_with_wu_blast", "Best seed alignment found in - strand");
		} else {
			print_debug(5, "generate_seed_alignment_with_wu_blast", "Best seed alignment found in + strand");
		}
		if (scalar(@{$seed{"hsps"}}) != 0) {
			$alignable = 1;
		}
		$cdna =~ s/ .*//;
		$global_seed{$cdna} = \%seed;
		#print_seed(\*STDOUT, $cdna, 10000, %seed);
	}
	if ($alignable != 1) {
		%global_seed = ();
	}
	return %global_seed;
}

# Get a seed alignment for each cDNA sequence in the blast output file
# Return value is a hash containing details about each seed alignment
sub process_blast_output_file {
	my ($file, $strand) = @_;
	open(BLAST, "<$file") || die "Cannot open blast output file: $!";
	my $blast = new BPlite::Multi(\*BLAST);
	my %seed;
	print_debug(5, "process_blast_output_file", "Reading $strand strand BLASTN output file");
	while (my $report = $blast->nextReport) {
		my $cdna = $report->query;
		my $queryLength = $report->queryLength;
		my $best_subject = "";
		my $best_sub_len = 0;
		my $best_score = 0;
		my @best_hsps = ();
		print_debug(5, "process_blast_output_file", "Query = $cdna");
		while (my $subject = $report->nextSbjct) {
			my @hsps = ();
			my $score = 0;
			while (my $hsp = $subject->nextHSP) {
				my $hsp_object = pairagon::hsp::new($hsp->sb, $hsp->se, $hsp->qb, $hsp->qe);
				if ($strand eq "-") {
					my $qb = $queryLength - $hsp->qb + 1;
					my $qe = $queryLength - $hsp->qe + 1;
					$hsp_object = pairagon::hsp::new($hsp->sb, $hsp->se, $qb, $qe);
				}
				if (compatible_hsp($hsp_object, @hsps) && hsp_validity($hsp, $queryLength) != -1) {
					push(@hsps, $hsp_object);
					@hsps = sort_genomic(@hsps);
					$score += $hsp->score;
				}
			}
			print_debug(5, "process_blast_output_file", "Sbjct = ".$subject->name.", Score: $score");

			# Hack to make the pipeline choose non-random chromosomes if the seed alignments are 
			# identical in chrN and chrN_random

			if ($subject->name =~ /random/) {
				$score -= 10;
			}

			if ($score > $best_score) {
				$best_score   = $score;
				$best_subject = $subject->name;
				$best_sub_len = $subject->length;
				@best_hsps    = @hsps;
			}
		}
		$best_subject =~ s/>//;
		$best_subject =~ s/ .*//;
		$seed{$cdna}{"subject"} = $best_subject;
		$seed{$cdna}{"sublen"}  = $best_sub_len;
		$seed{$cdna}{"score"}   = $best_score;
		$seed{$cdna}{"strand"}  = $strand;
		$seed{$cdna}{"hsps"}    = \@best_hsps;
		$seed{$cdna}{"gb_start"}= 1;
		$seed{$cdna}{"gb_end"}  = 0;
		print_debug(3, "process_blast_output_file", "Best Sbjct = $best_subject, Score: $best_score");
		if (@best_hsps > 0) { # NOT(Empty BLAST file, no alignment found)
			$seed{$cdna}{"gb_start"}= $best_hsps[0]->sb;
			$seed{$cdna}{"gb_end"}  = $best_hsps[$#best_hsps]->se;
		}
	}
#foreach my $s (keys %seed) {
#	print_seed(\*STDOUT, $s, 10000, %{$seed{$s}});
#}
	return %seed;
}

#################################################
#  New seed alignment program handler goes here
#################################################

####################
# Pairagon Instance
####################

package pairagon::instance;
sub new {
	my $class = shift;
	my $self     = {@_};
	bless($self, $class);
	$self->{PAIRAGON}        = $self->program_dir."/pairagon";
	$self->{PAIRAGON2ESTGEN} = $self->program_dir."/pairagon2estgen";
	return $self
}

#####################################################################
# Check for executable files
#####################################################################
sub check_executables {
	print "# Checking Pairagon executables\n";
	my $pairagon = shift;
	my $prog_dir = $pairagon->program_dir;
	my $hmm_file = $pairagon->hmm_file;
	if (! -f "$prog_dir/pairagon") {
		die "Cannot find Pairagon executable (pairagon) in $prog_dir";
	}
	if (! -f "$prog_dir/pairagon2estgen") {
		die "Cannot find Pairagon output conversion executable (pairagon2estgen) in $prog_dir";
	}
	if (! -f "$hmm_file") {
		die "Cannot find Pairagon HMM parameter file $hmm_file";
	}
}

#####################################################################
# run Pairagon
#####################################################################
sub run {
	print "# Running Pairagon\n";
	my $pairagon = shift;
	my $command_line = $pairagon->get_alignment_command_line;
	print $PROGRESS_FH "$command_line\n";
	system($command_line) == 0 || die "Could not run ".$pairagon->pairagon;
	$command_line = $pairagon->get_conversion_command_line("estgen"); # I might make GTF/GFF output some day!
	print $PROGRESS_FH "$command_line\n";
	system($command_line) == 0 || die "Could not run ".$pairagon->pairagon2estgen;
}

#####################################################################
# command line to run executable using the options of this instance
#####################################################################
sub get_alignment_command_line {
	my $pairagon       = shift;
	my $executable     = $pairagon->pairagon;
	my $parameter_file = $pairagon->hmm_file;
	my $cdna_file      = $pairagon->cdna;
	my $genomic_file   = $pairagon->genomic;
	my $output_file    = $pairagon->state_output_file;
	my $seed_option    = "";
	my $alignment_mode = "";
	my $splice_mode    = "";
	my $optimized_mode = "";
	my $output_mode    = "";
	if ($pairagon->seed) {
		$seed_option = "--seed=".$pairagon->seed;
	}
	if ($pairagon->alignment_mode) {
		$alignment_mode = "--alignment_mode=".$pairagon->alignment_mode;
	}
	if ($pairagon->splice_mode) {
		$splice_mode = "--splice_mode=".$pairagon->splice_mode;
	}
	if ($pairagon->optimized_mode) {
		$optimized_mode = $pairagon->optimized_mode;
	}
	if ($pairagon->output_mode) {
		$output_mode = $pairagon->output_mode;
	}
	return "$executable $parameter_file $cdna_file $genomic_file $seed_option $alignment_mode $splice_mode $optimized_mode $output_mode > $output_file";
}

#####################################################################
# conversion to user friendly output
#####################################################################

sub get_conversion_command_line {
	my $pairagon = shift;
	my ($style) = @_;
	if ($style eq "estgen") {
		return $pairagon->get_estgen_conversion_command_line;
	} else {
		die "Cannot convert the state sequence to $style style output.";
	}
}

#####################################################################
# conversion to est2genome style output
#####################################################################
sub get_estgen_conversion_command_line {
	my $pairagon     = shift;
	my $executable   = $pairagon->pairagon2estgen;
	my $cdna_file    = $pairagon->cdna;
	my $genomic_file = $pairagon->genomic;
	my $input_file   = $pairagon->state_output_file;
	my $output_file  = $pairagon->estgen_output_file;
	return "$executable $input_file -cdna=$cdna_file -genomic=$genomic_file > $output_file";
}

# Hardcoded

sub pairagon        {shift->{PAIRAGON}}
sub pairagon2estgen {shift->{PAIRAGON2ESTGEN}}

# Set via new() call

sub program_dir     {shift->{PROGRAM_DIR}}
sub hmm_file        {shift->{HMM_FILE}}
sub cdna            {shift->{CDNA}}
sub genomic         {shift->{GENOMIC}}
sub seed            {shift->{SEED}}
sub alignment_mode  {shift->{ALIGNMENT_MODE}}
sub splice_mode     {shift->{SPLICE_MODE}}
sub optimized_mode  {shift->{OPTIMIZED_MODE}}
sub output_mode     {shift->{OUTPUT_MODE}}
sub output_dir      {shift->{OUTPUT_DIR}}
sub instance_name   {shift->{INSTANCE_NAME}}

sub state_output_file {
	my $self = shift;
	$self->output_dir."/".$self->instance_name.".pair";
}

sub estgen_output_file {
	my $self = shift;
	$self->output_dir."/".$self->instance_name.".estgen";
}

################
# hsp
################

package pairagon::hsp;
use overload '""' => '_overload';

sub new {
	my $hsp = bless {};
	($hsp->{SB}, $hsp->{SE}, $hsp->{QB}, $hsp->{QE}) = @_;
	return $hsp;
}

sub _overload {
	my $hsp = shift;
	return "(".$hsp->sb.", ".$hsp->qb.") (".$hsp->se.", ".$hsp->qe.")";
}

sub qb              {shift->{QB}}
sub qe              {shift->{QE}}
sub sb              {shift->{SB}}
sub se              {shift->{SE}}

