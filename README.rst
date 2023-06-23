TWINSCAN/N-SCAN 3.0 DISTRIBUTION
--------------------------------

CONTENTS:

This README file is split into the following three parts:

  I. TWINSCAN/N-SCAN Documentation
 II. Twinscan specific Documentation
III. N-SCAN specific Documentation

Address questions or comments to twinscan@cse.wustl.edu.


I. TWINSCAN/N-SCAN Documentation
--------------------------------

1. Twinscan/N-SCAN Executable
==============================

Twinscan/N-SCAN 3.5 is written in C and the source code, if included, is found in 
the src/ directory.  To create the executable for linux, type 

make linux

from the main directory.  The Twinscan/N-SCAN executable will be placed in the bin/ 
directory. 

If you wish to build for a different archetecture/OS than Linux Intel, several
special Makefile targets have been made and tested for this: alpha, macosx, 
macosxg5, and sparc.  To build Pairagon, use one of the following make targets: 
pairagon-linux, pairagon-sparc, pairagon-macosx or pairagon-macosxg5.

You may wish to add the TWINSCAN environment variable to your login script.
For example, under Linux with a bash shell, type

export TWINSCAN="/home/bart/Twinscan"

If this is a binary only distribution, the Twinscan/N-SCAN executable will already
be in the bin/ directory.  To test the executable, type  

./test-executable   

in the main directory.

Twinscan/N-SCAN has been tested under Linux with gcc version 2.96.
HOWEVER THERE IS NO GUARANTEE THAT TWINSCAN/N-SCAN WILL COMPILE OR RUN WITH YOUR 
PARTICULAR COMBINATION.  

2. Twinscan/N-SCAN Parameter Files
==================================

Each Twinscan/N-SCAN zoe HMM parameter file is specific for a particular target-informant pair
such as human-mouse. The reason for this is that the Twinscan/N-SCAN parameters are
sensitive to evolutionary distance.  Development of new parameter files is an 
ongoing project. If your favorite pairs of genomes are not present, you may be 
able to use something similar and still get satisfactory results. At this time we 
are only distributing the following Twinscan/N-SCAN zoe HMM parameter files:

Twinscan parameter files:
 worm_iscan_est-5584-genes-6-28-2005.zhmm                for C. elegans with C. briggsae as informant 
                                                         (uses explicit intron length model upto 4000 bases);
                                                         includes parameters for the est mode of Twinscan 2.03.
 crypto_iscan-1208-genes-09-15-2003.zhmm                 for Cryptococcus.
 arabidopsis_iscan-7815-genes-10-04-2004.zhmm            for A. thaliana with B. oleracea as informant
                                                         (uses explicit intron length model upto 2000 bases)
 human_iscan_est-9993-genes-07-07-2005.zhmm and
 human_iscan-9993-genes-09-13-2004.zhmm                  for H. sapiens with M. musculus as informant. The latter does not
                                                         have parameters for the EST mode of Twinscan 2.03. These
                                                         are also the generic mammalian parameter files, and 
                                                         hence can be used for M. musculus with H. sapiens as 
                                                         informant, or R. norvegicus with H. sapiens as informant.

N-SCAN parameter files:
 human_nscan-10-03-2005.zhmm                             for human (hg17) as target with mouse (mm5) as informant.


3.  Version History
===================
3.0  build 20051110RB, Nov. 2005
           Merged Twinscan/N-SCAN in a single distribution. Documentation updated to reflect the merge.
2.03 build 20051004CW, Oct. 2005
           Code optimization and several bugfixes.
           Added HMM parameter files for mammalian (human/mouse/rat) gene prediction with EST evidence.
        Dec 2004
           Added HMM and blast parameter files for mammalian (human/mouse/rat) gene prediction.
           The mammalian HMM parameter file does not include the EST mode.
2.02 build 20041011CW, Oct. 2004
           Added gene prediction with EST support.
           Added the ability to use WWAMs for INITCONS and TERMCONS models.
           Fixed a bug that counted signal scores (start, stop codon and 
           splice sites) twice for the conservation sequence model.
           Added parameter file for Arabidopsis gene prediction.
2.01 build 20040819MA, Aug. 2004
           Patch release for version 2.01. August 2004.
2.01, July 2004
           Explicit intron length model was implemented. New parameter
           file for C. elegans with explicit intron length was added.
2.0 beta, Feb. 2004
           Parameter file for Cryptococcus added February 2004.
2.0 beta, Dec. 2003  
           Released version 2.0 beta.



II. TWINSCAN SPECIFIC DOCUMENTATION
-----------------------------------

This part of the file contains the following sections:

1.  Quick Start Guide
2.  Twinscan Overview
3.  Running Twinscan - Basic Instructions
4.  Known Limitations


1.  Quick Start Guide
=====================

An example script (described in detail below) for the Twinscan analysis 
pipeline is included.  To access run

bin/runTwinscan2.pl

2.  Twinscan Overiew
=====================

Twinscan finds genes in a "target" genomic sequence by simultaneously
maximizing the probability of the gene structure in the target and the
evolutionary conservation dervied from "informant" genomic sequences.

The target sequence (i.e. the sequence to be annotated) should generally be
of draft or finished quality.  The informant can range from a single sequence 
to a whole genome in any condition from raw shotgun reads to finished assembly.  
Details about how the quality of the informant database effects predictive 
accuracy can be found in Flicek, et. al. 

Information complementary to this file can be found in the following:

P. Hu and M.R. Brent.  Using Twinscan to predict gene structures in
genomic DNA sequence.  Current Protocols in Bioinformatics (in press).

If you use Twinscan in your research, please cite the following
references:

P. Flicek, E. Keibler, P. Hu, I. Korf, M.R. Brent. Leverging the mouse
genome for gene prediction in human: from whole genome shotgun reads
to a global synteny map.  Genome Research 13. 46-54.

I. Korf, P. Flicek, D. Duan, M.R. Brent. 2001.  Integrating genomic
homology into gene-structure prediction.  Bioinformatics 17. S140-S148.

In order to run Twinscan you will need the following components:

  (1) Twinscan 3.0 executable
  (2) Twinscan zoe HMM parameter file
  (3) DNA sequence
  (4) Conservation sequence
  (5) EST sequence (optional)


(1) Twinscan Executable
-----------------------

See Section I.1

(2) Twinscan Parameter File
---------------------------

See Section I.2

(3) DNA Sequence
----------------

The target sequence must be in FASTA format, must be longer that 500 bp, and should
have the repetitive elements and low complexity sequence masked. We normally do 
this with RepeatMasker (http://ftp.genome.washington.edu/cgi-bin/RepeatMasker) 
with the MaskerAid (http://sapiens.wustl.edu/maskeraid/) improvements.  Neither
program is included with this distribution.


(4) Conservation Sequence
-------------------------

Conservation sequence is a symbolic representation of the the best alignments
between the target and informant sequences. The format of the conservation
sequence file is very simple: a definition line that includes the BLAST
database name and a second line of conservation symbols (which are
just numbers). The conseq.pl script included in the distribution creates
conservation sequence when given a DNA sequence and a BLAST report of the
target vs. the informant.  We generally use WU-BLAST (http://blast.wustl.edu)
to create the BLAST report.  NCBI BLAST works with our software, but the 
BLAST parameters found in the runTwinscan2.pl example script need to be 
changed.

(5) EST Sequence
-------------------------

EST sequence is a symbolic representation of evidence from ESTs that align to
the target sequence. The format is similar to the Conservation sequence, but 
the possible values for each position are E, I, N (to represent Exon, Intron and
not known). The estseq.pl script included in the distribution creates
EST sequence when given a DNA sequence and a (set of) BLAT reports of the
the ESTs aligned to the target.

3.  Running Twinscan 2.03 - Basic instructions
===========================================

In general there are five steps required to run Twinscan.  These are all 
contained in the example script runTwinscan2.pl, which you may have to tweak 
for your particular environment.  

The runTwinscan2.pl script sets a path to find the other files included with
this distribution.  Supporting programs from other sources such as WU-BLAST
and RepeatMasker should be in your path.  

Step 1: Mask target sequence with RepeatMasker
----------------------------------------------

While this step is not required to run Twinscan, it will improve performance
by reducing false-positive predictions.  The runTwinscan2.pl script will use 
RepeatBlaster (part of the MaskerAid package) if it is detected.  If not,
RepeatMasker is used. 

Step 2: Create informant BLAST database
---------------------------------------

Several methods for creating BLAST databases exist.  We describe ours in
Flicek et. al.

Step 3: Run BLAST
-----------------

The choice of BLAST parameters is an important consideration and will
affect both the time required for the Twinscan analysis pipeline and the 
performance of the gene-prediction algorithm.  See Flicek et. al. for
the BLAST parameters we chose to annotate the human genome.  The default 
parameters used by runTwinscan2.pl are defined in parameters/Celegans.blast.param 
file which we chose to annotate C. elegans genome. This can be changed 
with the -B option to the runTwinscan2.pl script. We have also included  
files specific to Cryptococcus, Arabidopsis and human annotation in the parameters 
directory. 

Step 4: Create conservation sequence
------------------------------------

Conseq.pl is used to create the conservation sequence required for 
Twinscan.  Conseq.pl requires BPlite.pm which is included with this distribution
in the lib/ directory.

Step 5: Create EST sequence
------------------------------------

Estseq.pl is used to create the EST sequence required for 
Twinscan. 

Step 6: Run Twinscan 2.03
------------------------

Twinscan 2.03 takes a number of command-line parameters.  One parameter
file (e.g. worm_iscan-3255-genes-12-12-2003.zhmm) and two sequence 
files (the target sequence and the conservation sequence) are required.

Twinscan may be run in "Genscan-compatible" mode by skipping the "-c=<conseq>" 
option.  In this case only the zoe HMM parameters
and the target sequence are required.  

In practice, Twinscan's memory requirements are approximately linear 
with the length of the target sequence.  A rough guideline is 1 GB of memory
for 1 Mb of input sequence.  Twinscan's native output is generally able
to be read by Genscan parsers.  An included program, zoe2gff 
converts the Twinscan 2.03 output to GTF2 (http://genes.cs.wustl.edu/GTF2.html).
Both outputs are saved by runTwinscan2.pl.

The file example.output contains the output from runTwincan2.pl using the 
BLAST parameters found in the script.  


4.  Known Limitations
======================

Genscan-compatible mode does not produce predictions that are identical
Genscan predictions.  Specifically promoters are often predicted in 
different places and exons may be slightly different near very long introns.


III. N-SCAN SPECIFIC DOCUMENTATION
----------------------------------

OVERVIEW
========

N-SCAN performs gene prediction on a target genome using information from DNA 
sequence modeling and from single or multiple genome alignments to the target. A 
parameter file with human as target and mouse as informant is included with this release. 
Programs for a simplified parameter estimation method that will allow phylogenetic 
parameters to be estimated for human as target and a limited number of other mammalian 
informants are also included.

The tar file contains the following N-SCAN specific files and this README 
file: 

	refseqs_hg17.tar		             our version of UCSC-derived refseqs
	bin/*.pl and lib/*.pm		             file manipulation programs
	chr22.align.gz			             human (mouse) alignment example for chr22
	human_nscan-10-03-2005.zhmm	     human (mouse) parameters

The gene prediction process can be divided into two sections:
	1. Parameter estimation and
	2. Gene prediction.

1. PARAMETER ESTIMATION
=======================

Parameter estimation requires 4 inputs:

	(1) DNA sequence files for both target and informant 
	(2) gene annotation files for each target chromosome in GTF format 
	(3) an old parameter file from which to copy the DNA models.  
	(4) informant-alignment files for each chromosome

(1) DNA Sequence File
---------------------

The input alphabet is {A, C, G, T, N}. Upper and lower case are treated the same. DNA 
sequence files for many organisms can be downloaded from UCSC at 
http://hgdownload.cse.ucsc.edu/downloads.html. Our best gene-prediction performance is 
obtained by downloading DNA sequence files from UCSC and masking all lower-case 
sequence to N except for low-complexity and simple repeats which are translated to 
upper-case sequence.

(2) GTF Annotation
------------------

The Brent Lab's GTF annotation file for each human chromosome (with naming 
convention chr<chr#>.gtf) is included in the download (N_SCAN requires annotation of 
5'UTR features). This annotation is a subset of the hg17 RefSeqs available for download 
from UCSC at  http://hgdownload.cse.ucsc.edu/downloads.html. A description of GTF 
format is given at http://genome.cs.wustl.edu.edu/~bio/ . Click on "resources", click on 
"software" and click on "GTF (gene transfer format)." 


(3) Old Parameter File
----------------------

To estimate a new parameter file, the simple method is to copy DNA sequence model 
parameters from human and reestimate the phylogenetic tree parameters for a new set of 
informants. For other mammals, human DNA sequence model parameters are a good 
approximation and phylogenetic  parameters can be estimated for organisms with an 
available MULTIZ alignment that includes human as target and the informant genome. 
Currently, the best gene-prediction performance is generated from UCSC's MULTIZ 
alignments.

(4) Informant-Alignment File
----------------------------

The informant-alignment file consists of a FASTA header line and one line for the target 
and one line for each informant. The length of the target line and each informant line is 
the same as the length of the DNA sequence. For each character in the DNA sequence, 
there is a corresponding character in each informant sequence. The informant sequence 
alphabet is {A, C, G, T, _, .} where an informant character from the set {A, C, G, T} 
means the informant character either aligns or mismatches the corresponding target 
characters, {_} is used for informant gaps within aligned regions, and {.} is used for 
target regions for which the given informant does not align. The program 
maf_to_align.pl, for generating chromosome-alignment files from MULTIZ alignments 
in UCSC's MAF format, is included in this package.  

MULTIZ alignments in MAF form, generated by the bioinformatics group at UCSC, are 
available for download at http://hgdownload.cse.ucsc.edu/downloads.html. Information 
about and programs for generating MULTIZ alignments directly from DNA sequence are 
available from Webb Miller's lab and can be downloaded at 
http://www.bx.psu.edu/miller_lab/ . (MAF files for pair alignments can be generated 
from BLASTZ pair alignments. Check the UCSC download page for specific examples 
of BLASTZ parameter settings. BLASTZ and programs for converting from BLASTZ 
output format (LAV) to MAF format are also available for download from Webb Miller's 
lab. Pair alignments provide a greater degree of flexibility in choosing targets and 
informants compared to the available MULTIZ alignments at the cost of a slight decrease 
in performance.)

New Parameter File
------------------

Generating a new parameter file for a based on human requires chromosome-alignment 
files for each chromosome (with naming convention chr<chr#>.align) and GTF 
annotation files. In the simplified method, there are three steps to parameter estimation:

Step 1
 
Program step_1_get_ss.pl is a Perl program that uses two Perl modules, 
GTF_Parser_UTR_SU.pm and MultiRandom.pm. These modules can be included with 
the -I option on the command line. It collects counts on alignment patterns for different 
types of feature models. As command line input it requires

	-I GTF_Parser_UTR_SU.pm
 	-I MultiRandom.pm
	<chromosome-alignment directory path>
	<output directory path>
	<GTF annotation directory path>
	<number of genomes, including target>
	<target and informant names>

It produces files containing counts for the each of the feature models (i.e. Markov chain, 
WAM,...) with the extension .ss in the output directory.

Step 2

Program step_2_train_nscan.pl is a Perl program that invokes a C program that uses the 
counts produced in Step 1 and a phylogenetic tree to train the phylogenetic parameters 
using EM. As command line input it requires

            <train command (the compiled C program)>
	<phylogenetic tree>
	<model type>
	<output directory path>

N.B. In the Brent lab system a C program is submitted to our batch queue for each model 
of conservation. Currently there are 130 such separate jobs that are launched. Some code 
modification will probably be necessary to submit these jobs to your queueing system. 

Step 3

Program step_3_CNC-phylo-zhmm.pl copies one of the 5'UTR exon models to use as a 
model for a Conserved non-coding (CNC) state, collects all of the phylogenetic feature 
models and creates a new parameter file with the new phylogenetic model and copies the 
other models from an old N-SCAN parameter file. As command line input it requires

	<output directory path>:
   	<old parameter file to use as template>
	<new parameter file name>
	<directory where downloaded programs are stored>

A parameter file with the new parameter file name will be created and placed in the 
output directory.


2. GENE PREDICTION
==================

Running N-SCAN requires a DNA sequence fragment, an informant-alignment fragment, 
and a parameter file as input for gene prediction. If the DNA sequence is 1 Mbp, then N-
SCAN requires slightly less than 2000 Mb of memory to run. This limits the length of 
DNA fragments on which gene prediction can be run.

DNA Sequence Fragment

The DNA sequence fragment on which gene prediction will be run must be in a FASTA 
file. 

Informant-alignment Fragment
----------------------------

The informant-alignment file consists of a FASTA header line and one line for the target 
and one line for each informant. The length of the target line and each informant line is 
the same as the length of the DNA sequence. The informant-alignment fragment consists 
of a FASTA header line and one line for each informant (note that the DNA sequence is 
present in the informant-alignment file, but not the informant-alignment fragment). The 
length of each informant line in the informant-alignment fragment is equal to the length 
of the DNA sequence fragment to which it corresponds.   

The program select_align.pl can be run to select the target and informants from an 
alignment file if fewer informants are desired than are present in the MULTIZ-based 
alignment file.

Parameter file
--------------

The parameter file contains the initial-state, state-transition, length-distribution, DNA-
model, and phylogenetic-model probabilities converted to log-likelihood scores. A 
parameter file (see Section I.2) for human (hg17) as target with mouse (mm5) as informant is included. An 
example informant-alignment file for human chromosome 22 with mouse as informant is 
included with this download.  

Running N-SCAN
--------------

See Section I.1. The -c and -a options are mutually exclusive. Supplying a 
informant-alignment fragment with the -a option will invoke N-SCAN. 

FILES
=====

Perl programs

step_1_get_ss.pl
step_2_train_nscan.pl
step_3_CNC-phylo-zhmm.pl 
	(fake_CNC.pl, combine_params.pl, split_zhmm.pl, combine_zhmm.pl)

maf_to_align.pl
select_align.pl

GTF_Parser_UTR_SU.pm
MultiRandom.pm

C programs
==========

iscan
train

Data sets
=========

refseqs_hg17
human (hg17) parameter file with mouse (mm5) as informant
human (hg17) alignment for chromosome 22 with mouse (mm5) as informant
