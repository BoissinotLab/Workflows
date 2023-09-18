#!/usr/bin/perl

###############################################################################################
###############################################################################################
#	This script reads a list of sets of orthologs from a simple comma separated value (CSV)
#	file, loads each protein sequence from its corresponding input file, and outputs all 
#	orthologs to a common fasta file. It does this for every line in the CSV file, i.e. each
#	line in the CSV is a group of orthologs, each ortholog in turn separated by a comma.
#
#	Orthologs do not need to be in order in the CSV, but note that whatever gene ID is in the 
#	first column will be used to name the orthologs file.
#	
#	Once the orthologs file has been created, the script can use CLUSTAL-Omega, PRANK, or MAFFT to
#	produce alignments, if the user specifies these options (and provided the packages are
#	installed in the user's server to run via command line)
#
#	INPUT
#	(+)A CSV list of groups of orthologs, where line corresponds to an orthologous group,
#		and each ortholog within each group (i.e. line) is separated by commas.
#	(+)Corresponding DNA and Protein sequences for each ortholog, in FASTA format.
#		Note that these may be generated using our script <process_gff_cds_proteins.pl>.
#	(!)This script works on the assumption that you have not yet generated the groups of orthologs,
#		but you have the list. i.e., you only have single sequences and a simple list of orthologs.
#		If you already have files with the multiple FASTA seqeuences for each orthologous group,
#		you will need to modify this script to start directly at the alignment step.
#
#	OUTPUT
#	Up to four files will be generated for each desired ID:
#	(+)<gene_ID>.orthologs.dna.fasta	| A multiple-sequence FASTA file of the Nucleotide orthologs. Not aligned.
#	(+)<gene_ID>.orthologs.prot.fasta	| A multiple-sequence FASTA file of the protein orthologs. Not aligned.
#	(+)<gene_ID>.prot_aln.fasta			| An alignment file of the protein orthologs, FASTA format.
#	(+)<gene_ID>.pal2nal.paml			| A codon-based nucleotide ortholog alignment, for use with PAML
#
#   REQUIREMENTS:
#   (+) BioPerl (www.bioperl.org)
#   For multiple-sequence alignments, you need one of the following
#   (+) CLUSTAL-Omega
#   (+) PRANK
#   (+) MAFFT
#	
#	For USAGE please run without any options or specifying flag -h or --help
#
#   ########################################################################################
#	# DISCLAIMER                                                                           #
#	# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING,#
#   # BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A      #
#   # PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE  #
#   # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL    #
#   # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF #
#   # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF #
#   # LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
#   # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF    #
#   # THE POSSIBILITY OF SUCH DAMAGE.                                                      #
#   ########################################################################################
#
#	Please send bug reports or suggestions for improvements to:
#   v.sojo.11@ucl.ac.uk
#      or
#   d.jeffares@ucl.ac.uk (or dj@katipo.org)
#
#	CREATED BY
#	Victor Sojo, Daniel Jeffares, Bartlomiej Tomiczek, Mario dos Reis
#	(C) University College London. 2014 (C)
#
#############################################################################################
#############################################################################################

# LIBRARIES
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Copy;
use Bio::SeqIO;
use Bio::Seq;

# CONSTANTS
use constant { TRUE => 1, FALSE => 0 };
my $script_function = "\n*********\nThis script aligns lists of orthologs using CLUSTAL-Omega, PRANK, or MAFFT, 
and produces codon-based nucleotide alignments using pal2nal. 
If DNA and Protein alignments are already available, the initial alignment step can be skipped 
by not specifying any alignment method.\n*********\n";

my $usage = "\n$script_function".
			"\n*********\n* Usage *\n*********".
			"\nperl $0".
			"\n\t -l (--list)       list_orthologs_file (CSV of orthologs, separated by commas, a group of orthologs in each line)".
			"\n\t[-i (--indir)      input_seqs_folder (the folder where input FASTA files are to be found)]".
			"\n\t[-o (--outdir)     output_folder (the folder where output files are to be created)]".
			"\n\t[-j (--namecolumn) column_to_name_file (default: 1. i.e. the first orthologue will be used for the filename)]".
			"\n\t[-f (--forcenums)  (use numbers instead of using a gene ID to name the orghologs file, e.g. \"OG_00001.prot_aln.fasta\")]".
			"\n\t[-c (--doclustal)  (align using CLUSTAL-Omega. This is the default alignment method)]".
			"\n\t[-C (--clustalbin) clustal_bin (If your CLUSTAL binary is not the default 'clustalo', specify it this way)]".
			"\n\t[-p (--doprank)    (align using PRANK)]".
			"\n\t[-P (--prankbin)   prank_bin (If your PRANK binary is not the default 'prank', specify it this way)]".
			"\n\t[-m (--domafft)    (align using MAFFT)]".
			"\n\t[-M (--mafftbin)   mafft_bin (If your MAFFT binary is not the default 'mafft', specify it this way)]".
			"\n\t[-A (--alnoptions) \"aln extra options\" (default \"-once -quiet\" for PRANK and \"--quiet\" for MAFFT)]".
			"\n\t[-n (--dopal2nal)  (calculate pal2nal codon alignments)]".
			"\n\t[-N (--pal2nalbin) pal2nal_bin (If your PAL2NAL binary is not the default 'pal2nal', specify it this way)]".
			"\n\t[-u (--pal2nalopt) \"<options>\" (pal2nal options. default: \"-output paml -nogap -nomismatch\")]".
			"\n\t[-x (--fileprefix)  file_prefix (eg.: \"MRSA_\", used only if general files are printed, like error files)]".
			"\n\t[-e (--errors)     (print error files. If -x is not used, files will get replaced every time)]".
			"\n\t[-s (--simplify)   (simplify gene names, e.g. change \"RP.345.d_85A\" to \"RP345d85A\")]".
			"\n\t[-q (--quiet)      (quiet, print errors but not warnings)]".
			"\n\t[-h (--help)       (display usage information and exit)]".
			"\nNote: you can use either the short version (e.g. '-l') or the long one (e.g. '--list') interchangeably.".
			"\n\nexample:\n\t./align_orthologs.pl -l orthologs.csv -i single_fastas/ -o orthos/ -c --dopal2nal -s --quiet".
			"\n\n";

# COMMAND-LINE INPUT
my ($help, $orthologs_file, $in_dir, $out_dir, $name_column, $forcenumbers, $do_clustal, $do_prank, $do_mafft,
	$clustal_bin, $prank_bin, $mafft_bin, $aln_options, $do_pal2nal, $pal2nal_bin, $pal2nal_opts, $file_prefix, $print_error_files, $simplify, $quiet);
Getopt::Long::Configure ("bundling"); # to allow case sensitivity in single-letter options
GetOptions(
	'list|l=s'			=> \$orthologs_file,
	'indir|i=s'			=> \$in_dir,
	'outdir|o=s'		=> \$out_dir,
	'namecolumn|j=i'	=> \$name_column,
	'forcenums|f'		=> \$forcenumbers,
	'doclustal|c'	    => \$do_clustal,
	'clustalbin|C=s'	=> \$clustal_bin,
	'doprank|p'		    => \$do_prank,
	'prankbin|P=s'		=> \$prank_bin,
	'domafft|m'		    => \$do_mafft,
	'mafftbin|M=s'		=> \$mafft_bin,
	'alnoptions|A=s'	=> \$aln_options,
	'dopal2nal|n'		=> \$do_pal2nal,
	'pal2nalbin|N=s'	=> \$pal2nal_bin,
	'pal2nalopt|u=s'	=> \$pal2nal_opts,
	'fileprefix|x=s'	=> \$file_prefix,
	'errors|e'			=> \$print_error_files,
	'simplify|s'		=> \$simplify,
	'quiet|q'			=> \$quiet,
	'help|h'			=> \$help
	);

if ($help || !(-e $orthologs_file)) { die $usage; }

# DEFAULT USER INPUT VALUES
$in_dir			||= '';
$out_dir		||= '';
$clustal_bin	||= 'clustalo';
$prank_bin		||= 'prank';
$mafft_bin		||= 'mafft';
$aln_options	||= '';
$pal2nal_bin	||= 'pal2nal.pl';
$pal2nal_opts	||= '-output paml -nogap -nomismatch';
$file_prefix	||= '';
$name_column	||= 1;
--$name_column; # because in perl the first index is zero

if ($in_dir  && !((substr $in_dir,  -1, 1) eq "/")) { $in_dir  .= "/"; }
if ($out_dir && !((substr $out_dir, -1, 1) eq "/")) { $out_dir .= "/"; }
# remove a dot at the end if it's there
if ($out_dir && ((substr $out_dir, -1, 1) eq ".")) { $out_dir = substr $out_dir, 0, -1; }
if ($in_dir  && ((substr $in_dir,  -1, 1) eq ".")) { $in_dir  = substr $in_dir,  0, -1; }

# CHECK OPTION DEPENDENCIES
my $warnings = "";
my $missing_IDs = "";

# don't let the user pick two or more alignment methods
if ($do_clustal && $do_prank || $do_clustal && $do_mafft || $do_prank && $do_mafft)
{
	die "\n********************\n*   FATAL ERROR!   *\n********************\n"
		 ."You have selected more than one alignment method. Please use only one alignment flag:\n"
		 ."\t-c for CLUSTAL-omega\n"
		 ."\t-p for PRANK\n"
		 ."\t-m for MAFFT\n"
		 ."\n";
}


if ($do_prank && $aln_options eq '') { $aln_options = "-once -quiet"; }
if ($do_mafft && $aln_options eq '') { $aln_options = "--quiet"; }

my $do_alignment = $do_clustal||$do_prank||$do_mafft;

if ($do_pal2nal && !$do_clustal && !$do_prank && !$do_mafft)
{
	# user wants to align codons with pal2nal but forgot to select either of the alignment tools
	print "\n****************\n*   WARNING!   *\n****************\n".
		"You've chosen to calculate codon alignments using PAL2NAL, but ".
		"you have not selected any alignment method.\n".
		"We assume you already have aligned files in the <output> (not input!) folder with name format:\n".
		"\t<GENE_ID>.prot_aln.fasta\n".
		"where <GENE_ID> is the ID of the ortholog in the first ID of the orthologous group in the orthologs file, or:\n".
		"\tOG_00037.prot_aln.fasta\n".
		"where in \"OG_00037\" the number 37 corresponds to the position (line number, starting at 1) in the orthologs file.\n\n";
	continue_or_die();
}
################
# TO DO: check if the output directory exists!!!! Otherwise create it
################

my ($sp_gene,$so_gene,$sc_gene,$sj_gene) = "";
my @orthologs = ();
my $dna_seqIn;
my $dna_seqIO_object;
my $dna_seqOut;
my $prot_seqIn;
my $prot_seqIO_object;
my $prot_seqOut;
my $gene_id = "";
my $missing = 0;
my $something_missing = FALSE;
my $found = 0;
my $file_name = "";
my $key;
my $percentCut = 5;	# progress will be printed every 5%
	
if ($do_clustal || $do_prank || $do_mafft)
{
	print "\n****************\n*   WARNING!   *\n****************\n".
			"You've chosen to calculate protein sequence alignments.\n".
			"This may take several HOURS!!!\n".
			"Percentage indicators will be printed every ${percentCut}\%, ".
			"but if you think it's taking too long, hit Ctrl+Z to cancel.\n".
			"(IMPORTANT NOTE: Ctrl+C might not stop the script, use Ctrl+Z!).\n";
	continue_or_die();
}

# open and load file into array of lines, each containing a set of orthologs
open FILE, $orthologs_file or die $!;
my @ortho_file_lines = <FILE>;
close FILE;
my $line;

my $n_orthogroups = @ortho_file_lines; # total number of genes minus those that weren't multiples of 3
my $i = 0;
my $percentStep = int($n_orthogroups * $percentCut / 100);
if ($n_orthogroups < 20)
{
	$percentStep = 1;
	$percentCut = int(100/$n_orthogroups);
}
my $curPctMult = 1;

my $empty_files = 0;
my $good_files = 0;
my $MAX_EMPTY_FILES = 20;

my $orthocounter = 0;
my $testCounter = 10;

# iterate the lines extracted from the orthologs file (each line is a group of orthologs separated by commas)
foreach (@ortho_file_lines)
{
####################
--$testCounter;
#if (!$testCounter) { last; }
####################
	$line = $_;
	chomp($line);
	
	# if ($line =~ /^[\w\.]+,[\w\.]+/) # line is valid CSV of orthologs: has a succession of letters/numbers, then a comma, then another succession
	if ($line !~ /^\s/ && $line =~ tr/,//) # line is not empty and there's at least one comma in it
	{
		# each lines has the orthologs seprarated by commas, so split them
		@orthologs = split(/,/, $line);
		# for each ortholog, load the protein sequence (if it exists), and put it into a common FASTA file
		# note: ortholog files will be named after the id of the first ortholog ("column") unless user chose to use numbers
		++$orthocounter;
		if ($forcenumbers)
		{
			$file_name = sprintf("OG_%05i", $orthocounter);
		}
		else
		{
			# get the gene id of the column chosen to name the output file
			$file_name = get_gene_id($orthologs[$name_column]);
		}
		# create the common DNA and Protein FASTA files that will hold all the sequences for this group of orthologs
		$dna_seqOut  = Bio::SeqIO->new(-file => ">${out_dir}${file_name}.orthologs.dna.fasta",  -format=>"Fasta", -alphabet=>'DNA');
		$prot_seqOut = Bio::SeqIO->new(-file => ">${out_dir}${file_name}.orthologs.prot.fasta", -format=>"Fasta", -alphabet=>'protein');
		$warnings = "";
		$found = 0;
		$something_missing = FALSE;
		
		foreach (@orthologs)
		{
			$gene_id = get_gene_id($_);
			# check if both the DNA and Protein files exist, otherwise something is missing
			if ($_ && (-e "${in_dir}${gene_id}.dna.fasta") && (-e "${in_dir}${gene_id}.prot.fasta"))
			{
				# read the DNA file for this ortholog
				$dna_seqIO_object = Bio::SeqIO->new(-file => "${in_dir}${gene_id}.dna.fasta", -format=>"fasta", -alphabet=>"DNA");
				$dna_seqIn   = $dna_seqIO_object->next_seq();
				# write it to the common DNA file for the orthologous group
				$dna_seqOut->write_seq($dna_seqIn);
				# read the Protein file for this ortholog
				$prot_seqIO_object = Bio::SeqIO->new(-file => "${in_dir}${gene_id}.prot.fasta", -format=>"fasta", -alphabet=>"protein");
				$prot_seqIn   = $prot_seqIO_object->next_seq();
				# write it to the common Protein file for the orthologous group
				$prot_seqOut->write_seq($prot_seqIn);
				++$found;
			}
			else
			{
				$warnings .= "${gene_id},";
				my $temp = $_;
				substr($temp,4,0,"_");
				$missing_IDs .= $temp."\n";
				$something_missing = TRUE;
			}
		}
		
		# count how many genes are missing
		if (!$quiet && $warnings)
		{
				print "***WARNING*** Protein and/or DNA file(s) not found for gene(s): "
					. (substr $warnings, 0, -1) . "\nWere ignored!\n";
		}
		
		if (!(-s "${out_dir}${file_name}.orthologs.prot.fasta") || !(-s "${out_dir}${file_name}.orthologs.dna.fasta"))
		{
			# files are empty! perhaps user specified the wrong input directory. Delete both files
			unlink ("${out_dir}${file_name}.orthologs.prot.fasta");
			unlink ("${out_dir}${file_name}.orthologs.dna.fasta");
			++$empty_files;
		}
		else { ++$good_files; }
		
		
		if (!$good_files && $empty_files >= $MAX_EMPTY_FILES)
		{
			my $die_message = "\n******************** FATAL ERROR ********************\n" .
				"Code aborted because the first $MAX_EMPTY_FILES lists of orthologs that ".
				"we attempted to process had all members missing and until ".
				"that point not a single list had been processed successfully.\n";
			if ($in_dir) {$die_message .= "It's very likely that you haven't specified the right input folder!";}
			else {$die_message .= "Did you forget to specify an input folder? Use the -i option.";}	
			$die_message .= "\n*****************************************************\n\n";
			die $die_message;
		}
		# align protein sequences if user wants to
		if ($found >1) # no sense in aligning only one sequence
		{
			if ($do_clustal)
			{
				# align the proteins using clustal omega via command line
				system("${clustal_bin} ${aln_options} -i ${out_dir}${file_name}.orthologs.prot.fasta".
						" -o ${out_dir}${file_name}.prot_aln.fasta");
			}
			elsif ($do_prank)
			{
				system("${prank_bin} -d=\"${out_dir}${file_name}.orthologs.prot.fasta\""
						." -o=\"${out_dir}${file_name}.prot_aln.prank\" ${aln_options}");
				move("${out_dir}${file_name}.prot_aln.prank.1.fas", "${out_dir}${file_name}.prot_aln.fasta");
			}
			elsif ($do_mafft)
			{
				system("${mafft_bin} ${aln_options} ${out_dir}${file_name}.orthologs.prot.fasta"
							."> ${out_dir}${file_name}.prot_aln.fasta");
			}
		}
		if ($something_missing) { ++$missing; }
		
		if ((!$do_alignment && $do_pal2nal) || ($found>1 && $do_alignment && $do_pal2nal))
		{
			# do pal2nal if requested by the user
			system("${pal2nal_bin} ${out_dir}${file_name}.prot_aln.fasta ".
					"${out_dir}${file_name}.orthologs.dna.fasta ${pal2nal_opts} ".
					">${out_dir}${file_name}.pal2nal.paml");
		}
	} # end if: the line matches a regex to find a comma-separated list of orthologs
	else # the line did not match the regex
	{
		print "***ERROR*** this line could not be understood! Orthologs must be alphanumeric (may include underscores or dots), and separated by commas.\n";
		print "\t".$line;
	}
	
	++$i;
	# print the percentage of files processed
	if (!($i % $percentStep) && $percentCut*$curPctMult<100)
	{
		print $percentCut*$curPctMult . "\%: ${i} of ${n_orthogroups} lists of orthologs have been processed...\n";
		++$curPctMult;
	}
	if ($i == $n_orthogroups)
	{
		print "100\% of ${n_orthogroups} lists of orthologs processed!\n";
	}
}

if ($missing)
{
	print "\n***WARNING*** ${missing} of the ${n_orthogroups} groups of orthologs specified in file <${orthologs_file}>".
			" had at least one missing member, i.e.: the protein sequences file could not be found!!!\n\n";

	if ($print_error_files)
	{
		open(MISSINGIDSFILE, ">${out_dir}${file_prefix}.missing_IDs.txt"); # will overwrite!
		#print MISSINGIDSFILE "Timestamp: "; #write text, no newline
		#print MISSINGIDSFILE &timestamp(); #write text-returning fcn
		#print MISSINGIDSFILE "\n"; #write newline
		#*** Print freeform text, semicol required ***
		print MISSINGIDSFILE $missing_IDs;
		close(MISSINGIDSFILE);
	}
}

sub continue_or_die
{
	# DESC: prints a request for the user to hit ESCAPE/q to stop the run,
	# or any other key to continue
	print "\nPress 'q' or ESCAPE to quit now or any other key to continue...";
	use Term::ReadKey;
	ReadMode('cbreak');
	$key = ReadKey(0);
	ReadMode('normal');
	if ((lc($key) eq 'q') || ord($key) == 27 ) { die "\n\nYou've cancelled this run!\n\n"; }
	
	print "\n\n";
}

sub get_gene_id
{
	# DESC: extracts and cleans the gene ID from a dirty entry
	# initially, the ID can be something such as: "ID=SPOG_05762;Name=NAD dependent epimerase/dehydratase"
	# so need to clean it
	my $uglyID = $_[0];
	$uglyID =~ s/ID.//;	# remove "ID=" or "ID " or "ID#"
	$uglyID =~ s/ .*//;	# remove everything from the first space to the end
	if ($simplify)
	{	
		$uglyID =~ s/\W//g; # remove any non-word characters
		$uglyID =~ s/_//g;	# remove any underscores (don't get removed by \W)
	}
	return $uglyID; 	# no longer ugly
}
