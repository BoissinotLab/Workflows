#!/usr/bin/perl

###############################################################################################
###############################################################################################
#	This script reads annotation information from a GFF 3 file, and extracts the ID,
#	coding sequence (CDS) positions, supercontig, and strand (+/-) for each gene.
#	It then uses that information to read the corresponding DNA from a supercontigs FASTA file
#	and the protein sequence from another FASTA file, first checking whether the ID is in a
#	list that may be specified in a separate file. If no list is given, all IDs are processed
#	
#	The script then assembles together all CDS sequences in the right order and outputs two files
#	for each gene: the protein FASTA sequence, and the CDS(s).
#
#	Any genes that could not be processed (e.g., CDS present but no protein, or vice-versa)
#	will be detailed in the log file and printed out to a file.
#	This script translates the full CDS and checks for equality with the protein,
#   so another typical reason for failure is a mismatch between the translated CDS and the protein sequence
#
#	INPUT
#	(+)A GFF3 file with annotations specifying, for each desired protein:
#		(.)The Gene's unique identifier (ID)
#		(.)The DNA block (contig, supercontig, etc.) that contains the coding sequences
#		(.)The start and end positions of all CDS fragments of a protein-coding gene
#		(.)The strand in which the sequence is (+ or -). If (-), need to get reverse complement
#	(+)A file containing all the supercontigs with the corresponding DNA in FASTA format
#	(+)A file containing all the protein sequences, also in FASTA format
#	( )A file specifying which IDs are desired. This will be a simple list with each ID in a 
#	   new line. If no file is specified, all gene IDs will be processed wherever possible
#
#	OUTPUT
#	A pair of files for each desired ID:
#	(+)<gene_ID>.dna.fasta		
#	(+)<gene_ID>.prot.fasta		
#	Additionally, a log file:
#	(+)<dateTime>.log			contains information on how the script ran and details of any errors
#
#	REQUIREMENTS
#	(+)BioPerl (www.bioperl.org)
#
#	USAGE AND INSTRUCTIONS
#   For usage and instructions please run without any options or specifying flag -h or --help
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
#	CREATED BY
#	Victor Sojo, Daniel Jeffares, Bartlomiej Tomiczek, Mario dos Reis
#	University College London (2014)
#
#   Please send bug reports or suggestions for improvements to:
#   v.sojo.11@ucl.ac.uk
#      or
#   d.jeffares@ucl.ac.uk (or dj@katipo.org)
#
#############################################################################################
#############################################################################################

=pod
	0)	Check that the three mandatory files have been provided, and whether the optional values were specified.
	1)	Load all the DNA blocks (i.e. supercontigs) into a hash of DNA sequence objects: |key=dnaBlockID|DNAseq|
	2)	Load all desired IDs from the optional -l file into an array. If no -l/-i, define a boolean to read
		all. Alternatively, user may have specified unwanted IDs via -i. In that case, define another
		corresponding boolean called ignore=TRUE
	3)	Run through the GFF file.
		Treat first sequence as a special case. Afterwards, for every 'gene' found, output the previous sequence,
		check if ID is in the wanted list, empty/re-declare a seq object, and gather every CDS from the
		appropriate supercont, taking care of reversing and pasting in the right order if strand is (-).
		Once another GENE is found put the necessary info into the following hash:
		|{key=ID}|dnaSeq|
		Dump original DNA data to release memory
	4)	Load the proteins file. Traverse sequence by sequence; if ID is in selected list/in the hash, translate
		the DNA sequence, check for equality, and if all good output both DNA and protein as fasta files.
=cut

# LIBRARIES & MODULES
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::Seq;

# CONSTANTS
use constant { TRUE => 1, FALSE => 0 };

# text to display usage information
my $script_function = "\n*********\nThis script reads gene annotation information from a GFF file, 
extracts the corresponding DNA sequences from a contig/genome FASTA file, 
and outputs FASTA sequences of the coding sequence (DNA) and corresponding proteins for each gene, 
potentially validating against provided amino acid sequences.\n*********\n";

my $usage = "\n$script_function".
			"\n*********\n* Usage *\n*********\n".
			"perl $0".
			"\n\t -g (--gff)        gff_file".
			"\n\t -d (--dna)        dna_file (FASTA format by default)".
			"\n\t[-p (--proteins)   proteins_file (FASTA format by default)]".
			"\n\t[-G (--gffversion) gff_version (default: 3)]".
			"\n\t[-l (--list)       file_of_ids (If only some genes in the GFF need to be processed, otherwise all will be)]".
			"\n\t[-f (--format)     seq_format (default: \"FASTA\")]".
			"\n\t[-o (--outdir)     out_dir]".
			"\n\t[-x (--fileprefix) file_prefix (eg.: \"Sp\". Used only if general files are printed, like error files)]".
			"\n\t[-e (--errors)     (print error files. If -x is not used, files will get replaced every time)".
			"\n\t[-s (--simplify)   (simplify gene names, e.g. change \"RP.345.d_85A\" to \"RP345d85A\")]".
			"\n\t[-h (--help)       (print usage information and exit)".
			"\nNote: you can use either the short version (e.g. '-g') or the long one (e.g. '--gff') interchangeably.".
			"\n\nexample:\n\t./process_gff_cds_proteins.pl -g input/Rp_genes.gff3 -d input/Rp_DNA.fasta -p input/Rp_proteins.fasta -o single_fastas/".
			"\n\n";

# COMMAND-LINE INPUT
my ($gff_file,$gff_version,$dna_file,$prot_file,$seq_format,$list_file,
	$out_dir,$file_prefix,$print_error_files,$simplify,$help);
Getopt::Long::Configure ("bundling"); # to allow case sensitivity in single-letter options
GetOptions(	
	'gff|g=s' 		=> \$gff_file,
	'gffversion|G=i'=> \$gff_version,
	'dna|d=s'  		=> \$dna_file,
	'proteins|p=s'	=> \$prot_file,
	'format|f=s'	=> \$seq_format,
	'list|l=s'		=> \$list_file,
	'outdir|o=s'	=> \$out_dir,
	'fileprefix|x=s'=> \$file_prefix,
	'errors|e'		=> \$print_error_files,
	'simplify|s'	=> \$simplify,
	'help|h'		=> \$help
	);

# check if user only wants to display the usage information
if ($help) { die $usage; }

##############################################
# check whether necessary files have been specified, else print usage information
##############################################
(($gff_file && -e $gff_file) && ($dna_file && -e $dna_file)) || die $usage;

# DEFAULT USER INPUT VALUES
$gff_version	||= 3;
$seq_format		||= 'fasta';
$out_dir		||= '';
$file_prefix	||= '';
$simplify		||= TRUE;

if ($out_dir && ((substr $out_dir, -1, 1) eq ".")) { $out_dir = substr $out_dir, 0, -1; }
if ($out_dir && !((substr $out_dir, -1, 1) eq "/")) { $out_dir .= "/"; }


##############################################
# load list of sequences from the dna_file
##############################################
print "\nLoading DNA blocks...\n";
# make a hash of DNA seq objects %dna_blocks
my %dna_blocks = ();
# and another hash %dna_genes
my %dna_genes = ();
my $check_IDs = FALSE;
my %id_list = ();
my $gene_count = 0;
		
my $dna_inObjs = Bio::SeqIO->new('-format'=>"$seq_format", -file=>"$dna_file", -alphabet=>'dna');
while (my $seq = $dna_inObjs->next_seq())
{
	$dna_blocks{$seq->id} = $seq;
}
undef $dna_inObjs;
print "done!\n";

##############################################
# process the list of desired IDs, if any
##############################################
if ($list_file)
{
	if(-e $list_file)
	{
		print "\nLoading list of IDs to process...\n";
		# user has specified a list of desired IDs
		# set boolean to check the list later in the code
		$check_IDs = TRUE;
		
		# load the file into the %id_list hash
		open FILE, $list_file or die $!;
		my $key;
		while (my $line = <FILE>)
		{
			chomp($line);
			if ($line !~ /^\s/)
			{
				#($key) = $line =~ /^\S+/g; # this trims (removes whitespace from both ends), but is actually unnecessary
				$key = get_gene_id($line); #clean the gene ID
				$id_list{$key} = 1; # that 1 is irrelevant, it's just to have something there. Searching hashes is much quicker than arrays in Perl
			}
		}
		close FILE;
		print "done!\n";
	}
	else { print_file_error($list_file); die "\nNote that this file is not mandatory; if you don't specify a list of IDs all the genes will be processed.\n" }
}

##############################################
# Load and process GFF FILE annotations
##############################################
print "\nLoading GFF annotations and extracting DNA gene fragments...\n";

# ??????????????????????????
# TODO: THIS CODE NEEDS TO BE IMPROVED!!! IT'S CURRENTLY SPLITTING EACH GFF LINE
# INTO AN ARRAY AS [VS] COULDN'T FIND A WAY TO MAKE Bio::Tools::GFF WORK!!!
# ??????????????????????????
my ($contig,$tag,$start,$end,$strand,$id);
if ($gff_version==3)
{
	# GFF3 example:
	# ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
	($contig,$tag,$start,$end,$strand,$id)=(0,2,3,4,6,8);
}
elsif ($gff_version==2)
{
	# GFF2 example:
	# Chr1  curated  CDS 365647  365963  .  +  1  Transcript "R119.7"
	($contig,$tag,$start,$end,$strand,$id)=(0,2,3,4,6,8);
}

my $current_gene = '';
my $previous_gene = '';
my $seq;
my $gff_obj = Bio::Tools::GFF->new(-file => "$gff_file", -gff_version => "$gff_version");

# TODO: change the code to receive GFFs via standard input
# my $gff_obj = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => "$gff_version");
# run the script like this: "cat input.gff | perl myscript.pl "

my $ignore = FALSE; # this variable will be used to ignore IDs that aren't on the list

my $n_block_errors=0;
use constant { MAX_BLOCK_ERRORS => 10 };
my $n_IDs = 0;
if ($list_file) { $n_IDs = keys %id_list;}
my $non_mult_3_errors = 0;
# a hash to hold all non-multiple-of-3 sequence IDs
my %non_mult_3_IDs;

# 
while (my $gffFeature = $gff_obj->next_feature())
{
	my @gffLine = split(/\t/, $gffFeature->gff_string());
	# TODO: ??????????????????????????
	# if ($gffFeature->has_tag('Gene')) DOES NOT WORK!
	# ??????????????????????????
	my $current_tag = lc($gffLine[$tag]);
	if ($current_tag eq "gene" || $current_tag eq 'locus_tag')
	{
		# It's a new gene, so all CDSs for the previous one have been concatenated together
	
		# Increase the counter
		$gene_count++;
	
		# Output previous fasta seq to the genes hash, but first check if it is a multiple of 3
		$non_mult_3_errors += create_dna_seq($seq, $current_gene);
		
		$ignore = TRUE;
		# empty the seq to start over
		$seq = '';
		
		# ??????????????????????????
		# TODO: ($current_gene) = $gffFeature->get_tag_values('gene'); DOES NOT WORK!!!
		# ??????????????????????????
		$current_gene = get_gene_id($gffLine[$id]);
		
		if (!$check_IDs)
		{
			$ignore = FALSE; # we're not checking IDs, so don't ignore anything
		}
		elsif (exists $id_list{$current_gene})
		{
			$ignore = FALSE; # id is in the list, we don't want to ignore it 
			# remove the gene ID from the list
			delete $id_list{$current_gene};
		}
	}
	elsif ($current_tag eq 'cds' && !$ignore)
	{
		if (exists $dna_blocks{$gffLine[$contig]})
		{
			#my $subseq = $dna_blocks{$gff->seq_id}->subseq($gff->start, $gff->end); DOES NOT WORK!!!
			my $subseq = $dna_blocks{$gffLine[$contig]}->subseq($gffLine[$start], $gffLine[$end]);
			
			if ($gffFeature->strand < 0) # SOMEHOW THIS WORKS, BUT OTHER GFF METHODS DON'T
			{
				# gene is in negative strand, need to calculate the reverse complement of the sequence
				my $subseqObj = Bio::Seq->new(-seq => "$subseq", id => "irrelevant");
				$subseq = $subseqObj->revcom->seq;
			}
			$seq .= $subseq;
		}
		else #  DNA block does not exist... this might be a nomenclature error
		{
			++$n_block_errors;
			if ($n_block_errors >= MAX_BLOCK_ERRORS)
			{
				die "\n******************** FATAL ERROR ********************\n" .
					"Code aborted after failing to find ".MAX_BLOCK_ERRORS." DNA blocks as specified in the GFF file!\n" .
					"It is likely that there are discrepancies between the names of your DNA blocks in the GFF annotations and those in the DNA file\n". 
					"e.g., they are called \"supercontig\" in one and \"supercont\" or \"Chromosome\" in the other.\n".
					"You may have to make this change with a Find+Replace operation in a text editor." .
					"\n*****************************************************\n\n";
			}
		}
	}
	# else: ignore any other tags
}

# there might still be a sequence floating around, because they only get printed upon finding a new gene (which
# might not happen in the last case). Print it if so, just like the others
$non_mult_3_errors += create_dna_seq($seq, $current_gene);

# we're done with the %dna_blocks, dump them, since any useful info is now in the %dna_genes hash
undef %dna_blocks;



###################################################
# Print info about any DNA- seq processing errors
###################################################

# check that at least one of the IDs has been processed
if ($list_file && $n_IDs == keys %id_list)
{
	die "\n******************** FATAL ERROR ********************\n" .
		"Code aborted because none of the $n_IDs genes specified in IDs file <$list_file> " .
		"were found in the DNA blocks file <$dna_file>!!!\n" .
		"It is likely that you've specified a wrong file name." .
		"\n*****************************************************\n\n";
}

# any gene still in the id_list was not removed above, so wasn't processed
foreach (keys %id_list)
{
	print "***WARNING*** Gene $_ was not found in the GFF annotations, could not be processed.\n";
}

# the following three blocks print a summary of successes and errors
# how many genes were processed successfully
if ($list_file) { printf("\ndone!: %.0f%% (%d/%d) of the genes specified in IDs file <$list_file> were processed successfully!\n",(($n_IDs-(keys %id_list)-$non_mult_3_errors)*100/$n_IDs), ($n_IDs-(keys %id_list)-$non_mult_3_errors), $n_IDs); }
else { printf("\ndone!: %.1f%% (%d/%d) of the genes in the file <$gff_file> were processed successfully!\n",(($gene_count-(keys %id_list)-$non_mult_3_errors)*100/$gene_count), ($gene_count-(keys %id_list)-$non_mult_3_errors), $gene_count); }

if ($non_mult_3_errors || (keys (%id_list)))
{
	print "However, please note that ";
}
# how many had lengths that were not multiples of 3
if ($non_mult_3_errors>1)
{
	print "$non_mult_3_errors sequences had lengths that weren't multiples of 3, thus were ignored.\n"
}elsif ($non_mult_3_errors==1){print "1 sequence length was not a multiple of 3, so it was ignored.\n"}
# how many were not found in the GFF annotations
if (keys(%id_list) > 1)
{
	print "". (keys %id_list) . " genes were not found in the GFF annotations, so could not be processed.\n";
}elsif (keys(%id_list)==1){print "1 gene was not found in the GFF annotations, so could not be processed.\n"}

# Print files with the lists of missing & non-mult-3 IDs, if user wants to
if ($print_error_files)
{
	my $id = "";
	
	#missing annotations:
	my $missing_annotations = "";
	foreach (keys (%id_list))
	{
		$id = $_;
		substr($id,4,0,"_");
		$missing_annotations .= "$id\n";
	}
	if ($missing_annotations)
	{
		open(MISSINGIDSFILE, ">$out_dir"."$file_prefix"."_missing_annotations.txt"); # will overwrite!
		print MISSINGIDSFILE $missing_annotations;
		close(MISSINGIDSFILE);
	}
	undef $missing_annotations;
	# sequences with non-multiple-of-3 number of bases:
	my $non_mult_3_list = "";
	foreach (keys (%non_mult_3_IDs))
	{
		$id = $_;
		substr($id,4,0,"_");
		$non_mult_3_list .= "$id\n";
	}
	if ($non_mult_3_list)
	{
		open(MISSINGIDSFILE, ">$out_dir"."$file_prefix"."_non_mult_3_IDs.txt"); # will overwrite!
		print MISSINGIDSFILE $non_mult_3_list;
		close(MISSINGIDSFILE);
	}
	undef $non_mult_3_list;
	undef %non_mult_3_IDs;
}


###################################################
# Translate sequences and output DNA & Prot files
###################################################
my $n_genes = $gene_count-(keys %id_list)-$non_mult_3_errors;# total number of genes minus those that weren't multiples of 3
my $i = 1;
my $tenPercentStep = int($n_genes / 10);
my $curPctMult = 1;

# NOW LOAD PROCESS THE PROTEINS FILE, IF SPECIFIED
# load list of sequences from the $prot_file if any
my $gene_id = "";
if ($prot_file)
{
	if (-e $prot_file)
	{
		print "\nCalculating translations, writing output files, and checking whether DNA translations are identical to sequences in the proteins file...\n";
		my $unmatched = 0; # how many CDS sequences do not match the provided aminoacid sequence
		
		my $prot_inObjs = Bio::SeqIO->new('-format'=>"$seq_format", -file=>"$prot_file", -alphabet=>'protein');
		while (my $prot_seqObj = $prot_inObjs->next_seq())
		{
			$gene_id = get_gene_id($prot_seqObj->id);
			
			# WARNING!!!: this is not a universal fix, it is specific to this problem because proteins have ID SJAG_06636T0
			# whereas genes have ID SJAG_06636 (without the T0)
			if (substr $gene_id, -2, 2 ) { $gene_id = substr $gene_id, 0, -2; }
			
			# check if the dna block is in the list of desired genes with a valid sequence
			if (exists $dna_genes{$gene_id})
			{
				if ($dna_genes{$gene_id}->translate->seq() ne $prot_seqObj->seq())
				{
					print "***WARNING*** Protein and CDS sequences for gene $gene_id do not match. DNA translation will be used!\n";
					++$unmatched;
				}
				# Output both to respective files
				output_dna_prot_files($dna_genes{$gene_id});
			}
			delete($dna_genes{$gene_id}); # to avoid repetitions
		}
		print "\ndone! Finished processing protein translations.\n\n";
		if ($unmatched == 1)
		{
			print("However, please note that 1 protein amino-acid sequence did not match the translation of the CDS, so the latter was used.\n\n");
		}
		elsif ($unmatched > 1)
		{
			print("However, please note that $unmatched protein amino-acid sequences did not match the translation of the coresponding CDSs, so the latter were used.\n\n");
		}
	}
	else
	{
		print_file_error($prot_file);
		die "\nPlease note that you only need to specify a proteins file if you want to compare " .
			  "DNA tranlations to otherwise-obtained protein sequences\n";
	}
}
else # user does not want to check the protein sequences or does not have them. Just translate the DNA.
{
	print "\nCalculating protein translations and writing output files...\n";
	my $key;
	
	# TO DO: ADD ERROR FLAG DETECTOR AND CHECK FOR FINAL OUTPUT ACCORDINGLY
	# AND OUTPUT ALL THE ERRORS TO A LOG FILE
	
	foreach $key(keys %dna_genes)
	{
		output_dna_prot_files($dna_genes{$key});
	}
	print "\ndone! Finished processing protein translations.\n\n";
}

##################################
#          SUB-ROUTINES          #
##################################
sub create_dna_seq
{
	# DESC: takes a DNA text seq and a gene_id, and creates a BioPerl Seq object at position gene_id of
	# the dna_genes hash
	my $seq = $_[0];
	my $gene_id = $_[1];
	if ($seq && $gene_id) # check that there's something in there
	{
		if (!(length($seq) % 3))
		{
			$dna_genes{$gene_id} = Bio::Seq->new(-seq =>"$seq", id=>"$current_gene");
		}
		else
		{
			print "***ERROR*** Number of bases in CDS regions for <${current_gene}> is not a multiple of 3! Will be ignored.\n";
			$non_mult_3_IDs{$current_gene}=1;
			return 1;
		}
	}
	return 0;
}

sub get_gene_id
{
	# DESC: extracts and cleans the gene ID from a messy entry
	# initially, the ID can be something as ugly as: "ID=SPOG_05762;Name=NAD dependent epimerase/dehydratase"
	# and we need just "SPOG05762"
	my $uglyID = $_[0];
	$uglyID =~ s/ID.//;	# remove "ID=" or "ID " or "ID#"
	$uglyID =~ s/ .*//;	# remove everything from the first space to the end
	$uglyID =~ s/\W//g; # remove any non-word characters
	$uglyID =~ s/_//g;	# remove any underscores (don't get removed by \W)
	return $uglyID; 	# no longer ugly
}

sub output_dna_prot_files
{
	# DESC: receives a BioPerl DNA sequence object and outputs that DNA sequence to a (fasta) file
	####### and the corresponding protein translation to another (fasta) file
	my $out = Bio::SeqIO->new(-file => ">$out_dir".$_[0]->id.".dna.$seq_format", -format=>"$seq_format");
	$out->write_seq($_[0]);
	$out = Bio::SeqIO->new(-file => ">$out_dir".$_[0]->id.".prot.$seq_format", -format=>"$seq_format");
	$out->write_seq($_[0]->translate);
	
	# print a percentage of files processed
	if (!($i % $tenPercentStep) && 10*$curPctMult<100)
	{
		print 10*$curPctMult . "\%: $i of $n_genes valid genes have been processed ".
								"and DNA/Protein files created.\n";
		++$curPctMult;
	}
	
	++$i;
	if ($i == $n_genes)
	{
		print "100\% of $n_genes valid genes processed!\n";
	}
}

sub print_file_error
{
	# DESC: prints an error message for a file that the user specified but does not exist
	print "*****ERROR***** File \"$_[0]\" does not exist!!!\n";
}
