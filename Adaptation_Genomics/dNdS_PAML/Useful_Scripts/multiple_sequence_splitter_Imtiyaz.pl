#!/usr/bin/perl

###############################################################################################
###############################################################################################
#	This script Splits a FASTA file with many sequences into one file per sequence.
#	It would be wise to output all the files to an new directory using the --output (-o) flag.
#	INPUT
#	A fasta format file with many sequences (ususally genes).
#	OUTPUT
#	One fasta file per sequence.
#	REQUIREMENTS
#	(+)BioPerl (www.bioperl.org)
#
#	USAGE AND INSTRUCTIONS
#   For usage and instructions please run without any options or specifying flag -h or --help
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

# LIBRARIES & MODULES
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;


# text to display usage information

my $script_function = "\n*********\nThis script splits a FASTA file with many sequences into one file per sequence
It also checks that each output is unique (skips those that are not)\n*********\n";

my $usage = 		"$script_function".
			"\n*********\n* Usage *\n*********\n".
			"perl $0 < in.fasta".
			"\n\t [-i (--input) fasta_file] (default STDIN)".
			"\n\t [-o (--output) directory] (default working dir)".
			"\n\t [-s (--simplify)*".
			"\n\t *takes only first space-delim field from sequence id (probably gene name), removes all but A-Z0-9 characters".
			"\n\n";

 


# COMMAND-LINE INPUT
my ($input_fasta_file,$out_dir,$simplify,$help);
Getopt::Long::Configure ("bundling"); # to allow case sensitivity in single-letter options
GetOptions(	
	'input|i=s' 	=> \$input_fasta_file,
	'out|o=s'		=> \$out_dir,
	'simplify|s'	=> \$simplify,
	'help|h'		=> \$help
	);

#defaults
$out_dir ||= ".";


# check if user only wants to display the usage information, or other die conditions
if ($help) { die $usage; }
#check if output die exists
if ($out_dir){
	unless (-e $out_dir){ die "$usage\nERROR! Output directory [$out_dir] not present!\n\n"}
}

#hash to recall unique list of outputs
my %outputlist = ();


#define input object
my $input_seqs = Bio::SeqIO->new ('-format'=> 'fasta','-fh'=>\*STDIN);
if ($input_fasta_file){
	#complain if the $input_fasta_file is not present
	unless (-e $input_fasta_file){
		die "$usage.\nERROR! Can't locate the input file [$input_fasta_file]!\n\n"
	}
	#if it is present change $input_seqs object
	$input_seqs = Bio::SeqIO->new ('-format'=> 'fasta','-file'=>"$input_fasta_file");

}
#define out locations & object
my $output_seqs = Bio::SeqIO->new ('-format'=> 'fasta','-fh'=>\*STDOUT);	

#do the input/output
while (my $s = $input_seqs->next_seq()){
    my $id = $s->id();

    #skip if this id seen before
    if ($outputlist{$id}){
    	warn "This id seen before [$id]! Skipping this one.";
    }
    
    #if not seen before, output the fasta
    else {
    	#simply if told to
    	if ($simplify){
    		$id =~ s/[^a-zA-Z0-9]//g;
    		$s->display_id($id);
			$s->description('');
    	}
    	my $output_seq = Bio::SeqIO->new('-format'=> 'fasta','-file'=>">$out_dir/$id.fasta");
    	$output_seq->write_seq($s);
    	$outputlist{$id} = 1;
    }
}

#report
my $number = scalar keys %outputlist;
warn "\nDone. Output $number files to directory [$out_dir].\n\n";






