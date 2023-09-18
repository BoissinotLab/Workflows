#!/usr/bin/perl

#############################################################################################
#############################################################################################
#	This script extracts protein-coding sequences (DNA), 
#	and the corresponding protein sequences from a Genbank or Embl format file, or files.
#	This script requires Bioperl.
#	
#	INPUT: Genbank or Embl format file(s). 
#	
#	OUTPUT: fasta format files named: $species.dna.fasta, $species.prot.fasta 
#	The script also outputs: pseudogene sequences named: $species.prot-pseudogenes.fasta, 
#	$species.dna-pseudogenes.fasta
#
#	Because Embl and Genbank files are not completely standardised, you will need to  
#	define the 'tag' name genes. In some species this is called 'systematic_name',
#	in others it is 'Gene' etc.
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
##############################################################################################
#############################################################################################
#############################################################################################

#Load modules
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Getopt::Long;

#Define global variables
my ($inputfiles,$format,$species,$primarytag,$help);
my $script_function = "\n*********\nThis script extracts CDS sequences and protein sequencees from a Genbank or Embl file.\n*********\n";
my $usage =		"$script_function".
			"\n*********\n* Usage *\n*********\n".
			"\nperl $0 -i inputfiles -s species_name -t tag*".
			"\n\t[-f genbank/embl (embl default)] [-help]".
			"\n\t*This tag will differ depending on the species/input files. 'systematic_id' is one example.\n";

#Define the command line options
GetOptions (
	"h" 		=> \$help,			#if you need it
	"i=s" 		=> \$inputfiles,	#eg: "*.embl"
	"s=s" 		=> \$species,		#eg: human
	"f=s" 		=> \$format,		#eg: genbank/embl
	"t=s"		=> \$primarytag		#eg: systematic_id
);

#check the command line options
die $usage if $help;
die "$usage\nERROR!: You must specify input files! (eg: \"*.genbank\").\n\n" unless ($inputfiles);
die "$usage\nERROR!: You must define species! (eg: human).\n\n" unless ($species);
die "$usage\nERROR!: You must define tag! (eg: systematic).\n\n" unless ($primarytag);

#default format for files
$format ||= 'embl';
my $pseudogene_count=0;;
my $gene_count=0;

#make a log file
my $date = timestring();
open (LOG, ">proteins_from_embl.$date.log");
print LOG "Command: $0 -i $inputfiles -s $species -t $primarytag -f $format\n";
print LOG "Date: $date\n";

#collect the list of input files
my @files = glob("$inputfiles");
print LOG "Genome files: @files\n";

#define the output files
my $dna_seq_out = Bio::SeqIO->new ('-format' => 'fasta', '-file' =>">$species.dna.fasta");
my $pro_seq_out = Bio::SeqIO->new ('-format' => 'fasta', '-file' =>">$species.prot.fasta");
my $dna_seq_out_pseud = Bio::SeqIO->new ('-format' => 'fasta', '-file' =>">$species.dna-pseudogenes.fasta");
my $pro_seq_out_pseud = Bio::SeqIO->new ('-format' => 'fasta', '-file' =>">$species.prot-pseudogenes.fasta");

#iterate over the input files, collect CDS and output DNA and protein seqs.
foreach my $file (@files){
	my $input_seqs = Bio::SeqIO->new ('-format' => "$format", '-file' =>"$file");
	warn "loading genome file $file\n";
	while (my $s = $input_seqs->next_seq()){

		my $seqid = $s->id();
		warn "processing CDS seqs in $seqid\n";
		
		#Runs through each sequence feature
    	foreach my $f ($s->all_SeqFeatures()){

			#Start of gene loop
			#Skip unless the features is a coding sequence (CDS)
			next unless ($f->primary_tag eq 'CDS');
			
			#collect the $systematic gene name (if present)
			my $systematic = 'NA';
			if ($f->has_tag($primarytag)){
				my @gene_names = $f->get_tag_values($primarytag);
				$systematic = join("|",@gene_names); 
			}
			#if the $systematic gene name is not present, use chromosome.start.stop as the name
			else {
				my $start = $f->start();
				my $end = $f->end();
				$systematic = "$seqid.$start.$end";
			}
			
			#collect the gene (spliced, if required), with correct strand
			my $dna = $f->spliced_seq();
			my $pro = $f->spliced_seq->translate();
			
			#rename seqs with $systematic
			$pro->id($systematic);
			$dna->id($systematic);

			#check for pseudogenes
			my $pseudogene = check_for_pseudo($pro);
			
			#if its a pseudo gene put into different file & warn
			if ($pseudogene){
				$pseudogene_count++;
				print LOG ">$systematic|pseudogene|$pseudogene\n";
				$dna_seq_out_pseud->write_seq($dna);
				$pro_seq_out_pseud->write_seq($pro);

			}
			
			#otherwise output seqs
			$gene_count++;
			$dna_seq_out->write_seq($dna);
			$pro_seq_out->write_seq($pro);
		}
	}
	warn "done with genome file $file\n";
}

#complete, and exit
warn "Done. Outputs are:
$species.dna.fasta
$species.prot.fasta
$species.prot-pseudogenes.fasta
$species.dna-pseudogenes.fasta

Genes: $gene_count
Pseudogenes: $pseudogene_count
\n\n";

print LOG "Done.
Genes: $gene_count
Pseudogenes: $pseudogene_count
\n\n";

exit;



#############################################################################################
#############################################################################################
# SUBROUTINES
#############################################################################################
#############################################################################################

sub timestring {
        my  $now = localtime;
        my @t=split(' ',$now);
        return "$t[2]-$t[1]-$t[4]";
}


sub check_for_pseudo {
	my ($pro) = @_;
	my $pseudogene;
	my $pro_seq = $pro->seq;
	unless ($pro_seq =~/^M/){
		$pseudogene .= 'incorrect_start ';
	}
	unless($pro_seq =~/\*$/){
		$pseudogene .= 'incorrect_end ';
	}
	my $stop_count = ($pro_seq =~ tr/\*/\*/);
	if ($stop_count > 1){
		$pseudogene .= 'internal_stop_codons ';
	}
	return $pseudogene;
}

