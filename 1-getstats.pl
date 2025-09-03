#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case);
use Array::Utils;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Sort::Fields;   # fieldsort function for sorting tab-delimited data in arrays
use Try::Tiny;      # simplified try/catch error catching
use File::chdir;    # allows use of the chdir function
use File::Path qw(make_path remove_tree);    # make and delete folders
use File::Copy;                              # move() and copy() functions
use Cwd;                                     # provides the getcwd function
use IO::Tee;           # allows writing a logfile using tee (see below)
use File::Basename;    # allows getting proper paths to files if needed

use FindBin;           #  allows script to find its own directory allowing packaging of other scripts, data etc
FindBin::again();      #  do not remove this 

##### correctly locate the lib folder with modules (if you have one)

use lib "$FindBin::Bin/lib";  # use if needed 
=pod
#remove the pod cut if you want to include a custom lib 
use mylib; # include ur functions in mylib.pm into the program, can even package BioLib2 here
=cut

use BioLib2; # use if you know its installed on server or package it in your lib folder 

my $BASEDIR = $FindBin::Bin;    # use findbin to get the folder where the script is invoked from

##############################################################################################################################

## using filehandles to read and create files

my $infile=shift;
my $outfile=shift;

my $USAGE=" $0 \n $0 infile outfile\n";
if (!$infile) {die "no input file\n$USAGE"}
if (!$outfile) {die "no output file\n$USAGE"}

$/="\/\/\n";

my $infh=readFile($infile);
my $outfh=createFile($outfile);
my $P1="#bin\to_numcontigs\to_n50\to_max\to_mean\to_length";
my $P2="b_numcontigs\tb_n50\tb_max\tb_mean\tb_length";
my $P3="a_numcontigs\ta_n50\ta_max\ta_mean\ta_length";
my $P4="e_numcontigs\te_n50\te_max\te_mean\te_length";
my $P5="v_numcontigs\tv_n50\tv_max\tv_mean\tv_length";
my $P6="u_numcontigs\tu_n50\tu_max\tu_mean\tu_length";

my $header="$P1\t$P2\t$P3\t$P4\t$P5\t$P6";
print $outfh "$header\n";

while (my $record = <$infh>) {
        #print Dumper $record;
        my @array = split(/\n/,$record);
        my $header=shift @array;
        my ($p1,$bin_name)=split(/=/,$header);

        my $original_fna = $array[0];
        my $bacteria_fna = $array[1];
        my $archaea_fna  = $array[3];
        my $euk_fna      = $array[5];
        my $virus_fna    = $array[7];
        my $unc_fna      = $array[9];
    
        my ($o_numcontigs,$o_n50,$o_max,$o_mean,$o_total_length)=compute_stats_on_fna($original_fna);
        my ($b_numcontigs,$b_n50,$b_max,$b_mean,$b_total_length)=compute_stats_on_fna($bacteria_fna);
        my ($a_numcontigs,$a_n50,$a_max,$a_mean,$a_total_length)=compute_stats_on_fna($archaea_fna);
        my ($e_numcontigs,$e_n50,$e_max,$e_mean,$e_total_length)=compute_stats_on_fna($euk_fna);
        my ($v_numcontigs,$v_n50,$v_max,$v_mean,$v_total_length)=compute_stats_on_fna($virus_fna);
        my ($u_numcontigs,$u_n50,$u_max,$u_mean,$u_total_length)=compute_stats_on_fna($unc_fna);


        my $P1="$o_numcontigs\t$o_n50\t$o_max\t$o_mean\t$o_total_length";
        my $P2="$b_numcontigs\t$b_n50\t$b_max\t$b_mean\t$b_total_length";
        my $P3="$a_numcontigs\t$a_n50\t$a_max\t$a_mean\t$a_total_length";
        my $P4="$e_numcontigs\t$e_n50\t$v_max\t$e_mean\t$e_total_length";
        my $P5="$v_numcontigs\t$v_n50\t$e_max\t$v_mean\t$v_total_length";
        my $P6="$u_numcontigs\t$u_n50\t$u_max\t$u_mean\t$u_total_length";
    
        my $result = "$bin_name\t$P1\t$P2\t$P3\t$P4\t$P5\t$P6";
        print $outfh "$result\n";
    
        #<stdin>;


​       
}


sub compute_stats_on_fna {
        my $file=shift;
        if (! -e $file) {
                die "expected to find file $file but its not there!\n";
        }

        my $num_contigs=0;
        my $n50 = 0;
        my $max = 0;
        my $mean = 0;
        my $total_length=0;
        if (! -z $file) {
           # if file isnt empty get values
                my $stats = `raw_n50 $file`; chomp $stats;
                my @array=split(/\s+/,$stats);
                $num_contigs=$array[1];
                $n50=$array[3];
                $max=$array[5];
                $mean=$array[7];
                $total_length=$array[10];
        }
    
        return ($num_contigs,$n50,$max,$mean,$total_length);
    
        # 0        1    2     3   4      5     6      7    8      9       10       11  12
        #contigs: 3725 n50: 26816 max: 176265 mean: 18556 total length: 69124227 n80: 11852


        # check if file is emptyy



}


__END__


  FOSMID_ID=bin1.fna.FNA.5KB.FNA
0 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA
1 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Bacteria.fna
2 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Bacteria.tsv
3 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Archaea.fna
4 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Archaea.tsv
5 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Eukaryota.fna
6 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Eukaryota.tsv
7 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Viruses.fna
8 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Viruses.tsv
9 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Unclassified.fna
10 genomes_mmetsp/results-MMETSP-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Unclassified.tsv
//

$infh->close;
$outfh->close;

##############################################################################################################################

## standard file read in and out: use when you dont want to use any biolib funxtions. usu. not used 

open (INP, "$infile") or die "Cannot open input file $infile! $!
";
open (OUT, ">$outfile") or die "Cannot open output file $outfile! $!
";
while (my $line = <INP>) {
        chomp $line;
        $line = trim ($line);
        next if $line eq "";
        next if $line =~ /^#/;
        my @array = tabsplit ($line);


}
close INP;
close OUT;

##############################################################################################################################

## bioperl file read in and out

my $seqin  = Bio::SeqIO->new(-file => "$infile", -format=>"fasta");  # other formats are allowed, embl, genbank etc 
my $seqout = Bio::SeqIO->new(-file =>">$outfile", -format=>"fasta");

while (my $inseq = $seqin->next_seq) {
  my $id   = $inseq->display_id;
  my $desc = $inseq->desc;
  my $len  = $inseq->length;
  my $seq  = $inseq->seq;

  # change the id, desc whatever 

  my $newid=$id."new";

  # change in the object

  $inseq->id($newid);

  # write object to output stream

  $seqout->write_seq($inseq);

  # or just print out as you prefer 

  print "$id    $len    $desc\n";
}


##############################################################################################################################

## getting a list of files (less than 80k )

# gets a list of files with .fasta extension from current directory (can check with getcwd)

my $cwd = getcwd();

my @files = <*.fasta>; chomp @files; 
foreach my $f (@files) {
	# do something ...
	my @data = `cat $f`;  # grab file contents in an array
        foreach my $line (@data) {
		# do something 
        }
}

##############################################################################################################################

## getting data from a file quickly without using filehandles, use backticks 

my $filename="infile.txt";
my @data=`cat $filename`; # each line goes in a separate array element, beware of the 
 that exists at the end of each line 
chomp @data ; # remove the \n at the end of each line in the array 
foreach my $line (@data) {
	# do something 
}

my $dataline=`cat $filename`; # entire file contents go in a single variable 

# can split it 

my @filedata = split(/\n/,$dataline); # now every line is an array element 
foreach my $line (@filedata) {
	# do something 
}


my @files = <*.fasta>; chomp @files; 
foreach my $f (@files) {
	# do something ...

}


##############################################################################################################################

## running a program on command line 

my $command = "cat infile.fasta | lenseq > lengths.txt";
system($command);
if (-e "lengths.txt") { print_info("Looks like the file lengths.txt exists now"); } # for more paranoia check if it has expected number of lines if you know what to expect 

# or just 

system ("cat infile.fasta | lenseq > lengths.txt")
if (-e "lengths.txt") { print_info("Looks like the file lengths.txt exists now"); }

# or even using backticks 

my @lengths=`cat infile.fasta | lenseq `; # instead of making a file you get the data from lengths.txt into an array 

# check how it looks like :

print Dumper @lengths;

