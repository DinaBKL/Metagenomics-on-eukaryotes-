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


my @fnas = <*.fna>; chomp @fnas;
my @metagenomes = <*.fasta.gz>; chomp @metagenomes;


foreach my $f (@fnas) {


=pod

time ./hisat2/hisat2 -x MAG-Haptophyta-Group1-Tiefwaren-Coassembly-12-bins.fna.h2 -1 GROS21-1_SAMEA9560375.fasta.gz.1.fa.gz -2 GROS21-1_SAMEA9560375.fasta.gz.2.fa.gz --sensitive  --no-spliced-align
ment -p 80 --summary-file sample_summary_sensitive2.txt -S sample_sensitive2.sam -f &

=cut

        foreach my $m (@metagenomes) {

                my $cmd = "./hisat2/hisat2 -x $f.h2 -1 $m.1.fa.gz -2 $m.2.fa.gz --sensitive --no-spliced-alignment -p 10 --summary-file $f--VS--$m.summary.txt -S /dev/null -f ";
                print "$cmd\n";

        }

}
