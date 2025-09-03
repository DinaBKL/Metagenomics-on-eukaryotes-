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


my @fnafiles = `ls -1 allgenomes/*FNA`; chomp @fnafiles;
#exit;
my %basefna ;

foreach my $f (@fnafiles) {
        my ($dir,$fname)=split(/\//,$f);
        # for mmetsep and uniprotkb these files have no results
        next if $fname eq "bin1057.fna.FNA.5KB.FNA"; # bin1057.fna.FNA.5KB.FNA
        next if $fname eq "bin48.fna.FNA.5KB.FNA";  #bin48.fna.FNA.5KB.FNA
        next if $fname eq "bin923.fna.FNA.5KB.FNA"; #bin923.fna.FNA.5KB.FNA
	next if $fname eq "bin745.fna.FNA.5KB.FNA";
        #push (@basefna,$fname);
        $basefna{$fname}=$f;
}

print Dumper \%basefna;

my @folders=qw/genomes_uniprotkbmm genomes_mmetsp genomes_uniprotkb/;


my $folder    = "genomes_mmetsp";
my $startname = "results-MMETSP";

my $folder = "genomes_uniprotkbmm";
my $startname = "results-uniprotkb_mm";

#check_results_in_folder($folder,$startname);


#my $folder    = "genomes_uniprotkb";
#my $startname = "results-uniprotkb";


#my $folder    = "genomes_uniref";
#my $startname = "results-uniref90";


my @taxfilenames = check_results_in_folder($folder,$startname); # $folder/startname-$files/taxResults_tax_per_contig.tsv


my $outfh = createFile("output_file_names_uniprotkbmm.txt");

foreach my $taxfile (@taxfilenames) {

        my ($folder1,$folder2,$taxname)=split(/\//,$taxfile);
        my ($p1,$p2,$p3)=split(/-/,$folder2);
        
        my $fnafile = "allgenomes/$p3";
        my $basefna = $p3;
        print "fnafile: $p3\n";
        print "fnafile: $fnafile\n";
        print "taxfile: $taxfile\n";
        my $outpath = "$folder1/$folder2/";
        print "outpath=$outpath\n";
    
        print $outfh "FOSMID_ID=$basefna\n";
        print $outfh "$outpath/$basefna\n";
    
        my ($cat_fna,$cat_tsv) = find_contigs ($fnafile,$basefna,$taxfile,$outpath,"d_Bacteria"); 
        print $outfh "$cat_fna\n$cat_tsv\n";
        my ($cat_fna,$cat_tsv) = find_contigs ($fnafile,$basefna,$taxfile,$outpath,"d_Archaea");
        print $outfh "$cat_fna\n$cat_tsv\n";
        my ($cat_fna,$cat_tsv) = find_contigs ($fnafile,$basefna,$taxfile,$outpath,"d_Eukaryota");
        print $outfh "$cat_fna\n$cat_tsv\n";
        my ($cat_fna,$cat_tsv) = find_contigs ($fnafile,$basefna,$taxfile,$outpath,"d_Viruses");
        print $outfh "$cat_fna\n$cat_tsv\n";
    
        my ($cat_fna,$cat_tsv) = find_unclassfied_contigs ($fnafile,$basefna,$taxfile,$outpath);
        print $outfh "$cat_fna\n$cat_tsv\n";
    
        print $outfh "//\n";
    
        #exit;

}


sub find_unclassfied_contigs {
        my ($fnafile,$basefna,$taxfile,$outpath) = @_;
        my $category_fna = "$outpath/$basefna-d_Unclassified.fna";
        my $category_tsv = "$outpath/$basefna-d_Unclassified.tsv";
        system ("cat $taxfile | grep -v d_Bacteria | grep -v d_Archaea | grep -v d_Eukaryota | grep -v Viruses | cut -f1 | faSomeRecords $fnafile stdin stdout > $category_fna");
        system ("cat $taxfile | grep -v d_Bacteria | grep -v d_Archaea | grep -v d_Eukaryota | grep -v Viruses > $category_tsv");
        print "Made $category_fna and $category_tsv\n";
        return ($category_fna,$category_tsv);

}

sub find_contigs {

        my ($fnafile,$basefna,$taxfile,$outpath,$category) = @_;
    
        my @results = `grep $category $taxfile`; chomp @results;
        my $numhits = scalar @results;
    
        my $category_fna = "$outpath/$basefna-$category.fna";
        my $category_tsv = "$outpath/$basefna-$category.tsv";
        system ("cat $taxfile | grep $category | cut -f1 | faSomeRecords $fnafile stdin stdout > $category_fna");
        system ("cat $taxfile | grep $category > $category_tsv");
        print "Made $category_fna and $category_tsv\n";
        #exit;
    
        return ($category_fna,$category_tsv);



}





sub check_results_in_folder {

	print "checking results...\n";
	    my $folder    = shift;
	    my $startname = shift;
	
	    my @taxfilenames;
	
	    foreach my $files (sort keys %basefna) {
	
	            my $expected_file = "$folder/$startname-$files/taxResults_tax_per_contig.tsv" ;
	            if (! -e $expected_file) {
	                    die "Cannot find $expected_file\n";
	            } else {
	                    print "Found $expected_file\n";
	            }
	            push (@taxfilenames,$expected_file);
	    }
	print "If we are here all files passed the check\n";
	    return @taxfilenames;

}




__END__

