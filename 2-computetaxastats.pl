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
use Sort::Hash;

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

my $USAGE=<<USAGE;

  * $0 *

   given a file containing a list of Eukaryota.tsv files summarizes the entire file 
   to find the most frequent taxon name, and the taxon name that recovers the maximum
   eukaryotic sequence

   input file format is standardized and like this
   output_file_names_uniprotkb.txt

FOSMID_ID=bin1.fna.FNA.5KB.FNA
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Bacteria.fna
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Bacteria.tsv
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Archaea.fna
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Archaea.tsv
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Eukaryota.fna
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Eukaryota.tsv
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Viruses.fna
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Viruses.tsv
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Unclassified.fna
genomes_uniprotkb/results-uniprotkb-bin1.fna.FNA.5KB.FNA//bin1.fna.FNA.5KB.FNA-d_Unclassified.tsv
//


  options:

     -i,-infile        input file


​    

      -cutoff          1000000 if a taxon length is less than 1MB
    
     -h,-help          show this help

   usage:

   $0 -infile output_file_names_uniprotkb.txt > uniprot_genome_taxa_lengths.txt 



USAGE

my $CUTOFF=1000000;

##############################################################################################################################

## declare option variables or others

my ($infile,$outfile,$help,$informat,$outformat,$format,$prefix,$minlength,$minidentity);
my $jobs=20;
GetOptions ("i|infile=s"        => \$infile,
            "o|outfile=s"       => \$outfile,
	    "c|cutoff=f"        => \$CUTOFF,
            "p|prefix=s"        => \$prefix,       # s for strings , values like infile.txt 
            "j|jobs=i"          => \$jobs,         # i for integers, values like 3,4,5 (perl will check for valid integers)
            "mi|minidentity=f"  => \$minidentity,  # f for floats, values like 0.1  (perl will check for valid float)
            "ml|minlength=f"    => \$minlength,
            "h|help"            => \$help          # simple switch if given value is 1 if omitted value is 0. can be tested. if ($help) { print "here it is $USAGE";}
     );

if (!$infile)  { die "$USAGE No input file (-i) "        };
#if (!$outfile) { die "$USAGE No output file (-o) "       };
#if (!$prefix)  { die "$USAGE No prefix specified (-p) "  };
if ($help)     { die $USAGE };


#$/="//\n";

my $filedata=`cat $infile`; 
#print Dumper $filedata;
my @data_array=split(/\/\/\n/,$filedata);
#print Dumper @data_array; exit;

##############################################################################################################################

## using filehandles to read and create files

my $infh=readFile($infile);
#my $outfh=createFile($outfile);

foreach my $record (@data_array) {
        next if $record !~ /^FOSMID/;
	#print Dumper $record;
	my @array=split(/\n/,$record); chomp @array;
	my @fnas = grep /Eukaryota.fna/, @array;
	my @tsvs = grep /Eukaryota.tsv/, @array;
	my $fnafile = shift @fnas; my $tsvfile = shift @tsvs;
	if (! -e $fnafile) { 	die "Cannot find fna $fnafile\n"; 	}
        if (! -e $tsvfile) {    die "Cannot find tsv $tsvfile\n";       }
	#next;
	# load lengths of all contigs
	#my %lenhash = load_fasta_lengths($fnafile,"hash");
	# load all taxon names
	#my %taxon_names = load_all_taxon_names ($tsvfile);
	

	my ($total_genome_length,$unidentified_genome_length)= compute_unidentified_genome_length($fnafile,$tsvfile);
	#print "$tsvfile\t$total_genome_length\t$unidentified_genome_length\n";
	my $genome_length_per_taxon=compute_taxa_length($fnafile,$tsvfile);
	$genome_length_per_taxon=~s/-_//g;
	
	    my $top_categories = compute_top_categories($fnafile,$tsvfile);
	    $top_categories=~s/-_//g;
	
	print "$fnafile\t$tsvfile\t$total_genome_length\t$unidentified_genome_length\t$top_categories\t$genome_length_per_taxon\n";

#	print "$fnafile\t$tsvfile\t$top_categories\n";

	#exit;

}

$infh->close;
#$outfh->close;


sub compute_top_categories {

	#-_cellular organisms;d_Eukaryota;p_Haptista;c_Haptophyta;-_Prymnesiophyceae;o_Prymnesiales;
	my $fnafile = shift;
	my $tsvfile = shift;
	
	    my %lenhash = load_fasta_lengths($fnafile,"hash");
	
	my $result=`cut -f1,9 $tsvfile`;
	my @data=split(/\n/,$result); chomp @data;
	
	my %categories;
	
	foreach my $line (@data) {
		my ($contig,$taxonomy)=split(/\t/,$line);
		next if $taxonomy eq "-_cellular organisms;d_Eukaryota";
		my @array=split(/;/,$taxonomy);
		shift @array; shift @array;
		my $cat1 = $array[0];
		my $cat2 = "$cat1;$array[1]";
		my $cat3 = "$cat2;$array[2]";
		my $cat4 = "$cat3;$array[3]";
		my $cat5 = "$cat4;$array[4]";	
		my $cat6 = "$cat5;$array[5]";
		my $cat7 = "$cat6;$array[6]";
		
		$categories{$cat1} += $lenhash{$contig};
		$categories{$cat2} += $lenhash{$contig};
		$categories{$cat3} += $lenhash{$contig};
		$categories{$cat4} += $lenhash{$contig};
		$categories{$cat5} += $lenhash{$contig};
	            $categories{$cat6} += $lenhash{$contig};
	            $categories{$cat7} += $lenhash{$contig};



	}
	
	my @lengths;
	foreach my $cat (sort keys %categories) {
		my $line = "$cat\t$categories{$cat}";
		push (@lengths,$line);
	}
	
	my @sorted = fieldsort '\t',['-2n'], @lengths;
	
	my $top1 = shift @sorted;
	my $top2 = shift @sorted;
	my $top3 = shift @sorted;
	my $top4 = shift @sorted;
	my $top5 = shift @sorted;
	    my $top6 = shift @sorted;
	    my $top7 = shift @sorted;


	my $result;
	if ($top1 eq "") { $top1 = "_none_\t0"; }
	    if ($top2 eq "") { $top2 = "_none_\t0"; }
	    if ($top3 eq "") { $top3 = "_none_\t0"; }
	    if ($top4 eq "") { $top4 = "_none_\t0"; }
	    if ($top5 eq "") { $top5 = "_none_\t0"; }
	    if ($top6 eq "") { $top6 = "_none_\t0"; }
	    if ($top7 eq "") { $top7 = "_none_\t0"; }



	my $result="[C1]$top1\t[C2]$top2\t[C3]$top3\t[C4]$top4\t[C5]$top5\t[C6]$top6\t[C7]$top7";
	return $result;

#	return \%categories;

}


sub compute_taxa_length {

	my $fnafile = shift;
	my $tsvfile = shift;
	
	#local $/="\n";
	
	    # load lengths of all contigs
	    my %lenhash = load_fasta_lengths($fnafile,"hash");
	    # load all taxon names, ignoring cellular organisms and d_Eukaryota
	    my %taxon_names = load_all_taxon_names ($tsvfile);
	#print Dumper \%taxon_names;
	
	# now go thru contig by contigs
	    my $result=`cat $tsvfile | cut -f1,9`; #chomp @data;
	    my @data=split(/\n/,$result);
	#print Dumper @data;exit;
	
	# set length of all taxa to 0
	my %taxon_lengths;
	foreach my $taxon (keys %taxon_names) {
		$taxon_lengths{$taxon}=0;
	}
	
	#print Dumper \%taxon_lengths; exit;
	
	    foreach my $line (@data) {
	            chomp $line;
	            my ($contig,$taxonomy)=tabsplit($line);
	            $contig=trim($contig); $taxonomy=trim($taxonomy);
	            #print "$contig,tax=$taxonomy\n";
	            #next if ($taxonomy eq "-_cellular organisms;d_Eukaryota");
		my @names = split (/;/,$taxonomy);
		foreach my $n (@names) {
			$n=trim($n);
		}
	
		foreach my $n (@names) {
			$taxon_lengths{$n} += $lenhash{$contig};
		}
	    }
	
	# print Dumper \%taxon_lengths; exit;
	# delete useless names
	
	delete $taxon_lengths{"-_cellular organisms"};
	delete $taxon_lengths{"d_Eukaryota"};
	
	# filter taxon lengths, if its less than 500,000 bp whats the point
	foreach my $t (sort keys %taxon_lengths) {
		if ($taxon_lengths{$t}<$CUTOFF) {
			delete $taxon_lengths{$t};
		}
	}
	
	# print Dumper \%taxon_lengths; exit;
	# now sort these names by contig lengths
	
	#my @sorted = sort_hash("numeric","desc",\%taxon_lengths);
	#foreach my $t (@sorted) {
	#	print "$t\t$taxon_lengths{$t}\n";
	#}
	#exit;
	#print Dumper @sorted; 
	
	my @t_lengths;
	foreach my $taxon (sort keys %taxon_lengths) {
		my $line = "[T]$taxon\t$taxon_lengths{$taxon}";
		#print "adding line=$line\n";
		push (@t_lengths,$line);
	}
	

#	exit;

#	print Dumper @t_lengths; exit;

	my @sorted_lengths = fieldsort '\t', ['-2n'], @t_lengths;
	

#	print Dumper @sorted_lengths; 

#	exit;

	my $final_result = join ("\t",@sorted_lengths);
	return $final_result;
	


}



sub compute_unidentified_genome_length {
	my $fnafile = shift;
	my $tsvfile = shift;
	

        my %lenhash = load_fasta_lengths($fnafile,"hash");
    my $total_genome_length=0;
    foreach my $c (%lenhash) { $total_genome_length += $lenhash{$c}; }
    my $numcontigs=scalar keys %lenhash;
    #print "total_length=$total_genome_length, numcontigs=$numcontigs\n";
    
    my $unidentified_genome_length=0;
        my $result=`cat $tsvfile | cut -f1,9`; #chomp @data;
    my @data=split(/\n/,$result);
    
    #print Dumper @data;
    my $num=scalar @data;
    #print "Found $num results \n"; 
    foreach my $line (@data) {
    	chomp $line;
    	my ($contig,$taxonomy)=tabsplit($line);
    	$contig=trim($contig); $taxonomy=trim($taxonomy);
    	#print "$contig , tax=$taxonomy\n";
    	if ($taxonomy eq "-_cellular organisms;d_Eukaryota") {
    		$unidentified_genome_length += $lenhash{$contig};
    	}
    	
    }
    
    return ($total_genome_length,$unidentified_genome_length);
    

}


sub load_all_taxon_names {
	my $file=shift;
	my $result=`cat $file | cut -f9`; 
	my @data=split(/\n/,$result); chomp @data;
	my %taxon_hash;
	#-_cellular organisms;d_Eukaryota
	foreach my $element (@data) {
		my @array=split(/;/,$element);
		foreach my $e (@array) {
			next if $e eq "-_cellular organisms";
			next if $e eq "d_Eukaryota";
			$e=trim($e);
			$taxon_hash{$e}++;
		}
	}
	return %taxon_hash;

}

__END__

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

