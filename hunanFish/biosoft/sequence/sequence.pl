#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my $sample="test.fq";
my $spece="rna";
my $adaptor="AGATCGGAAGAGCACACGTCT";
my $outdir="out";
my $pic_path="picture";
my $format ="fasta or fastq";
my $flag_adaptor;
GetOptions(
	"s|sample:s" =>	\$sample,
	"a|adaptor:s"=>	\$adaptor,
	"o|outdir:s" =>	\$outdir,
	"f|format:s" =>	\$format
);
my $id;
if($sample=~m/(.+?)\./){
	$id=$1;
}
if(not -d "out"){
	mkdir 'out';
}
chdir 'out';
mkdir $outdir;
chdir $outdir;
mkdir $pic_path;
my $fastqc="/var/www/html/microRNA/bin/fastqc";
my $cutadapt="/var/www/html/microRNA/software/cutadapt-1.4.2/bin/cutadapt";
my $q2a="/var/www/html/microRNA/pipeline/sequence/FastqToFasta.pl";
if($format eq "fastq"){
	system("ln -s ../../$sample ./");
#	print "$sample";
	system("$fastqc $sample") ;
	if($adaptor ne ""){
		system("$cutadapt -a $adaptor --discard-untrimmed -m 10 $sample -o reads_mv_Adaptor.fastq");
		system("$q2a -i reads_mv_Adaptor.fastq -o reads.fa");
	}else{
		system("$q2a -i $sample -o reads.fa") ;
	}
}
if($format eq "fasta"){
	system("ln -s  ../../$sample reads.fa");
}
system("perl ../../getSeqLength.pl reads.fa > length_$id.txt");
system("Rscript ../../reads_len.R $id $pic_path");
system("perl ../../collapse_reads_md.pl reads.fa $spece -a > reads_md.fa");
system("/var/www/html/microRNA/bin/blastall -p blastn -d ../../Rfam.fa -i reads_md.fa -o reads_md_Rfam.out -a 12 -m 8 -e 0.01 -F F");
system("perl ../../getBlastDB.pl -c reads_md.fa -b  reads_md_Rfam.out -id $id -iden 100 -cov 1.0 -o Rfam_$id");
system("perl ../../getBaseBias.pl miRNA_predict_$id.fa > FirstBaseBias_$id.txt");
system("perl ../../getBaseBiasAll.pl miRNA_predict_$id.fa > AllBaseBias_$id.txt");
system("Rscript ../../getBasebias.R $id $pic_path");



print "$sample\n$id\n";

my @files=glob "picture/*.pdf";
foreach my $tmp(@files){
	my $before=$tmp;
	$tmp=~s/pdf/png/;
	print "before $before, after $tmp \n";
	system("convert $before $tmp");
}
