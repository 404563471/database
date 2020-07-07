#!/usr/bin/perl -w 
# by Yan Pengcheng yanpc@bcc.ac.cn
# Perl wrapper of edgeR 
# 2013-01-23

=head1 Program Description

perl DGE.pl  [options]

        -d,--data	<string>        expression data file
        -g1,--group1	<string>	sample name in same group [s1 s2 ..]
        -g2,--group2	<string>	sample name in same group [s3 s4 ..]
	-p,--pvalue	<float>		pvalue of test	[0.01]
	-f,--fdr	<float>		FDR [0.01]
	-l,--logFC	<int>		logFC [1]
        -n,--normalize	<string>	normalization method [quantile] in RPKM,TMM,RLE,quantile
        -o,--output	<string>	output file
        -h,--hhelp               output help information to screen

=cut

use strict;
use Getopt::Long;
use Pod::Text;
use File::Basename qw(basename dirname);
use Data::Dumper;

my $EXEC = "/share/software";
my %method = ("TMM"=>1,"RLE"=>1,"FPKM"=>1,"quantile"=>1,"no"=>1);
my ($R_PATH, $HELP, $pvalue, $fdr, $logFC) = ("$EXEC/bin", "", 0.01, 0.01, 1);
my ($output, @g1, @g2, $norm, $data);
$norm = "TMM";

GetOptions(
        "d|data:s"=>\$data,
        "g1|group1:s{,}"=>\@g1,
        "g2|group2:s{,}"=>\@g2,
	"o|output:s"=>\$output,
	"p|pvalue:s"=>\$pvalue,
	"f|fdr:s"=>\$fdr,
	"l|logFC:s"=>\$logFC,
	"n|normalize:s"=>\$norm,
        #"v|verbose!"=>\$VERBOSE,
        "h|help!"=>\$HELP
);

die `pod2text $0` if($HELP || @g2 <= 0 || @g1 <= 0);
die `pod2text $0` if(!$method{$norm});

if(!-e $output){
	mkdir($output, 0755);
}

for(my $i=0;$i<@g1;$i++){
	$g1[$i] = "X".$g1[$i] if($g1[$i] =~ /^\d/);
	$g2[$i] = "X".$g2[$i] if($g2[$i] =~ /^\d/);
}

my $groups = join(",",@g1).",".join(",",@g2);
my $df = "c(\"".join("\",\"",@g1)."\",\"".join("\",\"",@g2)."\")";
$df =~ s/-/./g;
my $RG = join(",", map {"G1"} @g1).",".join(",", map {"G2"} @g2);
&edgeR($data, $output, $RG, $norm);

sub edgeR(){
	my ($infile, $title, $groups, $norm) = @_;
	$groups =~ s/,/","/g;
	if($norm eq "FPKM"||$norm eq "no"){
		$norm = "d";
	}else{
		$norm = "calcNormFactors(d, method=\"$norm\")";
	}

	my $Rscripts=qq|

options(digits=2);
library("limma",lib.loc = "/var/www/html/fish-v4/bin/lib64/R/library");
library("edgeR",lib.loc = "/var/www/html/fish-v4/bin/lib64/R/library");
x <- read.table("$infile", header=TRUE, sep="\\t", row.name=1);
df<-x[,$df];
#pdf(file="$title/dif.pdf");
png(file="$title/dif.png");
d<-DGEList(counts=df[apply(df,1,sum)>0,], group=c("$groups"));
f<-$norm;
cmdisp<-estimateCommonDisp(f);
if(is.na(cmdisp\$common.dispersion)){ cmdisp\$common.dispersion = 0.3 }
de<-exactTest(cmdisp);
h <- t(data.frame(GID=names(de\$table)));
write.table(h,file="$title/dif-all.xls",sep="\\t",quote=FALSE,col.names=F)
write.table(de\$table,file="$title/dif-all.xls",append=T, sep="\\t",quote=FALSE, col.names=F)

top<-topTags(de,n=sum(de\$table\$PValue < '$pvalue'));
topFDR<-subset(top\$table,FDR < $fdr & abs(logFC) >= $logFC);
topFDR\$trend = "UP";
topFDR[topFDR\$logFC < 0,]\$trend = "DOWN";
detags<-rownames(topFDR)
topFDR<-cbind(df[detags,],topFDR)

plotSmear(de, de.tags=detags, main="FC plot using common dispersion");
abline(h=c(-1,1), col="dodgerblue");

h <- t(data.frame(GID=names(topFDR)));
write.table(h,file="$title/dif.xls",sep="\\t",quote=FALSE,col.names=F);
mode(topFDR);
class(topFDR);
write.table(topFDR,file="$title/dif.xls",append=T,sep="\\t",quote=FALSE,col.names=F)
dev.off();
|;

#top<-topTags(de,n=sum(de\$table\$p.value < '0.01'));
#plotSmear(d, de.tags=detags, xlab="M",ylab="A", main="FC plot using common dispersion");
#, lib.size=colSums(x));
	open(RSRC,">$title.R")||die "Can not write R scripts\n";
	print RSRC $Rscripts;
	close RSRC;

	system("/var/www/html/fish-v4/bin/bin/R CMD BATCH --no-save $title.R");
	#system("rm $title.R*");
}
