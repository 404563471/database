#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage=qq{
Usage:perl $0  -c <reads_md_collapse.fa> -b <blast_out_m8> -id <ID> -ide <identity> -cov <coverage> -o <output> -h
	    -c|col <string>:the infile contain unique reads,cut the collapsed reads      collapse后的文件
	    -b|blast <string>:the output result form blast in m8 format
	    -id <string> the sample name
	    -iden <numeric> a numeric between 0 and 100(default 90)
	    -cov  <numeric> a numeric between 0 and 1.0(default 0.9)
	    -o|output 
	    -h|help print the help information
};

my $colFile;
my $blastFile;
my $outputFile;
my $help;
my $id;
my $my_cov=0.9;
my $my_iden=0.9;

GetOptions(
		"c|col=s" =>\$colFile,
		"b|blast=s" =>\$blastFile,
		"id=s" =>\$id,
		"cov" =>\$my_cov,
		"iden" =>\$my_iden,
		"o|output=s" =>\$outputFile,
		"h|help=s" =>\$help,
		);
die ($usage)if(!defined($colFile) ||!defined($blastFile) || !defined($id) || !defined($outputFile));
die ($usage) if(defined($help));
#print $colFile,"\n";
##################################################################################

open(COL,"<$colFile") or die ("Can not open the file.$!");
my ($seq_name,$NO,$len);
my ($total,$total_collapse)=(0,0);
my (%hash,%colla);
while(<COL>){
	chomp;
	if(/^>(.*)/){
		$seq_name=$1;
		$total++;
		$NO=(split(/_x/,$_))[1];
		$total_collapse+=$NO;
	}
	else{
		$len=length($_);
		$hash{$seq_name}=$len;
		$colla{$seq_name}=$_;
	}
}
close(COL);
###################################################################################

my ($miRNA,$rRNA,$scRNA,$snRNA)=(0,0,0,0);
my ($snoRNA,$tRNA,$other)=(0,0,0);
my %query;

my ($miRNA_collapse, $rRNA_collapse, $scRNA_collapse,$snRNA_collapse)=(0,0,0,0);
my ($snoRNA_collapse, $tRNA_collapse,$other_collapse)=(0,0,0);

open IN,"<$blastFile" or die "can't open the reads file.$!\n";
my @res_blast;
my $NO1;
my %miRNA_hash;
my %family;

while(<IN>){
	chomp;
	s/\s+/\t/g;
	@res_blast=split(/\t/,$_);

	my ($Q_name,$S_name,$identity,$Q_start,$Q_end)=@res_blast[0,1,2,6,7];
	$NO1=(split(/_x/,$Q_name))[1];
	if(exists $hash{$Q_name}){
		my $Q_length=$hash{$Q_name};
		my $coverage=sprintf("%.2f",($Q_end-$Q_start+1)/$Q_length);
		if($coverage>=$my_cov && $identity>=$my_iden && !(exists $query{$Q_name})){
#	print OUT $_,"\t",$coverage,"\n";	
			$query{$Q_name}=$NO1;
			my $Q_seq=$colla{$Q_name};
			if($S_name=~/mir|let/six){
				$miRNA++;
				$miRNA_collapse+=$NO1;
				my $mir=(split(/;/,$S_name))[1];
				if(exists $family{$mir}){$family{$mir}=$family{$mir}+$NO1;}
				else{$family{$mir}=$NO1;}
				$mir=$mir."\t".$Q_name;
				$miRNA_hash{$mir}=$Q_seq;
				delete $query{$Q_name};
			}
			elsif($S_name=~/.*rRNA/six){
				$rRNA++;
				$rRNA_collapse+=$NO1;
			}elsif($S_name=~/.*?scRNA/six){
				$scRNA++;
				$scRNA_collapse+=$NO1;
			}elsif($S_name=~/.*?snRNA/six){
				$snRNA++;
				$snRNA_collapse+=$NO1;
			}elsif($S_name=~/.*?snoRNA/six){
				$snoRNA++;
				$snoRNA_collapse+=$NO1;
			}elsif($S_name=~/.*?tRNA/six){
				$tRNA++;
				$tRNA_collapse+=$NO1;
			}else{
				$other++;
				$other_collapse+=$NO1;
			}
		}
	}
	else{print "error,error,error!!!";}
}
close(IN);

open(MIR,">miRNA_predict_$id.fa") or die("Can not open the file $!");
for my $mir (sort keys %miRNA_hash){print MIR ">",$mir,"\n",$miRNA_hash{$mir},"\n";}
close (MIR);

open(FAM,">family_$id.txt") or die("Can not open the file $!");
for my $fam (sort{$family{$b}<=>$family{$a}} keys %family){print FAM $fam,"\t",$family{$fam},"\n";}
close (FAM);


my ($blast_total,$blast_collapse,$blast_uniq_rate,$blast_collapse_rate)=(0,0,0,0);

$blast_total=$miRNA+$rRNA+$scRNA+$snRNA+$snoRNA+$tRNA+$other;

$blast_collapse=$miRNA_collapse+$rRNA_collapse+$scRNA_collapse+$snRNA_collapse+$snoRNA_collapse+$tRNA_collapse+$other_collapse;

$blast_uniq_rate=sprintf("%.2f",$blast_total/$total*100);
$blast_collapse_rate=sprintf("%.2f",$blast_collapse/$total_collapse*100);

open(OUT,">$outputFile") or die("Can not open the file $!");

print  OUT "sample\tmiRNA\trRNA\tscRNA\tsnRNA\tsnoRNA\ttRNA\tother\tblast_total\ttotal\tblast_rate(%)\n";
print OUT $id,"_uniq\t$miRNA\t$rRNA\t$scRNA\t$snRNA\t$snoRNA\t$tRNA\t$other\t$blast_total\t$total($blast_uniq_rate)\n";
print OUT $id,"_collapsed\t$miRNA\t$rRNA_collapse\t$scRNA_collapse\t$snRNA_collapse\t$snoRNA_collapse\t$tRNA_collapse\t$other_collapse\t$blast_collapse\t$total_collapse($blast_collapse_rate)\n";


my $blastDB;
if($blastFile=~/Rfam/){$blastDB="Rfam";}
else{$blastDB="ncRNA";}
my $Rscripts=qq|
rgs<-commandArgs(T)
infile=rgs[1]
DB=rgs[2]
data<-read.delim2(infile,header=T,sep="\t",row.name=1)

data<-as.matrix(data)
	par(mar=c(4,3,5,3)+0.1)
	nam<-rownames(data)
	for(i in 1:nrow(data)){
		pdf(paste("picture",paste("/",DB,"_",nam[i],".pdf",sep=""),sep=""))
			a<-data[i,1:7]
			a1<-as.numeric(a[!(is.na(a))])
			c<-which(!(is.na(a)))
			labels=sprintf("%s =%3.1f%s",colnames(data)[c],100*a1/sum(a1),"%")
			yanse=rainbow(7)
			pie(as.numeric(as.vector(a1)),cex=0.8,ylim=c(0,1.1),xlim=c(0,1.1),labels=NA,border="white",col=yanse[c],radius=0.9,main=nam[i],cex.main=2)
			legend(0.5,1.1,legend=labels,cex=0.7,bty="n",fill=yanse[c])
			dev.off()
	}
|;

open(RSRC,">$blastDB.R")||die "Can not write R scripts\n";
print RSRC $Rscripts;
close RSRC;
system("Rscript $blastDB.R $outputFile $blastDB");


if (not -f "../$blastDB.xls"){`echo "sample\tmiRNA\trRNA\tscRNA\tsnRNA\tsnoRNA\ttRNA\tother\tblast_uniq\ttotal\tblast_uniq_rate" >../${blastDB}.xls`;}
`echo "${id}_uniq\t$miRNA\t$rRNA\t$scRNA\t$snRNA\t$snoRNA\t$tRNA\t$other\t$blast_total\t$total\t$blast_uniq_rate" >>../${blastDB}.xls`;
`echo "${id}_collapsed\t$miRNA\t$rRNA_collapse\t$scRNA_collapse\t$snRNA_collapse\t$snoRNA_collapse\t$tRNA_collapse\t$other_collapse\t$blast_collapse\t$total_collapse\t$blast_collapse_rate" >>../${blastDB}.xls`;

my $out="reads_".$id."_".$blastDB;
open (READ,">$out") or die ("Can not open the file $!");

for my $key(sort keys %query){
	print READ $key,"\n";
#	print READ ">",$key,"\n";
#	print READ $colla{$key},"\n";
}
close OUT;
close READ;

