#!/usr/bin/perl -w

use strict;

my $infile = $ARGV[0];
my $group = $ARGV[1];
$group = "\"".join("\",\"",split(",",$group))."\"";

my $Rscripts=qq|
expr<-read.table("$infile",header=T,row.names=1)
gene<-rownames(expr)
nrows <- nrow(expr)
group <- c($group)

if(length(group) != length(expr))
	stop("Length of 'group' must equal number of columns in 'expr'")

groupNumber <- length(unique(group))

library("ggplot2")
if(groupNumber < length(expr)){
	GE <- lapply(split(t(expr),group), FUN = function(u) matrix(u, nrow = nrows, byrow = TRUE))
	GS <- split(names(expr),group)
	for (s in names(GS)){
		if(length(group[group == s]) == 2){
			tmp <- as.data.frame(GE[[s]])
			names(tmp) <- c("sample1","sample2")
			sample1 <- GS[[s]][1]
			sample2 <- GS[[s]][2]
			fn <- paste(sample1,"AND",sample2,".pdf",sep="")
			pdf(file=fn);
			par(mar=par()\$mar+c(0,2,0,0))
			plot.new();
			c <- cor(tmp\$sample1,tmp\$sample2,method="spearman")
			myT <- paste("Biological Replicates (",s,")\\n",expression(R ^ 2 ), " = ",c, sep="")
			myXL <- paste("Sample", sample1, "log10(FPKM)")
			myYL <- paste("Sample", sample2, "log10(FPKM)")
			p <-ggplot(tmp, aes(x=log10(sample1+1),y=log10(sample2+1)))
			p <- p + geom_abline(intercept=0,slope=1,linetype=2) + geom_point(size=1.5,alpha=I(1/3))
			p <- p + stat_smooth(method="lm",fill="blue",alpha=0.2)
			p <- p + ylab(myYL) + xlab(myXL) + geom_rug(size=0.8,alpha=0.01) + labs(title=myT)
			print(p)
			dev.off();
		}
	}
}

ts<-stack(expr)
allE<-ts[ts\$values>1e-8,]
Sample <- factor(allE\$ind)
m<-ggplot(allE)
m <- m + geom_density(aes(x= log10(values),group=Sample,color=Sample,fill=Sample),alpha=I(1/3))
m <- m + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)

pdf(file="density.pdf");
par(mar=par()\$mar+c(0,2,0,0))
plot.new();
myT <- c("FPKM Density Distribution")
myY <- c("Density")
myX <- c("log10(FPKM)")
m + labs(title=myT) + xlab(myX) + ylab(myY)
dev.off()

pdf(file="boxplot.pdf");
par(mar=par()\$mar+c(0,2,0,0))
plot.new();
myT <- c("Boxplot of FPKM Distribution")
myY <- c("log10(FPKM)")
m <- ggplot(allE,aes(Sample,log10(values)))
m <- m + geom_boxplot(aes(fill=Sample))
m + labs(title=myT) + ylab(myY)
dev.off();
|;

	my $title = "tmp";
	open(RSRC,">$title.R")||die "Can not write R scripts\n";
	print RSRC $Rscripts;
	close RSRC;

	my $status = system("R CMD BATCH --no-save $title.R");

	if(($status >>=8) != 0){
		#croak "Couldn't run qsub to submit job ($@)";
		die "[ERROR] Couldn't run :R , please check it manually\n";
	}
	#system("rm $title.R*");
