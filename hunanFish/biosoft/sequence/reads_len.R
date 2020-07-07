rgs<-commandArgs(T)
id=rgs[1]
root=rgs[2]
id
root
paste(root,"/reads_length_",id,".pdf",sep="")
pdf(paste(root,"/reads_length_",id,".pdf",sep=""))
reads_len<-read.delim2(paste("length_",id,".txt",sep=""),sep="\t",header=F,row.name=1)
rNo<-nrow(reads_len)-2
data<-as.numeric(as.vector(reads_len[1:rNo,]))
data1<-data/sum(data)*100
data1<-round(data1,2)
#names(data1)<-rownames(reads_len)[1:rNo]

barplot(data1,col=rainbow(rNo),ylim=c(0,max(data1)*1.2),xlim=c(0,1*rNo),width=1,space=0,tcl=-0.1,ylab="the percent of the reads(%)",xlab="the length of reads",main=id,las=1,cex.main=1.5)
text(x=seq(0,(rNo-1),by=1),y=data1+0.5,adj=-0.1,labels=data1,cex=0.5,font=1)
axis(1,0.5:(rNo-0.5),rownames(reads_len)[1:rNo],cex.axis=0.5)
abline(h=0)
dev.off()

