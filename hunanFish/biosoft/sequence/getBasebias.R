rgs<-commandArgs(T)
id=rgs[1]
root=rgs[2]
pdf(paste(root,"/AllBaseBias_",id,".pdf",sep=""))

allbasebias<-read.delim2(paste("AllBaseBias_",id,".txt",sep=""),sep="\t",header=T,row.names=1)
data<-t(as.matrix(allbasebias[1:4]))
for(j in 1:ncol(data)){data[,j]<-100*data[,j]/sum(data[,j])}
barplot(data,width=1,col=rainbow(4),ylim=c(0,100),xlim=c(0,ncol(data)*1.5+5),space=0.5,bty="n",xpd=TRUE,xlab="length(nt)",ylab="percent(%)",main="miRNA All Base Bias")
a<-names(allbasebias)
legend("topright",lwd=5,legend=a,cex=0.8,lty=1,col=rainbow(4),text.font=2,text.width=1)
dev.off()

######################################################


pdf(paste(root,"/FirstBaseBias",id,".pdf",sep=""))
firstbasebias<-read.delim2(paste("FirstBaseBias_",id,".txt",sep=""),sep="\t",header=T,row.names=1)
data<-t(as.matrix(firstbasebias[,1:4]))
for(j in 1:ncol(data)){
	data[is.na(data[,j]),j]=0
	data[,j]<-100*data[,j]/sum(data[,j])
	}
barplot(data,width=1,col=rainbow(4),space=0.5,xlim=c(0,15),bty="n",xpd=TRUE,xlab="length(nt)",ylab="percent(%)",main="miRNA First Base Bias")
a<-names(firstbasebias)
legend("topright",lwd=5,legend=a[2:5],cex=0.8,lty=1,col=rainbow(4),text.font=2,text.width=1)
dev.off()

