argv <- commandArgs(TRUE);
gene=as.character(argv[1]);#gene name
path_output=as.character(argv[2]);#path of output data 
survival_type=as.character(argv[3]);#choose survival type  "OS" or "DFS" and so on
cnv_flag=as.character(argv[4]); # "cnv" or other

setwd(path_output); # set path

#path_output="./"
#gene="BRCA2"
#survival_type="OS"

if(cnv_flag=='cnv')
{
  cnv_deletion_matrix <- as.matrix(read.table(paste(path_output,"cnv_deletion_matrix.txt",sep = "")))
  cnv_amplificatioin_matrix <- as.matrix(read.table(paste(path_output,"cnv_amplification_matrix.txt",sep = "")))
  cnv_deletion_time <- as.matrix(as.numeric(cnv_deletion_matrix[,2]))
  cnv_deletion_event <- as.matrix(as.numeric(cnv_deletion_matrix[,3]))
  cnv_amplificatioin_time <- as.matrix(as.numeric(cnv_amplificatioin_matrix[,2]))
  cnv_amplificatioin_event <- as.matrix(as.numeric(cnv_amplificatioin_matrix[,3]))

download <- matrix(nrow=max(nrow(cnv_amplificatioin_time),nrow(cnv_deletion_time)),ncol=4)

download[1:nrow(cnv_deletion_time),1] <- cnv_deletion_time
download[1:nrow(cnv_deletion_time),2] <- cnv_deletion_event
download[1:nrow(cnv_amplificatioin_time),3] <- cnv_amplificatioin_time
download[1:nrow(cnv_amplificatioin_time),4] <- cnv_amplificatioin_event
colnames(download)<-c(paste(survival_type," time of cnv deletion group",sep = ""), paste(survival_type," event of cnv deletion group",sep = ""), paste(survival_type," time of cnv amplification group",sep = ""), paste(survival_type," event of cnv amplfication group",sep = ""))
filename<- paste(gene,"_CNV_",survival_type,".csv",sep="")

write.csv(download, paste(path_output, filename, sep=""), row.names=FALSE)

}else{
low_matrix <- as.matrix(read.table(paste(path_output,"low_matrix.txt",sep = "")))
high_matrix <- as.matrix(read.table(paste(path_output,"high_matrix.txt",sep = "")))
low_time <- as.matrix(as.numeric(low_matrix[,2]))
low_event <- as.matrix(as.numeric(low_matrix[,3]))
high_time <- as.matrix(as.numeric(high_matrix[,2]))
high_event <- as.matrix(as.numeric(high_matrix[,3]))

download <- matrix(nrow=max(nrow(high_time),nrow(low_time)),ncol=4)

download[1:nrow(low_time),1] <- low_time
download[1:nrow(low_time),2] <- low_event
download[1:nrow(high_time),3] <- high_time
download[1:nrow(high_time),4] <- high_event
colnames(download)<-c(paste(survival_type," time of low expression group",sep = ""), paste(survival_type," event of low expression group",sep = ""), paste(survival_type," time of high expression group",sep = ""), paste(survival_type," event of high expression group",sep = ""))
filename<- paste(gene,"_",survival_type,".csv",sep="")

write.csv(download, paste(path_output, filename, sep=""), row.names=FALSE)
}
