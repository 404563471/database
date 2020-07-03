argv <- commandArgs(TRUE);
gene=as.character(argv[1]);#gene name
path_input_data=as.character(argv[2]);#path of input data
path_output=as.character(argv[3]);#path of output data 
path_r_file=as.character(argv[4]);
analysis_type=as.character(argv[5]);#choose analysis type

#gene='DDR1';
#path_input_data='input/';
#path_r_file='D:/Program Files/R_workspace/yingxm_demo/'
#path_output='output/'
#analysis_type='rfs_survival_cnv';

setwd(path_r_file); # set path
if (grepl("survival_cnv",analysis_type))
{
  library(survival)
  time <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_time.txt',sep=""),header=FALSE,sep="\t"));
  event <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_event.txt',sep=""),header=FALSE,sep="\t"));
  cnv_data <- as.matrix(read.table(paste(path_input_data, gene, '_cnv.txt',sep=""),header=FALSE,sep="\t"));
  #exprs_data <- as.matrix(read.table(paste(path_input_data, gene, '_expression.txt',sep=""),header=FALSE,sep="\t"));
  #Y <- cutoff.survival(exprs_data, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95)
  #cutoff.survival.value <- Y[which(row.names(Y)=="optimal"),1]
  cnv_data_cutoff <- matrix(nrow=nrow(cnv_data),ncol=1)
  cnv_data_cutoff[which(as.numeric(cnv_data[,1]) <= -0.1),1] <- 0
  cnv_data_cutoff[which(as.numeric(cnv_data[,1]) >= 0.1),1] <- 1
  
  if(length(which(as.numeric(cnv_data[,1]) <= -0.1))<10  | length(which(as.numeric(cnv_data[,1]) >= 0.1))<10){
    cat("Warning: The sample number of an cnv data group is less than 10.", file=paste(path_output,"err.log",sep = ""))
  }else{
  
  #for download
  survival_type=toupper(unlist(strsplit(analysis_type,"_"))[1])

  down_cnv_deletion_time <- as.matrix(time[which(as.numeric(cnv_data[,1]) <= -0.1),1])
  down_cnv_deletion_event <- as.matrix(event[which(as.numeric(cnv_data[,1]) <= -0.1),1])
  down_cnv_amplificatioin_time <- as.matrix(time[which(as.numeric(cnv_data[,1]) >= 0.1),1])
  down_cnv_amplificatioin_event <- as.matrix(event[which(as.numeric(cnv_data[,1]) >= 0.1),1])

  download <- matrix(nrow=max(nrow(down_cnv_amplificatioin_time),nrow(down_cnv_deletion_time)),ncol=4)
  download[1:nrow(down_cnv_deletion_time),1] <- down_cnv_deletion_time
  download[1:nrow(down_cnv_deletion_time),2] <- down_cnv_deletion_event
  download[1:nrow(down_cnv_amplificatioin_time),3] <- down_cnv_amplificatioin_time
  download[1:nrow(down_cnv_amplificatioin_time),4] <- down_cnv_amplificatioin_event
  colnames(download)<-c(paste(survival_type," time of cnv deletion group",sep = ""), paste(survival_type," event of cnv deletion group",sep = ""), paste(survival_type," time of cnv amplification group",sep = ""), paste(survival_type," event of cnv amplfication group",sep = ""))
  filename<- paste(gene,"_CNV_",survival_type,".csv",sep="")
  write.csv(download, paste(path_output, filename, sep=""), row.names=FALSE)



  kaplan_meier <- survfit(Surv(time,event) ~ cnv_data_cutoff,conf.type = "log") #conf.type="log-log"?
  cnv_deletion_matrix <- matrix(nrow=kaplan_meier$strata[1],ncol=4)
  cnv_amplification_matrix <- matrix(nrow=kaplan_meier$strata[2],ncol=4)
  colnames(cnv_deletion_matrix) <- c("surv", "time", "event", "censor")
  colnames(cnv_amplification_matrix) <- c("surv", "time", "event", "censor")
  i <- kaplan_meier$strata[1]
  cnv_deletion_matrix[,1] <- kaplan_meier$surv[1:i]
  cnv_deletion_matrix[,2] <- kaplan_meier$time[1:i]
  cnv_deletion_matrix[,3] <- kaplan_meier$n.event[1:i]
  cnv_deletion_matrix[,4] <- kaplan_meier$n.censor[1:i]
  i <- i+1
  cnv_amplification_matrix[,1] <- kaplan_meier$surv[i:length(kaplan_meier$surv)]
  cnv_amplification_matrix[,2] <- kaplan_meier$time[i:length(kaplan_meier$surv)]
  cnv_amplification_matrix[,3] <- kaplan_meier$n.event[i:length(kaplan_meier$surv)]
  cnv_amplification_matrix[,4] <- kaplan_meier$n.censor[i:length(kaplan_meier$surv)]
  write.table(cnv_deletion_matrix, paste(path_output,"cnv_deletion_matrix.txt",sep = ""))
  write.table(cnv_amplification_matrix, paste(path_output,"cnv_amplification_matrix.txt",sep = ""))
  survival_median <- c(summary(kaplan_meier)$table[,"median"][[1]],summary(kaplan_meier)$table[,"median"][[2]])
  names(survival_median) <- c("cnv_deletion","cnv_amplification")
  write.table(as.matrix(survival_median), paste(path_output,"survival_median.txt",sep = ""), col.names=FALSE, row.names=FALSE)
  
  mycoxph <- coxph(Surv(time,event)~as.factor(cnv_data_cutoff)) ###cnv_data_cutoff对应CNV的生存，转录组生存相应的变为其cutoff向量即可
  tmp <- summary(mycoxph)
  p_value <- tmp$sctest[[3]]
  HR <- tmp$coefficients[,"exp(coef)"]
  file.create(paste(path_output,"p_HR.txt",sep = ""))
  cat(p_value, file=paste(path_output,"p_HR.txt",sep = ""))      #p-value
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(HR, file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #HR
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(kaplan_meier$n[1], file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #loss
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(kaplan_meier$n[2], file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #gain
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)

#  summary(kaplan_meier)$table[,"median"][[1]][1]
#  summary(kaplan_meier)$table[,"median"][[1]][2]
}
}
