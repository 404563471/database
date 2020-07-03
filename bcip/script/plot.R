argv <- commandArgs(TRUE);
gene=as.character(argv[1]);#gene name
#path_input_data=as.character(argv[2]);#path of input data
path_output=as.character(argv[2]);#path of output data 
path_r_file=as.character(argv[3]);
analysis_type=as.character(argv[4]);#choose analysis type

sqlcmd=as.character(argv[5]);  #查询命令
cutoff.survival.value=as.character(argv[6]);  #输入cutoff值，如果是optimal，则是首次计算，需要计算optimal值�?

#gene='DDR1';
#path_input_data='input/';
#path_r_file='D:/Program Files/R_workspace/yingxm_demo/'
#path_output='output/'
#analysis_type='rfs_survival';

setwd(path_r_file); # set path


#######数据库查询对�?#########
if (grepl("survival",analysis_type))
{
  library(survival)
  #定义函数，找出每个探针在所选样本中表达量的分界值（低和高，然后做生存分析，观察是否可将生存曲线分开�?
  cutoff.survival <- function(marker, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95) {
  cat(paste("Optimizing cutoff using method survival", method, sep="_"))  #类似于print，输出，paste用于拼接，sep定义连接�?
  index <- intersect(which(!is.na(marker)), intersect(which(!is.na(time)), which(!is.na(event))))
  marker <- marker[index]
  time <- time[index]
  event <- event[index]
  n <- length(marker)
  nmarker <- length(unique(marker))
  if (nmarker < 3) stop("insufficient data")
  index <- order(marker)
  marker <- marker[index]
  time <- time[index]
  event <- event[index]
  y <- Surv(time, event) 
  q <- 1-(1-conf.level/100)/2
  z <- qnorm(q) #q分位数函数（标准正态分布）
  #nlow <- nmin:(n-nmin)
  nlow<-(nmin+1):(n-nmin)
  Y <- matrix(nrow=length(nlow), ncol=2)
  colnames(Y) <- c("expression","p.value") 
  rownames(Y) <- nlow
  for (i in nlow) {
    cat(".")
    j <- i-nmin
    if (marker[i] != marker[i+1]) {
      x <- c(rep(0, i), rep(1, n-i))
      model <- summary(coxph(y ~ x))
      Y[j, "expression"] <- marker[i]
      Y[j, "p.value"] <- model[[surv.test]]["pvalue"]
    }
    else Y[j, "p.value"] <- 2
  }
  if (method == "significance") {
    index.optimal <- which.min(as.numeric(Y[, "p.value"]))
    rownames(Y)[index.optimal] <- "optimal"
  }
  cat("\n")
  return(Y)
}

  #path_input_data="./"
  #gene="RFC2"
  #analysis_type="os_survival"

  library(RMySQL)

  #sqlcmd="select expression_metabric.EXPRESSION_VALUE, sample_metabric.OS_SURVIVAL_EVENT, sample_metabric.OS_SURVIVAL_TIME from expression_metabric, sample_metabric where (expression_metabric.GENE_NAME='RFC2' and 1=1 and expression_metabric.SAMPLE_ID = sample_metabric.ID and sample_metabric.N_C_TYPES = 'C')"
  #######R的数据库连接，根据分组条件筛选出样本id并得到time、event、exprs_data，存到变量query_data�?#########
  conn <- dbConnect(MySQL(), dbname="bcdb", username="root", password="***")
  query_data = dbGetQuery(conn, sqlcmd)
  dbDisconnect(conn)

  #time1 <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_time.txt',sep=""),header=FALSE,sep="\t"));
  #event1 <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_event.txt',sep=""),header=FALSE,sep="\t"));
  #exprs_data1 <- as.matrix(read.table(paste(path_input_data, gene, '_expression.txt',sep=""),header=FALSE,sep="\t"));
  #time <- c(as.factor(as.matrix(query_data[,3],stringsAsFactors = FALSE)))
  #event <- c(as.factor(as.matrix(query_data[,2],stringsAsFactors = FALSE)))
  #exprs_data <- c(as.factor(as.matrix(query_data[,1],stringsAsFactors = FALSE)))
  time <- as.matrix(as.numeric(query_data[,3]))
  event <- as.matrix(as.numeric(query_data[,2]))
  exprs_data <- as.matrix(as.numeric(query_data[,1]))
  if(cutoff.survival.value=='optimal')
  {
    Y <- cutoff.survival(exprs_data, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95)
    cutoff.survival.value <- Y[which(row.names(Y)=="optimal"),1]
    file.create(paste(path_output,"params.txt",sep = ""))
    cat(sqlcmd, file=paste(path_output,"params.txt",sep = ""))
    cat("\n", file=paste(path_output,"params.txt",sep = ""), append=TRUE)
    cat(quantile(exprs_data, 1, na.rm=TRUE), file=paste(path_output,"params.txt",sep = ""), append=TRUE)
    cat("\n", file=paste(path_output,"params.txt",sep = ""), append=TRUE)
    cat(quantile(exprs_data, 0, na.rm=TRUE), file=paste(path_output,"params.txt",sep = ""), append=TRUE)
    cat("\n", file=paste(path_output,"params.txt",sep = ""), append=TRUE)
    cat(cutoff.survival.value, file=paste(path_output,"params.txt",sep = ""), append=TRUE)
  }
  exprs_data_cutoff <- matrix(nrow=nrow(exprs_data),ncol=1)
  if(length(which(as.numeric(exprs_data[,1]) < as.numeric(cutoff.survival.value)))<10  | length(which(as.numeric(exprs_data[,1]) >= as.numeric(cutoff.survival.value)))<10){
    cat("Warning: The sample number of an expression group is less than 10.", file=paste(path_output,"err.log",sep = ""))
  }else{
  exprs_data_cutoff[which(as.numeric(exprs_data[,1]) <= as.numeric(cutoff.survival.value)),1] <- 0
  exprs_data_cutoff[which(as.numeric(exprs_data[,1]) > as.numeric(cutoff.survival.value)),1] <- 1


  #for download
  survival_type=toupper(unlist(strsplit(analysis_type,"_"))[1])

  down_low_time <- as.matrix(time[which(as.numeric(exprs_data[,1]) <= as.numeric(cutoff.survival.value)),1])
  down_low_event <- as.matrix(event[which(as.numeric(exprs_data[,1]) <= as.numeric(cutoff.survival.value)),1])
  down_high_time <- as.matrix(time[which(as.numeric(exprs_data[,1]) > as.numeric(cutoff.survival.value)),1])
  down_high_event <- as.matrix(event[which(as.numeric(exprs_data[,1]) > as.numeric(cutoff.survival.value)),1])

  download <- matrix(nrow=max(nrow(down_high_time),nrow(down_low_time)),ncol=4)
  download[1:nrow(down_low_time),1] <- down_low_time
  download[1:nrow(down_low_time),2] <- down_low_event
  download[1:nrow(down_high_time),3] <- down_high_time
  download[1:nrow(down_high_time),4] <- down_high_event
  colnames(download)<-c(paste(survival_type," time of low expression group",sep = ""), paste(survival_type," event of low expression group",sep = ""), paste(survival_type," time of high expression group",sep = ""), paste(survival_type," event of high expression group",sep = ""))
  filename<- paste(gene,"_",survival_type,".csv",sep="")
  write.csv(download, paste(path_output, filename, sep=""), row.names=FALSE)





  kaplan_meier <- survfit(Surv(time,event) ~ exprs_data_cutoff,conf.type = "log") #conf.type="log-log"?
  low_matrix <- matrix(nrow=kaplan_meier$strata[1],ncol=4)
  high_matrix <- matrix(nrow=kaplan_meier$strata[2],ncol=4)
  colnames(low_matrix) <- c("surv", "time", "event", "censor")
  colnames(high_matrix) <- c("surv", "time", "event", "censor")
  i <- kaplan_meier$strata[1]
  low_matrix[,1] <- kaplan_meier$surv[1:i]
  low_matrix[,2] <- kaplan_meier$time[1:i]
  low_matrix[,3] <- kaplan_meier$n.event[1:i]
  low_matrix[,4] <- kaplan_meier$n.censor[1:i]
  i <- i+1
  high_matrix[,1] <- kaplan_meier$surv[i:length(kaplan_meier$surv)]
  high_matrix[,2] <- kaplan_meier$time[i:length(kaplan_meier$surv)]
  high_matrix[,3] <- kaplan_meier$n.event[i:length(kaplan_meier$surv)]
  high_matrix[,4] <- kaplan_meier$n.censor[i:length(kaplan_meier$surv)]
  write.table(low_matrix, paste(path_output,"low_matrix.txt",sep = ""))
  write.table(high_matrix, paste(path_output,"high_matrix.txt",sep = ""))
  survival_median <- c(summary(kaplan_meier)$table[,"median"][[1]],summary(kaplan_meier)$table[,"median"][[2]])
  names(survival_median) <- c("low","high")
  write.table(as.matrix(survival_median), paste(path_output,"survival_median.txt",sep = ""), col.names=FALSE, row.names=FALSE)

  mycoxph <- coxph(Surv(time,event)~as.factor(exprs_data_cutoff)) ###cnv_data_cutoff对应CNV的生存，转录组生存相应的变为其cutoff向量即可
  tmp <- summary(mycoxph)
  p_value <- tmp$sctest[[3]]
  HR <- tmp$coefficients[,"exp(coef)"]
  file.create(paste(path_output,"p_HR.txt",sep = ""))
  cat(p_value, file=paste(path_output,"p_HR.txt",sep = ""))      #p-value
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(HR, file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #HR
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(kaplan_meier$n[1], file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #low
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)
  cat(kaplan_meier$n[2], file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)      #high
  cat("\n", file=paste(path_output,"p_HR.txt",sep = ""), append=TRUE)

#  summary(kaplan_meier)$table[,"median"][[1]][1]
#  summary(kaplan_meier)$table[,"median"][[1]][2]
  }
}
