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
#analysis_type='rfs_survival';

setwd(path_r_file); # set path


#######数据库查询对应#########
#######time ---- sample表RFS_SURVIVAL_TIME列，按SAMPLE_ID排序
#######event ---- sample表RFS_SURVIVAL_EVENT列，按SAMPLE_ID排序
#######exprs_data ---- expression表该基因的所选样本的EXPRESSION_VALUE列，按SAMPLE_ID排序
#######查询后输出文件########
######输出文件命名方式：生存数据---“基因名+生存分析类型+数据”、基因表达数据---“基因名+数据”
######DDR1_survival_time.txt、DDR1_ds_survival_time.txt、DDR1_dfs_survival_time.txt、DDR1_rfs_survival_time.txt、DDR1_dmfs_survival_time.txt
######DDR1_survival_event.txt、DDR1_ds_survival_event.txt、DDR1_dfs_survival_event.txt、DDR1_rfs_survival_event.txt、DDR1_dmfs_survival_event.txt
######DDR1_expression.txt
if (grepl("survival",analysis_type))
{
  library(survival)
  #定义函数，找出每个探针在所选样本中表达量的分界值（低和高，然后做生存分析，观察是否可将生存曲线分开）
  cutoff.survival <- function(marker, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95) {
  cat(paste("Optimizing cutoff using method survival", method, sep="_"))  #类似于print，输出，paste用于拼接，sep定义连接符
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
  nlow <- nmin:(n-nmin)
  Y <- matrix(nrow=length(nlow), ncol=2)
  colnames(Y) <- c("expression","p.value") 
  rownames(Y) <- nlow
  for (i in nlow) {
    cat(".")
    j <- i-nmin+1   
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
  time <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_time.txt',sep=""),header=FALSE,sep="\t"));
  event <- as.matrix(read.table(paste(path_input_data, gene, '_', analysis_type, '_event.txt',sep=""),header=FALSE,sep="\t"));
  exprs_data <- as.matrix(read.table(paste(path_input_data, gene, '_expression.txt',sep=""),header=FALSE,sep="\t"));
  Y <- cutoff.survival(exprs_data, time, event, method="significance", thres.p=0.05, nmin=10, surv.test="sctest", endtime=max(time), conf.level=95)
  cutoff.survival.value <- Y[which(row.names(Y)=="optimal"),1]
  exprs_data_cutoff <- matrix(nrow=nrow(exprs_data),ncol=1)
  exprs_data_cutoff[which(as.numeric(exprs_data[,1]) < as.numeric(cutoff.survival.value)),1] <- 0
  exprs_data_cutoff[which(as.numeric(exprs_data[,1]) >= as.numeric(cutoff.survival.value)),1] <- 1
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
}

#add by 20150908输出共表达分析所要的数据
#######exprs_data ---- expression表所有基因的所选样本的EXPRESSION_VALUE列，按GENE_NAME和SAMPLE_ID排序
#######查询后输出文件########
######输出文件命名方式：表达数据---“表达数据”
######expression.txt#########
if (analysis_type == 'coexpression')
{
  sqlcmd <- as.character(argv[6]) # coexpressioin sql
  #sqlcmd <- "SELECT * FROM expression_gse11121_gpl96"
  
  #exprs <- read.table(paste(path_input_data, 'expression.txt',sep=""),header=TRUE,sep="\t");
 
  library(RMySQL)
  library(WGCNA)
  library(reshape2)
  conn <- dbConnect(MySQL(), dbname="bcdb1", username="root", password="***")
  exprs = dbGetQuery(conn, sqlcmd)
  dbDisconnect(conn)
  exprs_data <- acast(data = exprs, GENE_NAME ~ SAMPLE_ID, value.var = "EXPRESSION_VALUE")
  if (ncol(exprs_data) > 20){
    gene_list <- rownames(exprs_data)
    corrcoef <- cor(exprs_data[gene,],t(exprs_data),use ="pairwise.complete.obs", method = "pearson")
    p.value <- corPvalueFisher(corrcoef,nSamples = ncol(exprs_data),twoSided = T) 
    q.value <- p.adjust(p.value,method = "fdr",n = length(p.value))
    cor_test <- matrix(nrow=length(gene_list),ncol=3)
    rownames(cor_test)<- colnames(corrcoef)
    colnames(cor_test)<-c("gene","corrcoef","q_value")
    cor_test[,1] <- rep(gene,length(gene_list))
    cor_test[,2] <- corrcoef
    cor_test[,3] <- q.value
    coexpression_data <- cor_test[which(abs(as.numeric(cor_test[,2])) > 0.3 & as.numeric(cor_test[,3]) < 0.05),]
    if( class(coexpression_data)=="matrix"){
    coexpression_data <- coexpression_data[-which(rownames(coexpression_data)== gene),]
    cor_data <- coexpression_data[order(abs(as.numeric(coexpression_data[,2])),decreasing = T),]
    write.table(cor_data, paste(path_output,gene,"_coexpression.txt",sep = ""),col.names=FALSE)
    }else{
        #write("Warnning: no coexpression gene was found!","err.log")
        write("Warnning: no coexpression gene was found!",paste(path_output,"err.log",sep=""))
    }
  }
  else{
    write("Warnning: Samples less than 20!",paste(path_output,"err.log",sep=""))
  }
}