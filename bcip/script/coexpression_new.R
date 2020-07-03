#######2015-10-20修改，coexpression查询优化脚本，采用直接load表达的RData文件的方�?##########################################
#######运行命令示例：Rscript 【R脚本�? 【基因名�? 【工作目录�? 【分析类型�? 【sampleID查询条件�? 【DatasetID�?
##Rscript /home/yxm/coexpression_data/coexpression_new.R DDR1 /home/yxm/coexpression_data coexpression "select sample_metabric.ID from sample_metabric where (sample_metabric.N_C_TYPES='C' and sample_metabric.tnbc_subtypes='non_TNBC' and sample_metabric.pam50_subtypes='Basal')" metabric
#'select sample_gse1456_gpl96.ID from sample_gse1456_gpl96 where sample_gse1456_gpl96.tnbc_subtypes='non_TNBC' and sample_gse1456_gpl96.pam50_subtypes='Her2' and sample_gse1456_gpl96.N_C_TYPES = 'C'

argv <- commandArgs(TRUE);
gene=as.character(argv[1]);#gene name
path_r_file=as.character(argv[2]);
analysis_type=as.character(argv[3]);#choose analysis type
sqlcmd=as.character(argv[4]);  #sampleID查询条件
#dataset_RData=as.character(argv[5]);   #该数据集的表达数据R文件路径
dataset_RData=paste('/var/www/bcancer/Public/coexpression_data/',as.character(argv[5]),'.RData',sep='');   #该数据集的表达数据R文件路径
prognosis=as.character(argv[6])

#gene='DDR1';
#path_r_file='/home/yxm/coexpression_test'
#analysis_type='coexpression';
#sqlcmd="select sample_metabric.ID from sample_metabric where (sample_metabric.N_C_TYPES='C' and sample_metabric.tnbc_subtypes='non_TNBC' and sample_metabric.pam50_subtypes='Basal')"
#dataset_RData="/home/yxm/coexpression_data/metabric.RData";

#sqlcmd="select sample_metabric.ID,OS_SURVIVAL_TIME from sample_metabric where 1=1 and sample_metabric.N_C_TYPES = 'C'"

setwd(path_r_file); # set path

if (analysis_type == 'coexpression')
{
library(RMySQL)
library(WGCNA)
library(reshape2)

#######R的数据库连接，根据分组条件（non_TNBC且Basal）筛选出样本id（查询得�?80个样本符合条件），存到变量sample_ids#########
conn <- dbConnect(MySQL(), dbname="bcdb", username="root", password="***")
sample_ids = dbGetQuery(conn, sqlcmd)

if(prognosis=='good')
{
  time_median <- quantile(as.matrix(as.numeric(sample_ids[,2])), 0.5, na.rm=TRUE)
  sqlcmd=paste(sqlcmd, "and", unlist(strsplit(unlist( strsplit(sqlcmd, ","))[2], " "))[1], ">=", time_median, sep=" ")
  sample_ids = dbGetQuery(conn, sqlcmd)
}
if(prognosis=='poor')
{
  time_median <- quantile(as.matrix(as.numeric(sample_ids[,2])), 0.5, na.rm=TRUE)
  sqlcmd=paste(sqlcmd, "and", unlist(strsplit(unlist( strsplit(sqlcmd, ","))[2], " "))[1], "<", time_median, sep=" ")
  sample_ids = dbGetQuery(conn, sqlcmd)
}

dbDisconnect(conn)

if (length(sample_ids$ID)<20){
  write("Warning: Samples less than 20!","err.log")
}else{
  load(dataset_RData)
  exprs_data <- tumor_data[,sample_ids$ID]
  if(sum(is.na(exprs_data[gene,])) > 0.5*length(sample_ids$ID)){
    write("Warning: Missing values of the query gene more than half!","err.log")
  }else{
    gene_list <- rownames(exprs_data)
    corrcoef <- cor(exprs_data[gene,],t(exprs_data),use ="pairwise.complete.obs", method = "pearson",quick=1)
    p.value <- corPvalueFisher(corrcoef,nSamples = ncol(exprs_data),twoSided = T) 
    q.value <- p.adjust(p.value,method = "fdr",n = length(p.value))
    cor_test <- matrix(nrow=length(gene_list),ncol=3)
    rownames(cor_test)<- colnames(corrcoef)
    colnames(cor_test)<-c("gene","corrcoef","q_value")
    cor_test[,1] <- rep(gene,length(gene_list))
    cor_test[,2] <- corrcoef
    cor_test[,3] <- q.value
    tmp_genes <- names(which(rowSums(is.na(exprs_data)) > 0.5*length(sample_ids$ID)))
    cor_test[match(tmp_genes,rownames(cor_test)),2:3] <- NA
    coexpression_data <- cor_test[which(abs(as.numeric(cor_test[,2])) > 0.3 & as.numeric(cor_test[,3]) < 0.05),]
    coexpression_data <- coexpression_data[-which(rownames(coexpression_data)== gene),]
    cor_data <- coexpression_data[order(abs(as.numeric(coexpression_data[,2])),decreasing = T),]
    write.table(cor_data, paste(gene,"_coexpression.txt",sep = ""),col.names=FALSE)
    
    #for download

    download <- cbind(rownames(cor_data),cor_data)
    dim(download)
    colnames(download)<-c("Gene1", "Gene2", "Corrcoef", "Qvalue")
    filename<- paste(gene,"_coexpression",".csv",sep="")
    write.csv(download, paste(filename, sep=""), row.names=FALSE)
  }
}

}
