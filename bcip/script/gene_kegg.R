argv <- commandArgs(TRUE);
gene=as.character(argv[1]);#gene Entrez_ID
pathwayid=as.character(argv[2]); #pathway id
path_r_file=as.character(argv[3]);

#gene="10488"
#pathwayid="path:hsa05030"
#path_r_file="/home/yxm/kegg"

setwd(path_r_file)
pathwayid<- substring(pathwayid,6)
#library(Cairo)
library(pathview)
pathview(gene.data=c(gene),gene.idtype="entrez", pathway.id = pathwayid, species = "hsa", out.suffix = gene)
