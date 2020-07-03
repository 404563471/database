argv <- commandArgs(TRUE);
strinput=readLines(as.character(argv[1]));#input groups file

#strinput<-"25.6,22.2,28.0,29.8,23.3,45.3;24.4,33,22,35,45,65,30.0,29.0,27.5;25.0,27.7,23.0,32.2;28.8,,45,45,45,45,28.0,31.5,25.9;20.6,21.2,22.0,21.2"
#strinput<-readLines("data.txt") 
strarrays<-unlist(strsplit(strinput,";"))
#判断分组个数，若等于二，做t检验；若大于二，做方差分析
if(length(strarrays)==2){
  array1 <- as.numeric(unlist(strsplit(strarrays[1],",")))
  array2 <- as.numeric(unlist(strsplit(strarrays[2],",")))
  ttest <- t.test(array1,array2,alternative = "two.sided",  paired = FALSE)
  p_value <- ttest$p.value
}else if (length(strarrays)>2){
  strdata<-as.numeric(unlist(strsplit(strarrays,",")))
  strtype <- numeric()
  for(i in 1:length(strarrays)){
    k<-length(unlist(strsplit(strarrays[i],",")))
    strtype <- c(strtype,k)
  }
  anovadata <-data.frame(expression=strdata,subtype=factor(rep(c(1:length(strarrays)),strtype)))
  #检查数据是否满足方差齐性
  bartlett_test <- bartlett.test(expression ~ subtype, data= anovadata)
  if (bartlett_test$p.value > 0.05){
    anova_test<-aov(expression~subtype,data=anovadata)
    p_value <- summary(anova_test)[[1]]["Pr(>F)"][[1]][1]
  }else if (bartlett_test$p.value <= 0.05){
    kruskal_test <- kruskal.test(expression~subtype,data = anovadata)
    p_value <- kruskal_test$p.value
  }
}

cat(p_value)

