
#install.packages('plyr')
library(plyr)

base_path<-'/Users/chens22/Documents/miRNA/'

tumors<-c('ACC', 'LGG', 'COAD', 'BLCA', 
          'CESC', 'GBM', 'DLBC', 'KICH', 
          'THYM', 'TGCT', 'HNSC','READ', 
          'ESCA', 'BRCA', 'CHOL', 'KIRC',
          'KIRP', 'LAML', 'LIHC', 'LUAD',
          'LUSC', 'MESO', 'OV', 'PAAD',
          'PCPG', 'PRAD', 'SARC',
          'SKCM', 'STAD', 'THCA', 'UCEC',
          'UCS', 'UVM')

manifest_data<-data.frame(NAME=NA, DISEASE.ABBV=NA)
abbv_data<-read.csv(file=paste(base_path, 'disease-abbreviations.csv', sep=''), header=TRUE)

for(t in tumors){
  
  curr_file<-paste(base_path, t, '/', t, '-manifest.csv', sep='')
  curr_data<-read.csv(file=curr_file, stringsAsFactors=FALSE, header=TRUE)
  
  names(curr_data)<-toupper(names(curr_data))
  names(curr_data)<-gsub('_', '.', names(curr_data))
  
  for(i in 1:nrow(curr_data)){
    
    curr_name<-toString(factor(curr_data[i, 'NAME']))
    curr_name<-strsplit(curr_name, '.', fixed=TRUE)[[1]][1]
    curr_data[i, 'NAME']<-curr_name
  
    
  }
  curr_data[,'DISEASE.ABBV']<-t
  
  manifest_data<-rbind.fill(manifest_data, curr_data)
  
}
manifest_data<-unique(manifest_data)
manifest_file<-paste(base_path, 'all-tumor-manifest.csv', sep='')
write.table(manifest_data, file=manifest_file, sep=',', row.names=FALSE)
