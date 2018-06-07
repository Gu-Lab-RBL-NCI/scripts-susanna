
## PURPOSE: format manifest file from CGC and add to all-tumor-manifest.csv
## INPUT: manifest data from CGC 		tumor-manifest.csv
## OUTPUT: formatted manifest data 		tumor-manifest.csv
## 		   all tumor manifest data 	all-tumor-manifest.csv

tumor <- 'ACC'

manifest_file<-paste0('/Users/chens22/Documents/', tumor, '/', tumor, '-manifest.csv')
manifest_data<-read.csv(file = manifest_file, header= T, stringsAsFactors=F)
names(manifest_data)<-toupper(names(manifest_data))
names(manifest_data)<-gsub('_', '.', names(manifest_data))
for(i in c(1,5,6,7,8,11,12,13,20)){
  manifest_data[,i]<-toupper(manifest_data[,i])
  manifest_data[,i]<-gsub('-', '.', manifest_data[,i], fixed=T)
  manifest_data[,i]<-gsub(' ', '.', manifest_data[,i], fixed=T)
  manifest_data[,i]<-gsub('_', '.', manifest_data[,i], fixed=T)
}
for(i in 1:ncol(manifest_data)){
  manifest_data[,i]<-toupper(manifest_data[,i])
}
for(n in unique(manifest_data$NAME)){
  new_name<-strsplit(toString(n), '.', fixed=T)[[1]][1]
  manifest_data[which(manifest_data$NAME==n), 'NAME']<-new_name
}
write.table(manifest_data, file=manifest_file, sep=',', row.names=F)

all_manifest_file<-'/Users/chens22/Documents/TCGA/all-tumor-manifest.csv'
all_manifest_data<-read.csv(file=all_manifest_file, header=T, stringsAsFactors=F)
disease_abbv<-all_manifest_data[, c('DISEASE.TYPE', 'DISEASE.ABBV')]
manifest_data<-merge(manifest_data, disease_abbv, by='DISEASE.TYPE', all=F, na.rm=T)
all_manifest_data<-rbind(manifest_data, all_manifest_data, by='')
