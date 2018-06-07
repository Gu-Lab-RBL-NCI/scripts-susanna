
## PURPOSE: get effective number of normal tissue
## INPUT: miRBase data                                      miRBase21-master.tsv
##        manifest data                                   all-tumor-manifest.tsv
##        isomiR summary data with calculated effective number  tumor.isomir.tsv
## OUTPUT: effective number of normal tissue for each miRNA    normal.patients.isomir.effective.tsv

require('data.table')
require('stringr')

base_path<-'/Users/chens22/Documents/miRNA/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=T)
manifest_file<-'/Users/chens22/Documents/miRNA/all-tumor-manifest.csv'
manifest_data<-read.csv(file=manifest_file, header=T)

tumors<-c('THYM', 'SKCM', 'READ', 'PCPG', 'CESC', 'PAAD', 'GBM', 'COAD', 
          'CHOL', 'ESCA', 'BLCA', 'KICH', 'UCEC', 'KIRP', 'HNSC', 'STAD', 'LUSC', 
          'LUAD', 'LIHC', 'PRAD', 'THCA')#BRCA, KIRC, 

getNormalPatients<-function(tumor){
  patients<-manifest_data[which(manifest_data$SAMPLE.TYPE=='SOLID.TISSUE.NORMAL' & manifest_data$DISEASE.ABBV==tumor), 'NAME']
  return(patients)
}

get5PMirna<-function(mirna_list){
  mirna_5p<-mirna_data[which(mirna_data$MIRNA%in%mirna_list), c('MIRNA', 'STRAND')]
  mirna_5p<-mirna_data[which(mirna_data$STRAND=='5P'),'MIRNA']
  return(mirna_5p)
}


all_effective_number<-data.frame(MIRNA=character())
all_rpm<-data.frame(MIRNA=character())
for(t in tumors){
  effective_file<-paste0('/Users/chens22/Documents/miRNA/', t, '/summary_files/', t, '.isomir.tsv')
  if(!file.exists(effective_file)){next}
  effective_data<-read.table(file=effective_file, sep='\t', header=T)
  effective_data<-unique(effective_data)

  normal_patients<-getNormalPatients(t)
  if(length(normal_patients)<1){next}
  effective_data<-effective_data[which(effective_data$SAMPLE%in%normal_patients),]
  effective_data<-effective_data[which(effective_data$MIRNA%in%get5PMirna(effective_data$MIRNA)),]
  effective_data$RPM<-(effective_data$READS/effective_data$TOTAL.READS.IN.SAMPLE)*1000000
  effective_data[,c('SAMPLE', 'READS', 'TOTAL.READS.IN.SAMPLE', 'SEED')]<-NULL
  
  
  effective_number<-aggregate(effective_data$EFFECTIVE.NUMBER, by=list(MIRNA=effective_data$MIRNA), mean)
  colnames(effective_number)[which(names(effective_number)=='x')]<-paste0(t)
  rpm<-aggregate(effective_data$RPM, by=list(MIRNA=effective_data$MIRNA), mean, na.rm=T)
  colnames(rpm)[which(names(rpm)=='x')]<-paste0(t, '.RPM')
  
  all_rpm<-merge(all_rpm, rpm, by='MIRNA', all=T, na.rm=T)
  all_effective_number<-merge(all_effective_number, effective_number, by='MIRNA', all=T, na.rm=T)
}

all_effective_number$RPM<-rowMeans(all_rpm[, 2:ncol(all_rpm)], na.rm=T)
effective_file<-paste0(base_path, 'normal.patients.isomir.effective.tsv')
write.table(all_effective_number, file=effective_file, sep='\t', row.names=F)
