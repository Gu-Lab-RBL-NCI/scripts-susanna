
## PURPOSE: get average normal tissue fidelity of miRNA across certain tissues
## INPUT: normal tissue fidelity  tumor.SOLID.TISSUE.NORMAL.patient.fidelity.by.tumor.tsv
##        manifest data                                             all-tumor-manifest.csv
## OUTPUT: normal tissue fidelity                               normal-tissue-fidelity.tsv

base_path<-'/Users/chens22/Documents/miRNA/'
disease_abbv_file<-'/Users/chens22/Documents/all-tumor-manifest.csv'
disease_abbv_data<-read.csv(file=disease_abbv_file, header=T)

getTumorTissue<-function(abbv){
  return(disease_abbv_data[which(disease_abbv_data$ABBREVIATION==abbv), 'PRIMARY.SITE'])
}

result_file<-'/Users/chens22/Documents/miRNA/normal-tissue-fidelity.tsv'
data<-read.table(file=result_file, sep='\t', header=T, stringsAsFactors=F)
data<-data.frame(MIRNA=data$MIRNA)
target_mirna<-data$MIRNA

tumors<-c('ACC', 'BLCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
          'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
          'SARC', 'SKCM', 'STAD','TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM') 

for(t in tumors){
  fid_file<-paste(base_path, t, '/', t, '.SOLID.TISSUE.NORMAL.patient.fidelity.by.tumor.tsv', sep='')
  if(file.exists(fid_file)){
    fid_data<-read.table(file=fid_file, sep='\t', header=TRUE)
    fid_data<-fid_data[which(fid_data$MIRNA%in%target_mirna),]
    data<-merge(data, fid_data[,c('MIRNA', 'FIDELITY.AVERAGE')], by='MIRNA', all=T, na.rm=T)
    tumor_tissue<-getTumorTissue(t)
    number_of_patients<-ncol(fid_data)-3
    if(tumor_tissue=='KIDNEY' || tumor_tissue=='LUNG'){
      data$FIDELITY.AVERAGE<-data$FIDELITY.AVERAGE*number_of_patients
    }
    colnames(data)[colnames(data)=='FIDELITY.AVERAGE']<-paste(tumor_tissue, '(n=', number_of_patients, ')', sep='')
  }
}

number_of_patients<-0
for(i in grep('KIDNEY', colnames(data), fixed=TRUE)){
  curr_colname<-toString(colnames(data)[i])
  num_pos<-gregexpr('n', curr_colname, fixed=TRUE)[[1]]+2
  number_of_patients<-number_of_patients+as.numeric(substr(curr_colname, num_pos, num_pos+1))
}
data[,'KIDNEY']<-rowSums(data[,grep('KIDNEY', colnames(data))])/number_of_patients
data<-data[,-grep('KIDNEY(', colnames(data), fixed=TRUE)]
colnames(data)[colnames(data)=='KIDNEY']<-paste('KIDNEY(n=', number_of_patients, ')', sep='')

number_of_patients<-0
for(i in grep('LUNG', colnames(data), fixed=TRUE)){
  curr_colname<-toString(colnames(data)[i])
  num_pos<-gregexpr('n', curr_colname, fixed=TRUE)[[1]]+2
  number_of_patients<-number_of_patients+as.numeric(substr(curr_colname, num_pos, num_pos+1))
}
data[,'LUNG']<-rowSums(data[,grep('LUNG', colnames(data))])/number_of_patients
data<-data[,-grep('LUNG(', colnames(data), fixed=T)]
colnames(data)[colnames(data)=='LUNG']<-paste('LUNG(n=', number_of_patients, ')', sep='')


data[,'FIDELITY.AVERAGE']<-rowMeans(data[,2:ncol(data)], na.rm=TRUE)
data<-data[order(-data$FIDELITY.AVERAGE),]
data$FIDELITY.AVERAGE<-NULL

write.table(data, file=result_file, sep='\t', row.names=F)
