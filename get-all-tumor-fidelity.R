
## PURPOSE: compile miRNA fidelity 5P across all tumors into one table
## INPUT: fidelity 5P data of all tumors      tumor_fidelity_info.tsv
## OUTPUT: compiled fidelity data     all_tumors_fidelity_summary.tsv

tumors<-c('GBM','CHOL','DLBC','UCS','ACC','UVM','MESO','KICH',
          'LAML','THYM','TGCT','READ','PAAD','PCPG','ESCA',
          'SARC','CESC','KIRP','LIHC','LGG','BLCA','SKCM',
          'COAD','STAD','OV','LUSC','PRAD','LUAD','HNSC',
          'THCA','UCEC','KIRC', 'BRCA')
base_path<-'/Users/chens22/Documents/miRNA/'

fid_summary_data<-data.frame(MIRNA=character(), RPM.AVERAGE=integer())

for(curr_tumor in tumors){
  curr_file<-paste(base_path, curr_tumor, '/', curr_tumor, '_fidelity_info.tsv', sep='')
  NHIS<-read.table(file=curr_file, sep='\t', header=TRUE)

  NHIS<-NHIS[which(NHIS$RPM.AVERAGE>100),c('MIRNA', 'FIDELITY.AVERAGE', 'RPM.AVERAGE')]
  colnames(NHIS)[colnames(NHIS)=='RPM.AVERAGE']<-'TEMP.RPM'
  
  fid_summary_data<-merge(fid_summary_data, NHIS, by='MIRNA', all=TRUE, na.rm=TRUE)
  fid_summary_data[,'RPM.AVERAGE']<-rowSums(fid_summary_data[,c('RPM.AVERAGE', 'TEMP.RPM')], na.rm=TRUE)
  fid_summary_data[,'TEMP.RPM']<-NULL
  
  colnames(fid_summary_data)[colnames(fid_summary_data)=='FIDELITY.AVERAGE']<-paste(curr_tumor, 'FIDELITY.AVERAGE', sep='.')
}
fid_summary_data[,'RPM.AVERAGE']<-fid_summary_data[,'RPM.AVERAGE']/length(tumors)
fid_summary_data<-fid_summary_data[-nrow(fid_summary_data),]
fid_summary_data<-fid_summary_data[which(fid_summary_data$RPM.AVERAGE>100),]

result_filename<-paste(base_path, 'all_tumors_fidelity_summary.tsv', sep='')
write.table(fid_summary_data, file=result_filename, sep='\t', row.names = FALSE)
