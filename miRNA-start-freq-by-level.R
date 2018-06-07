
## PURPOSE: get start frequency of isomiRs by high and low fidelity patients
## INPUT: miRBase data  miRBase-master.tsv
##        isomiR summary data   sample.isomir.tsv
##        isomiR cleavage data  sample.isomir.cleavage.tsv
## OUTPUT: table start frequency  tumor_fidelity level_fidelity_info.tsv

base_path<-'/Users/chens22/Documents/miRNA/'
mirna_summary_file<-'/Users/chens22/Documents/miRNA/miRBase21-master.tsv'

tumors<-c('UCEC')



getLowFidelityPatients<-function(data){
  temp_data<-NULL
  temp_data<-data[which(data$FIDELITY.LEVEL=='LOW'),]
  print(temp_data)
  result<-c()
  
  for(i in 1:nrow(temp_data)){
    
    result<-c(result, strsplit(toString(temp_data[i, 'PATIENT.ID']), '.', fixed=TRUE)[[1]][1])
    
  }
  return(unique(result))
}


getHighFidelityPatients<-function(data){
  temp_data<-NULL
  temp_data<-data[which(data$FIDELITY.LEVEL=='HIGH'),]
  print(temp_data)
  result<-c()
  
  for(i in 1:nrow(temp_data)){
    
    result<-c(result, strsplit(toString(temp_data[i, 'PATIENT.ID']), '.', fixed=TRUE)[[1]][1])
    
  }
  return(unique(result))
}


#returns dataframe of MIRNA and FIDELITY.5P from summary files
get_fidelity<-function(file_name){ #returns dataframe containing vectors? of MIRNA and FIDELITY.5P
  NHIS<-read.table(file = file_name, sep = '\t', header = TRUE)
  result<-data.frame(MIRNA=NHIS$MIRNA, FIDELITY=NHIS$FIDELITY.5)
  
  curr_patient<-strsplit(toString(basename(file_name)), ".", fixed=TRUE)[[1]][1]
  colnames(result)[colnames(result)=='FIDELITY']<-paste('FIDELITY.5P-', curr_patient, sep='')
  
  return(result)
}

#returns dataframe of MIRNA and TOTAL.READS from summary files
get_rpm<-function(file_name){ #returns a dataframe containing vectors of MIRNA and RPM: (TOTAL.READS/TOTAL.READS.IN.SAMPLE)1000000
  NHIS<-read.table(file = file_name, sep = '\t', header = TRUE)
  result<-data.frame(MIRNA = NHIS$MIRNA, RPM = (NHIS$TOTAL.READS/NHIS$TOTAL.READS.IN.SAMPLE)*1000000)
  
  curr_patient<-strsplit(toString(basename(file_name)), ".", fixed=TRUE)[[1]][1]
  colnames(result)[colnames(result)=='RPM']<-paste('RPM-', curr_patient, sep='')
  return(result)
}

get_cleavage_ratio<-function(file_name){
  
  NHIS<-read.table(file = file_name, sep = '\t', header = TRUE)
  reads<-NHIS$EXPRESSION.READS
  
  ratio_3<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START..3/reads)
  ratio_2<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START..2/reads)
  ratio_1<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START..1/reads)
  ratio0<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START.0/reads)
  ratio1<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START.1/reads)
  ratio2<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START.2/reads)
  ratio3<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$START.3/reads)
  ratio_other<-data.frame(MIRNA=NHIS$MIRNA, RATIO=NHIS$OTHER/reads)
  
  curr_patient<-strsplit(toString(basename(file_name)), ".", fixed=TRUE)[[1]][1]
  colnames(ratio_3)[colnames(ratio_3) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio_2)[colnames(ratio_2) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio_1)[colnames(ratio_1) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio0)[colnames(ratio0) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio1)[colnames(ratio1) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio2)[colnames(ratio2) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio3)[colnames(ratio3) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  colnames(ratio_other)[colnames(ratio_other) == 'RATIO'] <- paste('RATIO-', curr_patient, sep='')
  
  return(list(ratio_3, ratio_2, ratio_1, ratio0, ratio1, ratio2, ratio3, ratio_other))
  
}

get_paralogs<-function(mirna_file){
  NHIS<-read.table(file=mirna_file, sep='\t', header=TRUE)
  paralogs<-NHIS$PARALOGS>0
  mirna<-NHIS$MIRNA
  return(data.frame(MIRNA=mirna, PARALOGS=paralogs))
}

get_strand<-function(mirna_file){
  NHIS<-read.table(file=mirna_file, sep='\t', header=TRUE)
  strand<-NHIS$STRAND
  mirna<-NHIS$MIRNA
  return(data.frame(MIRNA=mirna, STRAND=strand))
}


getStartFreqLevel<-function(patients, cleavage_dir, summary_dir, tumor, level){
  
  #get fidelity and RPM from each file in summary files
  fidelity<-data.frame(MIRNA=NA, x=NA, y=NA)
  rpm<-data.frame(MIRNA=NA, x=NA, y=NA)
  
  for(i in summary_dir){
    curr_patient<-strsplit(basename(i), '.', fixed=TRUE)[[1]][1]
    if(curr_patient%in%patients){
      curr_fid<-get_fidelity(i)
      curr_rpm<-get_rpm(i)
      
      fidelity<-merge(fidelity, curr_fid, by='MIRNA', na.rm=TRUE, all=TRUE)
      rpm<-merge(rpm, curr_rpm, by='MIRNA', na.rm=TRUE, all=TRUE)
    }
  }
  
  cleavage_freq_3<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq_2<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq_1<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq0<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq1<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq2<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq3<-data.frame(MIRNA=NA, x=NA, y=NA)
  cleavage_freq_other<-data.frame(MIRNA=NA, x=NA, y=NA)
  
  
  for(f in cleavage_dir){
  
    curr_patient<-strsplit(basename(f), '.', fixed=TRUE)[[1]][1]
  
    if(curr_patient%in%patients){
      list_freqdf<-get_cleavage_ratio(f)
      curr_cleavage_freq_3<-list_freqdf[[1]]
      curr_cleavage_freq_2<-list_freqdf[[2]]
      curr_cleavage_freq_1<-list_freqdf[[3]]
      curr_cleavage_freq0<-list_freqdf[[4]]
      curr_cleavage_freq1<-list_freqdf[[5]]
      curr_cleavage_freq2<-list_freqdf[[6]]
      curr_cleavage_freq3<-list_freqdf[[7]]
      curr_cleavage_freq_other<-list_freqdf[[8]]
      
      cleavage_freq_3<-merge(cleavage_freq_3, curr_cleavage_freq_3, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq_2<-merge(cleavage_freq_2, curr_cleavage_freq_2, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq_1<-merge(cleavage_freq_1, curr_cleavage_freq_1, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq0<-merge(cleavage_freq0, curr_cleavage_freq0, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq1<-merge(cleavage_freq1, curr_cleavage_freq1, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq2<-merge(cleavage_freq2, curr_cleavage_freq2, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq3<-merge(cleavage_freq3, curr_cleavage_freq3, by='MIRNA', na.rm=TRUE, all=TRUE)
      cleavage_freq_other<-merge(cleavage_freq_other, curr_cleavage_freq_other, by='MIRNA', na.rm=TRUE, all=TRUE)
    }
    
  }
  
  
  
  cleavage_freq_3$START..3.AVERAGE<-cbind(rowMeans(cleavage_freq_3[, 2:ncol(cleavage_freq_3)], na.rm=TRUE))
  cleavage_freq_3<-transform(cleavage_freq_3, START..3.SD=apply(cleavage_freq_3[,2:(ncol(cleavage_freq_3)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq_2$START..2.AVERAGE<-cbind(rowMeans(cleavage_freq_2[, 2:ncol(cleavage_freq_2)], na.rm=TRUE))
  cleavage_freq_2<-transform(cleavage_freq_2, START..2.SD=apply(cleavage_freq_2[,2:(ncol(cleavage_freq_2)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq_1$START..1.AVERAGE<-cbind(rowMeans(cleavage_freq_1[, 2:ncol(cleavage_freq_1)], na.rm=TRUE))
  cleavage_freq_1<-transform(cleavage_freq_1, START..1.SD=apply(cleavage_freq_1[,2:(ncol(cleavage_freq_1)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq0$START.0.AVERAGE<-cbind(rowMeans(cleavage_freq0[, 2:ncol(cleavage_freq0)], na.rm=TRUE))
  cleavage_freq0<-transform(cleavage_freq0, START.0.SD=apply(cleavage_freq0[,2:(ncol(cleavage_freq0)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq1$START.1.AVERAGE<-cbind(rowMeans(cleavage_freq1[, 2:ncol(cleavage_freq1)], na.rm=TRUE))
  cleavage_freq1<-transform(cleavage_freq1, START.1.SD=apply(cleavage_freq1[,2:(ncol(cleavage_freq1)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq2$START.2.AVERAGE<-cbind(rowMeans(cleavage_freq2[, 2:ncol(cleavage_freq2)], na.rm=TRUE))
  cleavage_freq2<-transform(cleavage_freq2, START.2.SD=apply(cleavage_freq2[,2:(ncol(cleavage_freq2)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq3$START.3.AVERAGE<-cbind(rowMeans(cleavage_freq3[, 2:ncol(cleavage_freq3)], na.rm=TRUE))
  cleavage_freq3<-transform(cleavage_freq3, START.3.SD=apply(cleavage_freq3[,2:(ncol(cleavage_freq3)-1)], 1, sd, na.rm=TRUE))
  
  cleavage_freq_other$START.OTHER.AVERAGE<-cbind(rowMeans(cleavage_freq_other[, 2:ncol(cleavage_freq_other)], na.rm=TRUE))
  cleavage_freq_other<-transform(cleavage_freq_other, START.OTHER.SD=apply(cleavage_freq_other[,2:(ncol(cleavage_freq_other)-1)], 1, sd, na.rm=TRUE))
  
  fidelity$FIDELITY.AVERAGE<-cbind(rowMeans(fidelity[, 2:ncol(fidelity)], na.rm=TRUE))
  fidelity<-transform(fidelity, FIDELITY.SD=apply(fidelity[,2:(ncol(fidelity)-1)], 1, sd, na.rm=TRUE)) #sort of not legit; subtracted one for new average column
  
  rpm$RPM.AVERAGE<-cbind(rowMeans(rpm[, 2:ncol(rpm)], na.rm=TRUE))
  rpm<-transform(rpm, RPM.SD=apply(rpm[,2:(ncol(rpm)-1)], 1, sd, na.rm=TRUE)) #sort of not legit; subtracted one for new average column
  
  paralogs<-get_paralogs(mirna_summary_file)
  
  strand<-get_strand(mirna_summary_file)
  
  
  
  fidelity_5p<-data.frame(MIRNA=fidelity$MIRNA)
  
  fidelity_5p<-merge(fidelity_5p, paralogs, by='MIRNA', all=TRUE, na.rm=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, strand, by='MIRNA', all=TRUE, na.rm=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, fidelity[,c('MIRNA', 'FIDELITY.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, fidelity[,c('MIRNA', 'FIDELITY.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, rpm[,c('MIRNA', 'RPM.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, rpm[,c('MIRNA', 'RPM.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_3[,c('MIRNA', 'START..3.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_3[,c('MIRNA', 'START..3.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_2[,c('MIRNA', 'START..2.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_2[,c('MIRNA', 'START..2.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_1[,c('MIRNA', 'START..1.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_1[,c('MIRNA', 'START..1.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq0[,c('MIRNA', 'START.0.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq0[,c('MIRNA', 'START.0.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq1[,c('MIRNA', 'START.1.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq1[,c('MIRNA', 'START.1.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq2[,c('MIRNA', 'START.2.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq2[,c('MIRNA', 'START.2.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq3[,c('MIRNA', 'START.3.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq3[,c('MIRNA', 'START.3.SD')], by='MIRNA', all=TRUE)
  
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_other[,c('MIRNA', 'START.OTHER.AVERAGE')], by='MIRNA', all=TRUE)
  fidelity_5p<-merge(fidelity_5p, cleavage_freq_other[,c('MIRNA', 'START.OTHER.SD')], by='MIRNA', all=TRUE)
  
  
  fidelity_5p<-fidelity_5p[complete.cases(fidelity_5p),]
  fidelity_5p<-fidelity_5p[!duplicated(fidelity_5p),]
  
  result_file_name<-paste(base_path, tumor, '/', tumor, '_', level, '_fidelity_info.tsv', sep='')
  write.table(fidelity_5p, file=result_file_name, quote=FALSE, sep='\t', col.names = NA)
  
}

for (curr_tumor in tumors){
  curr_tumor<-tumors[1]
  summary_files<-paste(base_path, curr_tumor, '/summary', sep='')
  patient_fid_file<-paste(base_path, curr_tumor, '/', curr_tumor, '_patient_fidelity_level.csv', sep='')
  
  summary_data<-list.files(path=summary_files, pattern='*.tsv', full.names=TRUE)
  patient_fid_level_data<-read.csv(file=patient_fid_file, header=TRUE)

  
  low_fid_patients<-getLowFidelityPatients(patient_fid_level_data)
  
  high_fid_patients<-getHighFidelityPatients(patient_fid_level_data)
  
  
  cleavage_files<-paste(base_path, curr_tumor, '/isomir_cleavage', sep='')
  cleavage_data<-list.files(path=cleavage_files, pattern='*.tsv', full.names=TRUE)
  
  getStartFreqLevel(low_fid_patients, cleavage_data, summary_data, curr_tumor, 'low')
  
  getStartFreqLevel(high_fid_patients, cleavage_data, summary_data, curr_tumor, 'high')
  
}
