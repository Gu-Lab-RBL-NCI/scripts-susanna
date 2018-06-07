## PURPOSE: calculate the fidelity 5p and cleavage frequency for all miRNA
## INPUT: miRBase21 data                 miRBase21-master.tsv
##        manifest data                all-tumor-manifest.csv
##        isomiR summary data                tumor.isomiR.tsv
##        isomiR expression data  tumor.isomir.expression.tsv
## OUTPUT: table containing average RPM and FIDELITY 5P and cleavage frequency for all miRNA in a tumor


install.packages('ggplot2')
install.packages("viridis") # dependency
install.packages("devtools")
devtools::install_github("ropensci/plotly")
install.packages('dplyr')


base_path<-'/Users/chens22/Documents/'

mirna_file<- paste0(base_path, '/miRBase21/miRBase21-master.tsv')
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors = F)
manifest_file<- paste0(base_path, 'all-tumor-manifest.csv')
manifest_data<-read.csv(file=manifest_file, header=T, stringsAsFactors=F)
summary_ext<-'.isomir.tsv'
expression_ext<-'.isomir.expression.tsv'
cleavage_ext<-'.isomir.cleavage.tsv'
fidelity_ext<-'.fidelity.tsv'


#return dataframe of MIRNA and TOTAl.READS from miRBase21-results-final
get_reads<-function(summ_data){
  return(data.frame(MIRNA = summ_data$MIRNA, TOTAL.READS=summ_data$TOTAL.READS))
}

#return dataframe of MIRNA and CONSENSUS from miRBase21-results-final
get_consensus<-function(summ_data){
  return(data.frame(MIRNA = summ_data$MIRNA, CONSENSUS=summ_data$CONSENSUS))
}

#return dataframe of MIRNA and MOTIF from miRBase21-final-results
get_motif<-function(summ_data){
  return(data.frame(MIRNA = summ_data$MIRNA, MOTIF=summ_data$MOTIF))
}

#return patient id from filename
get_patient_id<-function(sample){
  id<-strsplit(toString(sample), ".", fixed=TRUE)[[1]][1]
  return(id)
}

get_patient<-function(summ_data){
  return(unique(summ_data$SAMPLE))
}

getTemplated5P<-function(mirna){
  
  extended_seqs<-unique(mirna_data[which(mirna_data$MIRNA==mirna), c('EXTENDED.SEQUENCE')])
  mature_seq<-toString(unique(mirna_data[which(mirna_data$MIRNA==mirna), 'SEQUENCE']))
  
  templated<-c()
  for(e in extended_seqs){
    start<-gregexpr(mature_seq, e, fixed=T)[[1]]
    curr_templated<-substr(e, start-5, start-1)
    templated<-c(templated, curr_templated)
  }
  return(unique(templated))
  
}

getTemplated5PData<-function(){
  unique_mirna<-unique(mirna_data$MIRNA)
  data<-data.frame(MIRNA=unique_mirna, TEMPLATED.5P=NA)
  for(m in unique_mirna){
    data[which(data$MIRNA==m), 'TEMPLATED.5P']<-paste(getTemplated5P(m), collapse=',')
  }
  return(data)
}

#
get_start_pos<-function(templated, consensus, sequence){
  seq_start<-1
  temp_seq<-substr(sequence, seq_start, seq_start+7)
  found = regexpr(pattern=temp_seq, consensus)
  
  if(found!=-1){ #return found if sequence starts within consensus
    return(found)
  }
  
  while(found==-1 & seq_start < length(sequence)+8){
    seq_start<-seq_start+1
    temp_seq<-substr(sequence, seq_start, seq_start+7)
    
    found<-regexpr(pattern=temp_seq, consensus)
  }
  
  curr_templated<-substr(sequence, 1, seq_start-1) #finds the nucleotides of the sequence up until the consensus is found
  
  for(i in 1:length(templated)){
    templated[i]<-substr(templated[i], nchar(templated[i])-(nchar(curr_templated)-1), nchar(templated[i])) #edits the templateds so that their length matches curr_templated
  }
  templated<-curr_templated%in%templated
  
  if(found==1 & templated){
    return((found-seq_start)+1)
  }else{
    return(NA)
  }
}


get_start_freq<-function(expression_data, consensus){ #templated: dataframe(MIRNA, templated nucleotides), expression_data: isomir expression data of patient, consensus: dataframe(MIRNA, consensus)
  MIRNA<-unique(expression_data$MIRNA)  #get unique MIRNA expressed in patient
  freq<-data.frame(MIRNA=unique(expression_data$MIRNA)) #intialize dataframe(MIRNA, start positions) to hold counts of each start position at each miRNA
  freq[,'START.OTHER']<-0
  freq[,'START.-3']<-0
  freq[,'START.-2']<-0
  freq[,'START.-1']<-0
  freq[,'START.0']<-0
  freq[,'START.1']<-0
  freq[,'START.2']<-0
  freq[,'START.3']<-0
  
  for(m in MIRNA){ #iterates through every miRNA expressed in patient
    temp_data<-subset(expression_data, MIRNA==m) #creates subset of expression file where MIRNA = curr MIRNA
    curr_con<-consensus[which(consensus$MIRNA==m), 'CONSENSUS'] #gets consensus for curr MIRNA
    curr_templated<-getTemplated5P(m)
    total_reads<-0
    for(i in 1:nrow(temp_data)){ #goes through every sequence in subset
      start_position<-get_start_pos(curr_templated, curr_con, temp_data[i, 'SEQUENCE'])
      curr_reads<-temp_data[i, 'READS']
      total_reads<-total_reads+curr_reads
      
      if(start_position-1 < -3 || start_position-1 > 3 || is.na(start_position)){
        freq[MIRNA==m, 'START.OTHER']<-freq[MIRNA==m, 'START.OTHER']+curr_reads
      }else{
        freq[MIRNA==m, paste('START.', toString(start_position-1), sep='')]<-freq[MIRNA==m, paste('START.', toString(start_position-1), sep='')]+curr_reads
      }
      
    }
    freq[MIRNA==m, 'EXPRESSION.READS'] <- total_reads
  }
  names(freq)[names(freq) == 'START.-1'] <- 'START..1'
  names(freq)[names(freq) == 'START.-2'] <- 'START..2'
  names(freq)[names(freq) == 'START.-3'] <- 'START..3'
  return(freq)
}


#returns dataframe of MIRNA and FIDELITY.5P from summary files
get_fidelity<-function(summary_data){ #returns dataframe containing vectors? of MIRNA and FIDELITY.5P
  return(as.vector(summary_data$FIDELITY.5P))
}

#returns dataframe of MIRNA and TOTAL.READS from summary files
get_rpm<-function(summary_data){ #returns a dataframe containing vectors of MIRNA and RPM: (TOTAL.READS/TOTAL.READS.IN.SAMPLE)1000000
  return(as.vector((summary_data$TOTAL.READS/summary_data$TOTAL.READS.IN.SAMPLE)*1000000))
}

#returns list of dataframes (MIRNA, CLEAVAGE RATIO) from isomir-cleavage files (contains reads at each position)
get_cleavage_ratio<-function(cleavage_data){
  
  reads<-cleavage_data$EXPRESSION.READS
  
  ratio_3<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START..3/reads)
  ratio_2<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START..2/reads)
  ratio_1<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START..1/reads)
  ratio0<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START.0/reads)
  ratio1<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START.1/reads)
  ratio2<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START.2/reads)
  ratio3<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START.3/reads)
  ratio_other<-data.frame(MIRNA=cleavage_data$MIRNA, RATIO=cleavage_data$START.OTHER/reads)
  
  curr_patient<-strsplit(toString(cleavage_data$SAMPLE[1]), ".", fixed=TRUE)[[1]][1]
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

keep5PMirna<-function(mirna_list){
  list<-c()
  for(m in mirna_list){
    strand<-unique(mirna_data[which(mirna_data$MIRNA==m), 'STRAND'])
    if(strand%in%'5P'){
      list<-c(list, m)
    }
  }
  return(list)
}

tumor<-'ACC'

manifest_file<-paste0(base_path, '/TARGET/TARGET-manifest.csv')
manifest_data<-read.csv(file=manifest_file, header=T, stringsAsFactors = F)
#TCGA data
summary_file<-paste0(base_path, tumor, '/summary_files/', tumor, '.isomir.', '.tsv')#paste0(base_path, tumor, '/summary_files/', tumor, summary_ext)
expression_file<-paste(base_path, tumor, '/summary_files/', tumor, '.isomir.expression.', i, '.tsv', sep='')#paste(base_path, tumor, '/summary_files/', tumor, expression_ext, sep='')

summary_data<-read.table(file = summary_file, sep = '\t', header = TRUE, stringsAsFactors=F)
expression_data<-read.table(file=expression_file, sep='\t', header=TRUE, stringsAsFactors = F)

# remove extension from names of samples
unique_samples<-unique(expression_data$SAMPLE)
for(s in unique_samples){
  sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
  expression_data[which(expression_data$SAMPLE==s), 'SAMPLE']<-sample
}
write.table(expression_data, file=expression_file, sep='\t', row.names=F)

unique_samples<-unique(summary_data$SAMPLE)
for(s in unique_samples){
  sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
  summary_data[which(summary_data$SAMPLE==s), 'SAMPLE']<-sample
}
write.table(summary_data, file=summary_file, sep='\t', row.names=F)


#get patients
patients<-unique(summary_data$SAMPLE)

curr_patients<-patients[patients%in%manifest_data[which(manifest_data$DISEASE.TYPE==t), 'NAME']]
curr_patients<-curr_patients[curr_patients%in%expression_data$SAMPLE]
cleavage_file<-paste0(base_path, tumor, '/', tumor, '.', tumor, cleavage_ext)
if(file.exists(cleavage_file)){
  mirna_cleavage_freq<-read.table(file=cleavage_file, sep='\t', header=T, stringsAsFactors = F)
}else{
  mirna_cleavage_freq<-data.frame(c())
}
temp_summ_data<-summary_data[which(summary_data$SAMPLE%in%curr_patients),]
temp_exp_data<-expression_data[which(expression_data$SAMPLE%in%curr_patients),]

# get cleavage frequency for each miRNA
for(p in curr_patients){
  curr_summ_data<-subset(temp_summ_data, SAMPLE==p)
  curr_exp_data<-subset(temp_exp_data, SAMPLE==p)
  patient_id<-get_patient_id(p)

  curr_consensus<-get_consensus(curr_summ_data)
  curr_motif<-get_motif(curr_summ_data)
  curr_start_freq<-get_start_freq(curr_exp_data, curr_consensus)
  curr_reads<-get_reads(curr_summ_data)
  
  patient_data<-data.frame(SAMPLE=c(p))
  patient_data<-merge(patient_data, Reduce(function(x,y) {merge(x,y, by='MIRNA', all=FALSE, na.rm=TRUE)}, list(getTemplated5PData(), curr_consensus, curr_motif, curr_start_freq, curr_reads)))
       
  mirna_cleavage_freq<-rbind.data.frame(mirna_cleavage_freq, patient_data)
  write.table(mirna_cleavage_freq, file=cleavage_file, quote=FALSE, sep='\t', row.names=F)
}



cleavage_data<-read.table(cleavage_file, sep='\t', header=T, stringsAsFactors = F)
cleavage_data<-merge(cleavage_data, summary_data, by=c('SAMPLE', 'MIRNA'), all=F)
write.table(cleavage_data, file=cleavage_file, sep='\t', row.names=F)


unique_mirna<-as.vector(unique(cleavage_data$MIRNA))

# add average fidelity 5p and RPM to cleavage frequency data
fidelity<-data.frame(MIRNA=unique_mirna, FIDELITY.AVERAGE = NA, FIDELITY.SD = NA, stringsAsFactors = FALSE)
rpm<-data.frame(MIRNA=unique_mirna, RPM.AVERAGE = NA, stringsAsFactors = FALSE)

cleavage_freq_3<-data.frame(MIRNA=unique_mirna)
cleavage_freq_2<-data.frame(MIRNA=unique_mirna)
cleavage_freq_1<-data.frame(MIRNA=unique_mirna)
cleavage_freq0<-data.frame(MIRNA=unique_mirna)
cleavage_freq1<-data.frame(MIRNA=unique_mirna)
cleavage_freq2<-data.frame(MIRNA=unique_mirna)
cleavage_freq3<-data.frame(MIRNA=unique_mirna)
cleavage_freq_other<-data.frame(MIRNA=unique_mirna)

for(m in unique_mirna){
  temp_summ_data<-subset(cleavage_data, cleavage_data$MIRNA==m)
  curr_fid<-get_fidelity(temp_summ_data)
  curr_rpm<-get_rpm(temp_summ_data)

  fidelity[which(fidelity$MIRNA==m),'FIDELITY.AVERAGE'] <- mean(curr_fid)
  fidelity[which(fidelity$MIRNA==m),'FIDELITY.SD'] <- sd(curr_fid)
  rpm[which(rpm$MIRNA==m),'RPM.AVERAGE'] <- mean(curr_rpm)

}
  
unique_patient <- unique(cleavage_data$SAMPLE)

for(p in unique_patient){
  temp_cleavage_data<-subset(cleavage_data, cleavage_data$SAMPLE==p)
  list_freqdf<-get_cleavage_ratio(temp_cleavage_data)
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



summary_cleavage<-data.frame(MIRNA=unique_mirna)

summary_cleavage<-merge(summary_cleavage, fidelity[,c('MIRNA', 'FIDELITY.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, fidelity[,c('MIRNA', 'FIDELITY.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, rpm[,c('MIRNA', 'RPM.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq_3[,c('MIRNA', 'START..3.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq_3[,c('MIRNA', 'START..3.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq_2[,c('MIRNA', 'START..2.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq_2[,c('MIRNA', 'START..2.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq_1[,c('MIRNA', 'START..1.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq_1[,c('MIRNA', 'START..1.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq0[,c('MIRNA', 'START.0.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq0[,c('MIRNA', 'START.0.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq1[,c('MIRNA', 'START.1.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq1[,c('MIRNA', 'START.1.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq2[,c('MIRNA', 'START.2.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq2[,c('MIRNA', 'START.2.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq3[,c('MIRNA', 'START.3.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq3[,c('MIRNA', 'START.3.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)

summary_cleavage<-merge(summary_cleavage, cleavage_freq_other[,c('MIRNA', 'START.OTHER.AVERAGE')], by='MIRNA', all=TRUE, na.rm=TRUE)
summary_cleavage<-merge(summary_cleavage, cleavage_freq_other[,c('MIRNA', 'START.OTHER.SD')], by='MIRNA', all=TRUE, na.rm=TRUE)
filename<-paste0(base_path, tumor, '/', tumor, fidelity_ext)
write.table(summary_cleavage, file=filename, quote=FALSE, sep='\t', row.names = FALSE)

