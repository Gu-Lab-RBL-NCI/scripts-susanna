
## PURPOSE: get alternative and canonical reads from CGC data
## INPUT: miRBase data                                                miRBase21-master.ts
##        manifest data                                               all-tumor-manifest.csv
##        isomiR summary data                                         tumor.isomir.tsv
##        isomiR expression data                                      tumor.isomir.expression.tsv
## OUTPUT: table comparing alt and canonical reads of certain miRNA   tumor.canonical.alt.read.tsv

base_path<-'/Users/chens22/Documents/miRNA/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors=F)
manifest_file<-'/Users/chens22/Documents/miRNA/all-tumor-manifest.csv'
manifest_data<-read.csv(file=manifest_file, header=T)


tumor<-'LGG'
mirna_9<-c('hsa-miR-9-5p-1-2-3', 'hsa-miR-9-3p-1-2-3')


summary_file<-paste(base_path, tumor, '/summary_files/', tumor, '.isomir.tsv', sep='')
summary_data<-read.table(file=summary_file, sep='\t', header=T)
summary_data[,'RPM']<-(summary_data$TOTAL.READS/summary_data$TOTAL.READS.IN.SAMPLE)*1000000

getFirstQuartile<-function(data){
  data<-sort(data, decreasing=F)
  median<-median(data, na.rm=T)
  data<-data[data<=median]
  return(median(data))
}

getThirdQuartile<-function(data){
  data<-sort(data, decreasing=F)
  median<-median(data, na.rm=T)
  data<-data[data>=median]
  return(median(data))
}

getSequence<-function(mirna){
  sequence<-unique(mirna_data[which(mirna_data$MIRNA==toString(mirna) | mirna_data$PRIMIRNA==toString(mirna)), 'SEQUENCE'])
  return(sequence)
}

getReadsData<-function(data){
  median<-median(data)
  maximum<-max(data)
  minimum<-min(data)
  first_quartile<-getFirstQuartile(data)
  third_quartile<-getThirdQuartile(data)
  avg_rpm<-mean(data)
  return(data.frame(MEDIAN=median, 
                   FIRST.QUARTILE=first_quartile, 
                   THIRD.QUARTILE=third_quartile, 
                   MAXIMUM=maximum, 
                   MINIMUM=minimum, 
                   RPM=avg_rpm))
}

getNormalPatients<-function(tumor){
  patients<-manifest_data[which(manifest_data$DISEASE.ABBV==tumor & manifest_data$SAMPLE.TYPE=='SOLID.TISSUE.NORMAL'), 'NAME']
  return(patients)
}

unique_mirna<-unique(summary_data$MIRNA)
data<-data.frame(MIRNA=NA, MEDIAN=NA, FIRST.QUARTILE=NA, THIRD.QUARTILE=NA, MAXIMUM=NA, MINIMUM=NA, RPM=NA)


for(m in unique_mirna){
  #if(!m%in%mirna_9){
    curr_rpm<-summary_data[which(summary_data$MIRNA==m),'RPM']
    curr_data<-getReadsData(curr_rpm)
    curr_data[,'MIRNA']<-m
    data<-rbind(data, curr_data)
  #}
}

data<-data[-1,]


filename<-paste0(base_path, tumor, '/', tumor, '.reads.tsv')
write.table(data, file=filename, sep='\t', row.names=F)

normal_patients<-getNormalPatients(tumor)
expression_file<-paste0(base_path, tumor, '/summary_files/', tumor, '.isomir.expression.tsv')
expression_data<-read.table(file=expression_file, sep='\t', header=T, stringsAsFactors = F)
expression_data<-expression_data[which(expression_data$MIRNA%in%mirna_9),]
unique_samples<-unique(expression_data$SAMPLE)
for(s in unique_samples){
  curr_sample<-strsplit(toString(s), '[.]')[[1]][1]
  expression_data[which(expression_data$SAMPLE==s), 'NAME']<-curr_sample
}
unique_samples<-unique(summary_data$SAMPLE)
for(s in unique_samples){
  curr_sample<-strsplit(toString(s), '[.]')[[1]][1]
  summary_data[which(summary_data$SAMPLE==s), 'NAME']<-curr_sample
}
summary_data<-summary_data[which(summary_data$MIRNA%in%mirna_9), c('NAME', 'MIRNA', 'TOTAL.READS')]
#expression_data<-merge(expression_data, summary_data, by=c('MIRNA', 'SAMPLE'), all=T)

data<-data.frame(MIRNA=NA, CANONICAL=NA, MEDIAN=NA, FIRST.QUARTILE=NA, THIRD.QUARTILE=NA, MAXIMUM=NA, MINIMUM=NA, RPM=NA)

for(m in mirna_9){
  canonical_seq<-getSequence(m)
  canonical_reads<-expression_data[which(expression_data$MIRNA==m & expression_data$VAR_5P==0), ]
  canonical_rpm<-c()
  samples<-unique(canonical_reads$NAME)
  samples<-samples[!samples %in% normal_patients]
  for(s in samples){
    total_reads<-unique(summary_data[which(summary_data$NAME==s & summary_data$MIRNA==m), 'TOTAL.READS'])
    reads<-sum(canonical_reads[which(canonical_reads$NAME==s), 'READS'])
    rpm<-(reads/total_reads)*1000000
    canonical_rpm<-c(canonical_rpm, rpm)
  }
  canonical_data<-getReadsData(canonical_rpm)
  canonical_data[,'MIRNA']<-m
  canonical_data[,'CANONICAL']<-T
  data<-rbind(data, canonical_data)
  
  
  alt_reads<-expression_data[which(expression_data$MIRNA==m & expression_data$VAR_5P==1),]
  alt_rpm<-c()
  samples<-unique(alt_reads$NAME)
  samples<-samples[!samples%in%normal_patients]
  for(s in samples){
    total_reads<-unique(summary_data[which(summary_data$NAME==s & summary_data$MIRNA==m), 'TOTAL.READS'])
    reads<-sum(alt_reads[which(alt_reads$NAME==s), 'READS'])
    rpm<-(reads/total_reads)*1000000
    alt_rpm<-c(alt_rpm, rpm)
  }
  alt_data<-getReadsData(alt_rpm)
  alt_data[,'MIRNA']<-m
  alt_data[,'CANONICAL']<-F
  data<-rbind(data, alt_data)
}
data<-data[-1,]
filename<-paste0(base_path, tumor, '/', tumor, '.canonical.alt.reads.tsv')
write.table(data, file=filename, sep='\t', row.names = F)



