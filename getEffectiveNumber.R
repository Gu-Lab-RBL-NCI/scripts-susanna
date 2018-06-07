
## PURPOSE: calculate effective number of miRNA for specified tumors
## INPUT: miRBase data                    miRBase21-master.tsv
##        manifest data                 all-tumor-manifest.tsv
##        isomiR expresson data     tumor.isomir.expression.tsv
##        isomiR summary data      tumor.isomir.expression.tsv
## OUTPUT: table containing effective number of miRNA   tumor.isomir.effective.tsv

require('data.table')
require('stringr')

base_path<-'/Users/chens22/Documents/miRNA/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=T)

expression_ext<-'.isomir.expression.tsv'
summary_ext<-'.isomir.tsv'

tumors<-c('BRCA')

manifest_file<-'/Users/chens22/Documents/miRNA/all-tumor-manifest.csv'
manifest_data<-read.csv(file=manifest_file, header=T)

getNormalPatients<-function(tumor){
  patients<-manifest_data[which(manifest_data$SAMPLE.TYPE=='SOLID.TISSUE.NORMAL' & manifest_data$DISEASE.ABBV==tumor), 'NAME']
  return(patients)
}

getEffectiveNumber<-function(reads){
  reads_freq<-reads/sum(reads)
  effective<-1/sum(reads_freq^2)
  return(effective)
}

getGosolovNumber<-function(reads){
  reads_freq<-reads/sum(reads)
  reads_freq<-reads_freq[order(-reads_freq)]
  gosolov<-sum(reads_freq/(reads_freq+reads_freq[1]-(reads_freq)^2))
  return(gosolov)
}
                            

for(t in tumors){
  expression_file<-paste0(base_path, t, '/summary_files/', t, expression_ext)
  expression_data<-read.table(file=expression_file, header = T, fill = F, sep = '\t', stringsAsFactors = F)
  expression_data$SEED<-substr(expression_data$SEQUENCE, 2,8)
  unique_samples<-unique(expression_data$SAMPLE)
  for(s in unique_samples){
    new_sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
    expression_data[which(expression_data$SAMPLE==s), 'SAMPLE']<-new_sample
  }
  write.table(expression_data, file=expression_file, sep='\t', row.names=F)
  summary_file<-paste0(base_path, t, '/summary_files/', t, summary_ext)
  summary_data<-read.table(file=summary_file, header = T, fill = F, sep = '\t', stringsAsFactors = F)
  unique_samples<-unique(summary_data$SAMPLE)
  for(s in unique_samples){
    new_sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
    summary_data[which(summary_data$SAMPLE==s), 'SAMPLE']<-new_sample
  }
  write.table(summary_data, file=summary_file, sep='\t', row.names=F)
  # 
  unique_samples<-unique(summary_data$SAMPLE)
  unique_samples<-unique_samples[unique_samples%in%getNormalPatients(t)]
  final_seeds_data<-data.frame(MIRNA=character(), SAMPLE=character(), stringsAsFactors=F)
  
  #unique_samples<-getNormalPatients(t)
  count<-0
  for(s in unique_samples){
    count<-count+1
    print(count/length(unique_samples))
    #expression_file<-paste0(base_path, t, '/sequence_info/', s, '.converted.unpaired.fastq', expression_ext)
    #if(!file.exists(expression_file)){
    #  print(t)
     # next
    #}
    curr_exp_data<-subset(expression_data, SAMPLE==s)
    #curr_exp_data<-read.table(file=expression_file, sep='\t', header=T, stringsAsFactors = F)
    #curr_exp_data$SAMPLE<-s
    #curr_exp_data$SEED<-substr(curr_exp_data$SEQUENCE, 2,8)
    seeds_data<-aggregate(curr_exp_data$READS, by=list(MIRNA=curr_exp_data$MIRNA, SEED=curr_exp_data$SEED), sum)
    colnames(seeds_data)[which(names(seeds_data)=='x')]<-'READS'
    seeds_data<-seeds_data[order(-seeds_data$READS),] 
    seeds_data$SAMPLE<-s
    
    unique_mirna<-unique(seeds_data$MIRNA)
    
    for(m in unique_mirna){
      curr_seeds_data<-seeds_data[which(seeds_data$MIRNA==m),]
      curr_seeds_data<-curr_seeds_data[order(-curr_seeds_data$READS),] 
      
      curr_reads<-curr_seeds_data$READS
      curr_reads<-curr_reads[order(-curr_seeds_data$READS)]
      
      effective<-getEffectiveNumber(curr_reads)
      gosolov<-getGosolovNumber(curr_reads)
      seeds_data[which(seeds_data$MIRNA==m), 'EFFECTIVE.NUMBER']<-effective
      seeds_data[which(seeds_data$MIRNA==m), 'GOSOLOV.NUMBER']<-gosolov
    }
    final_seeds_data<-rbind(final_seeds_data, seeds_data)
  }
  unique_samples<-unique(final_seeds_data$SAMPLE)
  for(s in unique_samples){
    new_sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
    final_seeds_data[which(final_seeds_data$SAMPLE==s), 'SAMPLE']<-new_sample
  }
  unique_samples<-unique(summary_data$SAMPLE)
  for(s in unique_samples){
    new_sample<-strsplit(toString(s), '.', fixed=T)[[1]][1]
    summary_data[which(summary_data$SAMPLE==s), 'SAMPLE']<-new_sample
  }
  summary_data<-unique(summary_data[,c('SAMPLE', 'TOTAL.READS.IN.SAMPLE')])
  summary_data<-summary_data[which(summary_data$SAMPLE%in%final_seeds_data$SAMPLE),]
  final_seeds_data<-merge(final_seeds_data, summary_data, by='SAMPLE', all=T, na.rm=T)

  final_seeds_data<-unique(final_seeds_data)
  seeds_file<-paste0(base_path, t, '/summary_files/', t, '.isomir.effective.tsv')
  write.table(final_seeds_data, file=seeds_file, sep='\t', row.names=F)
}



# 
# colnames(output) <- c("MIRNA", "EiSimpson", "EGolosov")
# reads <-aggregate(file$READS, by=list(MIRNA=file$MIRNA), sum)
# reads <- merge(reads, output, by="MIRNA")
# reads <- reads[order(-reads$x),] 
# reads$strand5p <- regexpr('5p', reads$MIRNA)
# mir5p <- reads[which(reads$strand5p>0),]
# mir5p <- mir5p[c(1:100),]
# mir5p <- mir5p[order(mir5p$EiSimpson),] 
# write.table(mir5p, "mir5p.txt", sep="\t", append = FALSE, row.names = F)
# 
# seeds[1,3]+seeds[2,3]+seeds[3,3]+seeds[4,3]+seeds[5,3]+seeds[6,3]+seeds[7,3]
