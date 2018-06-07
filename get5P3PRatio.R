
## PURPOSE: get ratio between reads of 5P and 3P strands
## INPUT: fidelity data   tumor_fidelity_info.tsv
##        miRBase21 data    miRBase21-master.tsv
## OUTPUT: plot of 5P to 3P reads ratio

mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
base_path<-'/Users/chens22/Documents/miRNA/'
fidelity_file_ext<-'_fidelity_info.tsv'

tumors<-c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
          'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
          'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM') 

curr_mirna<-'hsa-mir-21'

get5PStrand<-function(mirna, mirna_data){
  
  temp_mirna_data<-mirna_data[which(mirna_data$PRIMIRNA==mirna),]
  strand<-toString(temp_mirna_data[which(temp_mirna_data$STRAND=='5P'), 'MIRNA'])
  
  return(strand)
  
}

get3PStrand<-function(mirna, mirna_data){
  
  temp_mirna_data<-mirna_data[which(mirna_data$PRIMIRNA==mirna),]
  strand<-toString(temp_mirna_data[which(temp_mirna_data$STRAND=='3P'), 'MIRNA'])
  
  return(strand)
  
}

getTumorReads<-function(mirna, tumor){
  
  tumor_file<-paste(base_path, tumor, '/', tumor, fidelity_file_ext, sep='')
  
  tumor_fid_data<-read.table(file=tumor_file, header=TRUE, sep='\t')
  
  mirna_reads<-tumor_fid_data[which(tumor_fid_data$MIRNA==mirna), 'RPM.AVERAGE']
  
  return(mirna_reads)
  
}

mirna_data<-read.table(file=mirna_file, header=TRUE, sep='\t')

reads_data<-data.frame()

for(t in tumors){
  
  mirna_strands<-c(get5PStrand(curr_mirna, mirna_data), get3PStrand(curr_mirna, mirna_data))
  
  for(m in mirna_strands){
    curr_reads<-getTumorReads(m, t)
    reads_data[t, m]<-curr_reads
  }
  
  total_reads<-sum(reads_data[t,])
  
  reads_data[t,]<-reads_data[t,]/total_reads
  
}

matplot(reads_data, type = c("b"),pch=19,col = 1:2, xaxt='n', xlab='Tumors', ylab='Reads Ratio', main=curr_mirna)
axis(1, at=1:nrow(reads_data), labels=rownames(reads_data))
legend("topleft", legend = colnames(reads_data), col=1:2, pch=19)
