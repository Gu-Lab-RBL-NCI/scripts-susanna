
## PURPOSE: get fidelity of miRNA across all tumors
## INPUT: miRBase data                     miRBase21-master.tsv
##        fidelity data for all tumors  tumor_fidelity_info.tsv
## OUTPUT: fidelity of miRNA across all tumors  miRNA-fidelity.tsv

require('plyr')

mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=T)

primirna<-'hsa-mir-21'

tumors<-c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
          'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
          'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM') 

get5PMirna<-function(primirna){
  mirna_5p<-mirna_data[which(mirna_data$PRIMIRNA==primirna | mirna_data$MIRNA==primirna),]
  mirna_5p<-mirna_5p[which(mirna_5p$STRAND=='5P'), 'MIRNA']
  return(toString(unique(mirna_5p)))
}

getFidelityFilename<-function(tumor){
  filename<-paste('/Users/chens22/Documents/miRNA/', tumor, '/', tumor, '_fidelity_info.tsv', sep='')
  return(filename)
}

getTumorFidelityByMirna<-function(mirna, tumor){
  fid_data<-read.table(file=getFidelityFilename(tumor), sep='\t', header=T)
  mirna_fidelity<-fid_data[which(fid_data$MIRNA==get5PMirna(mirna)),]
  return(mirna_fidelity)
}

result<-data.frame()
for(t in tumors){
  curr_fid_data<-getTumorFidelityByMirna(primirna, t)
  curr_fid_data$MIRNA<-t
  result<-rbind.fill(result, curr_fid_data)
}

result$X<-NULL
result_filename<-paste('/Users/chens22/Documents/miRNA/', primirna, '/', primirna, '-fidelity.tsv', sep='')
write.table(result, file=result_filename, sep='\t', row.names = F)

