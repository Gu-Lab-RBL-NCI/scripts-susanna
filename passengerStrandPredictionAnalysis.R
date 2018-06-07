
## PURPOSE: return prediction of undocumented mature sequence based on motifs generated from TCGA data
## INPUT: TCGA motifs   TCGA.motifs.tsv
## OUPUT: passenger strand prediction table   passenger.strand.prediction.analysis.tsv

base_path <- '/Volumes/2TB (MAC)/Susanna/'


tcga_motifs_data <- read.table(paste0(doc_path, 'TCGA.motifs.tsv'), sep='\t', header=T, stringsAsFactors = F)
tcga_motifs_data<-tcga_motifs_data[which(tcga_motifs_data$EXISTS==F),]
#tcga_motifs_data <- tcga_motifs_data[complete.cases(motifs_data), ]

passenger_data <- tcga_motifs_data[,c('PRIMIRNA', 'STRAND', 'PREDICTED.SEQUENCE', 'UNIQUE', 'MOTIF.LEN', 'MOTIF')]

for(i in 1:nrow(passenger_data)){
  
  motif <- toString(passenger_data[i, 'MOTIF'])
  primirna<- toString(passenger_data[i, 'PRIMIRNA'])
  
  motif_file <- paste0(doc_path, motif, '.', primirna, '.reads.tsv')
  if(file.exists(motif_file)){
    motif_data <- read.table(motif_file, sep='\t', header=T, stringsAsFactors = F, check.names = F)
  }else{
    
    next
  }

  max_reads <- max(as.numeric(motif_data$SUM), na.rm=T)
  total_reads <- sum(as.numeric(motif_data$SUM), na.rm=T)
  passenger_data[i, 'PATIENTS'] <- ncol(motif_data[,!colnames(motif_data)%in%c('SUM', 'SEQUENCE', 'DISTANCE.FROM.PREDICTED', 'TEMPLATED'), drop=F])
  
  
  motif_data <- motif_data[which(motif_data$TEMPLATED==T),]
  if(nrow(motif_data) < 1){ 
    passenger_data[i, 'TEMPLATED'] <- F
    next
  }
  passenger_data[i, 'TEMPLATED'] <- T
  
  motif_data <- motif_data[which(motif_data$SUM == (max(as.numeric(motif_data$SUM), na.rm=T)) & motif_data$SUM > 2),]
  
  # motif_data <- motif_data[which(motif_data$DISTANCE.FROM.PREDICTED == min(motif_data$DISTANCE.FROM.PREDICTED)),]  
  
  if(nrow(motif_data) > 1){
    motif_data <- motif_data[which(motif_data$SUM == max(motif_data$SUM)),]  
    motif_data <- motif_data[which(nchar(motif_data$SEQUENCE) == min(nchar(motif_data$SEQUENCE))),]  
  }else if(nrow(motif_data) < 1){
    tcga_motifs_data[i, 'READS'] <- 1
    next
  }
  
  reads <- motif_data[, 'SUM']
  seq <- motif_data[, 'SEQUENCE']
  dist <- motif_data[, 'DISTANCE.FROM.PREDICTED']
  
  passenger_data[i, 'READS'] <- reads
  passenger_data[i, 'MAX.READS'] <- max_reads
  passenger_data[i, 'TOTAL.READS'] <- total_reads
  passenger_data[i, 'RATIO'] <- reads/total_reads
  passenger_data[i, 'TCGA.SEQUENCE'] <- seq
  passenger_data[i, 'DISTANCE'] <- dist
  
}

correct_predict <- passenger_data[which(passenger_data$TEMPLATED==T &
                                          (passenger_data$READS/passenger_data$MAX.READS)>=0.5),]
write.table(passenger_data, file=paste0(base_path, 'passenger/passenger.strand.prediction.analysis.tsv'), sep='\t', row.names = F)
