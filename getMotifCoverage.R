## PURPOSE: get coverage of motif with different lengths
## INPUT: miRBase data  miRBase21-master.tsv
## OUTPUT: table containing motif coverage  

mirna_file <- '/Users/chens22/Documents/miRBase21/miRBase21-master.tsv'
mirna_data <- read.table(file=mirna_file, sep='\t', header=T, stringsAsFactors=F)


getMotif<-function(length, seq){ # length: motif length, 
  # seq: mature sequence to be identified with motif
  # sequences: all mature sequences excluding paralogs and duplimirna (miRNA indistinguishable from each otehr)
  mid<-trunc(nchar(seq)/2)
  start<-mid-(trunc(length/2))
  end<-mid+(trunc(length/2))
  motif<-substr(seq, start, end)
  
}

isUniqueMotif<-function(motif, sequences){
  for(i in 1:length(sequences)){
    s<-sequences[i]
    found<-gregexpr(motif, s, fixed=TRUE)[[1]]
    if(found!=-1){
      print(motif)
      print(s)
      return(F)
    }
  }
  return(T)
}



for(a in unique(mirna_data$ACCESSION)){
  
  seq <- unique(mirna_data[which(mirna_data$ACCESSION==a), 'SEQUENCE'])
  duplimirna <- c()
  
  sequences <- unique(unique_motif_data[which(unique_motif_data$SEQUENCE!=seq), 'SEQUENCE'])
  
  for(s in sequences){
    if(gregexpr(seq, s, fixed=T)[[1]][1] != -1){
      
      m <- unique_motif_data[which(unique_motif_data$SEQUENCE==s), 'MIRNA']
      duplimirna <- c(duplimirna, m)
      
    }
  }
  if(length(duplimirna) > 0){
    
    mirna_data[which(mirna_data$ACCESSION==a), 'DUPLI.ID'] <- paste(duplimirna, collapse=',')
    print(paste(duplimirna, collapse=','))
  }
  
}

# unique_motif_data <- data.frame('MIRNA'=character(),
#                                 'ACCESSION'=character(),
#                                 'UNIQUE.MOTIF'=character(),
#                                 'SEQUENCE'=character())

unique_motif_data <- unique(mirna_data[,c('MIRNA', 'ACCESSION', 'UNIQUE.MOTIF', 'SEQUENCE')])

#collapse duplimotifs and get motifs

for(a in unique(mirna_data$ACCESSION)){

  duplimirna <- unique(mirna_data[which(mirna_data$ACCESSION == a), 'DUPLI.ID'])

  if(nchar(duplimirna) < 1){ next }

  seq <- unique(mirna_data[which(mirna_data$ACCESSION==a), 'SEQUENCE'])
  duplimirna <- strsplit(duplimirna, ',', fixed=T)[[1]]

  for(d in duplimirna){
    acc <- unique(mirna_data[which(mirna_data$MIRNA==d), 'ACCESSION'])
    curr_seq <- unique(mirna_data[which(mirna_data$ACCESSION==acc), 'SEQUENCE'])
    if(nchar(curr_seq) > nchar(seq)){
      unique_motif_data <- unique_motif_data[unique_motif_data$ACCESSION!=acc,]
    }

  }


}


for(a in unique_motif_data$ACCESSION){
  
  motif <- unique_motif_data[which(unique_motif_data$ACCESSION==a), 'UNIQUE.MOTIF']
  if(nchar(motif) > 0){ next }
    
  seq <- unique_motif_data[which(unique_motif_data$ACCESSION==a), 'SEQUENCE']
  sequences <- unique(unique_motif_data[which(unique_motif_data$SEQUENCE!=seq), 'SEQUENCE'])
  len <- 7
  motif <- getMotif(len, seq)
  while(!isUniqueMotif(motif, sequences)){
    len <- len + 1
    motif <- getMotif(len, seq)
  }
  unique_motif_data[which(unique_motif_data$ACCESSION==a), 'UNIQUE.MOTIF'] <- motif
}

unique_motif_data$MOTIF.LEN <- nchar(unique_motif_data$UNIQUE.MOTIF)
motif_len <- unique_motif_data$MOTIF.LEN
motif_len <- aggregate(motif_len, by=list(motif_len), FUN=sum)
colnames(motif_len) <- c('MOTIF.LEN', 'COUNT')
motif_len$COUNT <- motif_len$COUNT/motif_len$MOTIF.LEN
motif_len$FREQUENCY <- motif_len$COUNT/sum(motif_len$COUNT)
motif_len <- motif_len[with(motif_len, order(MOTIF.LEN)), ]


for(i in 1:nrow(motif_len)){
  motif_len[i, 'CUMULATIVE.FREQ'] <- motif_len[i, 'FREQUENCY']
  if( i == 1 ){ next }
  motif_len[i, 'CUMULATIVE.FREQ'] <- motif_len[i, 'FREQUENCY'] + motif_len[i-1, 'CUMULATIVE.FREQ']
}

plot(motif_len[,c('MOTIF.LEN', 'CUMULATIVE.FREQ')])
lines(motif_len[,c('MOTIF.LEN', 'CUMULATIVE.FREQ')], type='o') 
