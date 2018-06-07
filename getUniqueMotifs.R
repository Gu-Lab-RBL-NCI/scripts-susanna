
## PURPOSE: generate unique motifs from miRBase data
## INPUT: miRBase data  miRNA.xls
## OUTPUT: fasta file containing unique motifs all-unique-motifs.fa
require(gdata)
require(plyr)

base_path<-'/Users/chens22/Documents/'
miRBase_file <- '/Users/chens22/Documents/miRNA.xls'
miRBase_data <- read.xls(miRBase_file, sheet = 1, header = TRUE)

motif_len <- 8

getMirnaData <- function(miRBase_file){
  miRBase_data <- read.xls(miRBase_file, sheet = 1, header = TRUE)
  miRBase_data <- na.omit(miRBase_data)
  miRBase_data <- Filter(function(x)!all(is.na(x)), miRBase_data)
  
  names(miRBase_data)[names(miRBase_data) == 'Accession'] <- 'PRI.ACCESSION'
  names(miRBase_data)[names(miRBase_data) == 'Sequence'] <- 'PRI.SEQUENCE'
  names(miRBase_data)[names(miRBase_data) == 'ID'] <- 'PRIMIRNA'
  
  mature2_names <- c('PRI.ACCESSION', 'PRIMIRNA', 'PRI.SEQUENCE', 'Mature2_Acc', 'Mature2_ID', 'Mature2_Seq')
  mature2 <- miRBase_data[, mature2_names]
  miRBase_data <- miRBase_data[,!names(miRBase_data) %in% c('Mature2_Acc', 'Mature2_ID', 'Mature2_Seq')]
  names(mature2)[names(mature2) == 'Mature2_Acc'] <- 'ACCESSION'
  names(mature2)[names(mature2) == 'Mature2_ID'] <- 'MIRNA'
  names(mature2)[names(mature2) == 'Mature2_Seq'] <- 'SEQUENCE'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_Acc'] <- 'ACCESSION'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_ID'] <- 'MIRNA'
  names(miRBase_data)[names(miRBase_data) == 'Mature1_Seq'] <- 'SEQUENCE'
  miRBase_data <- rbind.fill(miRBase_data, mature2)
  
  miRBase_data$SEQUENCE<-gsub('U', 'T', miRBase_data$SEQUENCE)
  miRBase_data$PRI.SEQUENCE<-gsub('U', 'T', miRBase_data$PRI.SEQUENCE)
  
  miRBase_data <- miRBase_data[, -grep("X", colnames(miRBase_data), fixed=T)]
  
  return(miRBase_data)
}

getDupliMirna<-function(mirna){
  sequence<-mirna_data[which(mirna_data$MIRNA==mirna), 'SEQUENCE'][1]
  unique_mirna<-unique(mirna_data$MIRNA)
  dupli_id<-c()
  for(m in unique_mirna){
    if(!identical(toString(m), toString(mirna))){
      curr_seq<-mirna_data[which(mirna_data$MIRNA==m), 'SEQUENCE']
      for(c in curr_seq){
        if(gregexpr(sequence, c, fixed=TRUE)[[1]]!=-1 || gregexpr(c, sequence, fixed=TRUE)[[1]]!=-1){
          dupli_id<-c(dupli_id, m)
        }
      }
    }
  }
  return(unique(dupli_id))
}

getUniqueMotif<-function(length, seq, primary_sequences){
  mid<-trunc(nchar(seq)/2)
  start<-mid-(trunc(length/2))
  end<-mid+(trunc(length/2))
  motif<-substr(seq, start, end)
  if(!inPriSeq(motif, primary_sequences)){
    return(motif)
  }else{
    return(NA)
  }
  
}

getShiftedMotif<-function(length, seq, primary_sequences){
  mid<-trunc(nchar(seq)/2)
  start<-mid-(trunc(length/2))
  end<-mid+(trunc(length/2))
  motif<-substr(seq, start, end)
  shift<-1
  while(inPriSeq(motif, primary_sequences) & shift<trunc(nchar(seq)/2)){
    mid<-mid+shift
    start<-mid-(trunc(length/2))
    end<-mid+(trunc(length/2))
    motif<-substr(seq, start, end)
    if(inPriSeq(motif, primary_sequences)){
      mid_5p<-mid+(shift*-2)
      start<-mid-(trunc(length/2))
      end<-mid+(trunc(length/2))
      motif<-substr(seq, start, end)
      if(!inPriSeq(motif, primary_sequences)){
        return(motif)
      }
    }
    shift<-shift+1
  }
  if(nchar(motif)==length){
    return(motif)
  }else{
    return(NA)
  }
  
}

getExtendedMotif<-function(length, seq, primary_sequences){

  mid<-trunc(nchar(seq)/2)
  start<-mid-(trunc(length/2))
  end<-mid+(trunc(length/2))
  motif<-substr(seq, start, end)

  temp_end<-end
  while(inPriSeq(motif, primary_sequences)){
    temp_end<-temp_end+1
    if(temp_end>length(seq)){
      motif<-''
      break
    }
    motif<-substr(seq, start, temp_end)
  }
  motif_3p<-motif
  
  motif<-substr(seq, start, end)
  
  temp_start<-start
  while(inPriSeq(motif, primary_sequences)){
    temp_start<-temp_start-1
    if(temp_start<1){
      motif<-''
      break
    }
    motif<-substr(seq, temp_start, end)
  }
  motif_5p<-motif
  
  if(nchar(motif_5p)>nchar(motif_3p)){
    return(motif_5p)
  }else if(nchar(motif_5p)<nchar(motif_3p)){
    return(motif_3p)
  }else if(nchar(motif_5p)>0){
    return(motif_5p)
  }else{
    motif<-getTwoEndExtendedMotif(motif_len, seq, primary_sequences)
    return(motif)
  }
  return(NA)
}

getTwoEndExtendedMotif<-function(length, seq, primary_sequences){
  
  mid<-trunc(nchar(seq)/2)
  start<-mid-(trunc(length/2))
  end<-mid+(round(length/2))
  motif<-substr(seq, start, end)

  while(inPriSeq(motif, primary_sequences)){
    length<-length+1
    if(length%%2==0){
      start<-start-1
    }else{
      end<-end+1
    }
    if(start==1 && end==nchar(seq)){
      return(seq)
    }else if(start<1 && end>nchar(seq)){
      return(NA)
    }
    motif<-substr(seq, start, end)
  }
  return(motif)
}

inPriSeq<-function(seq, primary_sequences){
  for(i in 1:length(primary_sequences)){
    p<-primary_sequences[i]
    found<-gregexpr(seq, p, fixed=TRUE)[[1]]
    if(found!=-1){
      return(TRUE)
    }
  }
  return(FALSE)
}

getMotifPosition<-function(motif, primary_seq){
  
  start<-gregexpr(motif, primary_seq, fixed=TRUE)[[1]]
  end<-start+nchar(motif)
  return(c(start, end))
  
}

getPrimarySequences<-function(mirna){
  accession<-mirna_data[which(mirna_data$MIRNA==mirna),'ACCESSION']
  pri_accession<-mirna_data[which(mirna_data$MIRNA==mirna),'PRI.ACCESSION']
  dupli_id<-as.vector(strsplit(mirna_data[which(mirna_data$MIRNA==mirna), 'DUPLI.ID'], ',', fixed=TRUE)[[1]])
  primary_seq<-c()
  #get one primary sequence per accession
  for(i in 1:nrow(mirna_data)){
    curr_acc<-mirna_data[i,'ACCESSION']
    curr_pri_acc<-mirna_data[i,'PRI.ACCESSION']
    curr_mirna<-mirna_data[i, 'MIRNA']
    curr_seq<-mirna_data[i,'SEQUENCE']
    if(!(curr_acc%in%accession) && !(curr_pri_acc%in%pri_accession) && !(curr_mirna%in%dupli_id)){
      primary_seq<-c(primary_seq, curr_seq)
    }
  }
  primary_seq<-unique(primary_seq)
  primary_seq<-primary_seq[!is.na(primary_seq)]
  return(primary_seq)
}

mirna_data <- getMirnaData(miRBase_file)

mirna_data <- mirna_data[mirna_data$PRIMIRNA %like% "hsa", ]

mirna_data <- mirna_data[1:100,]

unique_mirna<-unique(mirna_data$MIRNA)
unique_mirna<-unique_mirna[!is.na(unique_mirna)]
for(m in unique_mirna){
  mirna_data[which(mirna_data$MIRNA==m),'DUPLI.ID']<-paste(getDupliMirna(m), collapse=',')
}


for(m in unique_mirna){
  mirna_data[which(mirna_data$MIRNA==m),'DUPLI.ID']<-paste(getDupliMirna(m), collapse=',')
  temp_data<-mirna_data[which(mirna_data$MIRNA==m),]
  curr_seq<-temp_data[1, 'SEQUENCE']
  pri_seq<-getPrimarySequences(m)
  
  unique_motif<-getUniqueMotif(motif_len, curr_seq, pri_seq)
  mirna_data[which(mirna_data$MIRNA==m),paste0('MOTIF.', motif_len)]<-unique_motif
  if(!is.na(unique_motif)){
    motif_position<-getMotifPosition(unique_motif, curr_seq)
    new_mirna_id<-paste(m, '[', motif_position[1], ',', motif_position[2], ']', sep='')
    mirna_data[which(mirna_data$MIRNA==m),'POS.MIRNA']<-new_mirna_id
  }
}

mirna_data$DUPLI.MIRNA<-NA

for(i in 1:nrow(mirna_data)){
  dupli_id<-mirna_data[i, 'DUPLI.ID']
  dupli_id<-strsplit(dupli_id,',', fixed=TRUE)[[1]]
  mirna<-mirna_data[i, 'MIRNA']
  motif<-mirna_data[i, 'MOTIF.13']
  if(length(dupli_id>0)){
    new_mirna_id<-c(mirna)
    new_mirna_id<-paste(new_mirna_id, collapse=',')
    for(d in dupli_id){
      temp_dupli<-strsplit(d, '-', fixed=TRUE)[[1]]
      temp_dupli<-paste(temp_dupli[3:length(temp_dupli)], collapse='-')
      new_mirna_id<-paste(new_mirna_id, temp_dupli, sep='-', collapse='-')
    }
    for(d in dupli_id){
      if(is.na(mirna_data[which(mirna_data$MIRNA==d),'DUPLI.MIRNA'])){
        mirna_data[which(mirna_data$MIRNA==d),'DUPLI.MIRNA']<-new_mirna_id
        mirna_data[which(mirna_data$MIRNA==d),'MOTIF.13']<-motif
        curr_seq<-mirna_data[which(mirna_data$MIRNA==d),'SEQUENCE']
        motif_start<-gregexpr(motif, curr_seq, fixed=TRUE)[[1]][1]
        motif_end<-motif_start+nchar(motif)
        mirna_data[which(mirna_data$MIRNA==d),'POS.MIRNA']<-paste(new_mirna_id, '[', motif_start, ',', motif_end, ']', sep='')
      }
    }
    if(is.na(mirna_data[which(mirna_data$MIRNA==mirna),'DUPLI.MIRNA'])){

      mirna_data[which(mirna_data$MIRNA==mirna), 'DUPLI.MIRNA']<-new_mirna_id
      curr_seq<-mirna_data[i,'SEQUENCE']
      motif_start<-gregexpr(motif, curr_seq, fixed=TRUE)[[1]][1]
      motif_end<-motif_start+nchar(motif)
      mirna_data[which(mirna_data$MIRNA==mirna),'POS.MIRNA']<-paste(new_mirna_id, '[', motif_start, ',', motif_end, ']', sep='')
    }
  }
}



#generate fasta file with mirna id w/ position and new unique motifs

#motif_file <- '/Volumes/2TB (MAC)/Susanna/miRBase21/all-unique-motifs.fa'
#motif_data <- readLines(motif_file)
data<-c()
unique_mirna<-unique(mirna_data$POS.MIRNA)
unique_mirna <- na.omit(unique_mirna)
for(m in unique_mirna){
  temp_data<-mirna_data[which(mirna_data$POS.MIRNA==m),]
  temp_data<-temp_data[1,]
  sequence<-temp_data[,'SEQUENCE']
  motif<-temp_data[,'MOTIF.13']
  #if(length(grep(toString(paste0('\t', motif)), motif_data, value = TRUE))<1){
   # print(m)
   # print(motif)
  #}
  
  first_line<-paste(paste('>', toString(m), sep=''), toString(motif), sep='\t')
  second_line<-toString(sequence)
  
  data<-c(data, first_line, second_line)
}

motifs_file<-paste(base_path, 'all-unique-motifs.fa', sep='')
file<-file(motifs_file)
writeLines(data, file)
close(file)
 