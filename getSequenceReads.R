
## PURPOSE: acquire reads of miRNA by filtering through collapsed files with a motif
## INPUT: miRBase data                                 miRBAse21-master.tsv
##        manifest data                              all-tumor-manifest.tsv
##        collapsed files         sample.converted.unpaired.fastq.collapsed
## OUTPUT: reads of primiRNA                        motif.primirna.reads.tsv
##         motifs generated from TCGA data  TCGA.motifs.include.paralogs.tsv

require('stringr')
require('plyr')
require('stringdist')

base_path <- '/Volumes/2TB (MAC)/Susanna/'#'/media/user/2TB (MAC)/Susanna/'
#base_path <- '/Users/chens22/Documents/'

doc_path <- '/Users/chens22/Documents/seq_reads/'

mirna_file <- paste0(base_path, 'miRBase21/miRBase21-master.tsv')
mirna_data <- read.table(file=mirna_file, sep='\t', header=T, stringsAsFactors = F)

manifest_file <- paste0(base_path, 'all-tumor-manifest.csv')
manifest_data <- read.table(file=manifest_file, sep=',', stringsAsFactors=F, header=T)

getManifestID <- function(name, tumor){
  
  return(toString(unique(manifest_data[which(manifest_data$NAME == name & manifest_data$DISEASE.ABBV == tumor), 'ID'][1])))
  
}


get3PMirna<-function(mirna){
  pri_mirna<-unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA'])
  mirna_3p<-c()
  for(p in pri_mirna){
    mirna_3p<-c(mirna_3p, toString(unique(mirna_data[which(mirna_data$PRIMIRNA==p & mirna_data$STRAND=='3P'), 'MIRNA'])))
  }
  if(length(mirna_3p)>0){
    return(toString(unique(mirna_3p)[1]))
  }else{
    print('ERROR ERROR ERROR ERROR ERROR')
    return('')
  }
}

get5PMirna<-function(mirna){
  pri_mirna<-unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA'])
  mirna_5p<-c()
  for(p in pri_mirna){
    mirna_5p<-c(mirna_5p, toString(unique(mirna_data[which(mirna_data$PRIMIRNA==p & mirna_data$STRAND=='5P'), 'MIRNA'])))
  }
  if(length(mirna_5p)>0){
    return(toString(unique(mirna_5p)[1]))
  }else{
    print('ERROR ERROR ERROR ERROR ERROR')
    return('')
  }
}

get5PSequence<-function(primirna){
  seq_5p<-mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='5P'), 'SEQUENCE']
  if(length(seq_5p)>1){
    print(paste('ERROR SEQ 5P', primirna))
    return('')
  }
  return(toString(seq_5p))
}

get3PSequence<-function(primirna){
  seq_3p<-mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='3P'), 'SEQUENCE']
  if(length(seq_3p)>1){
    print(paste('ERROR SEQ 3P', primirna))
    return('')
  }
  return(toString(seq_3p))
}

getSecondaryStructure<-function(sequence){
  input_data<-c(paste('> ', sep=''), sequence)
  input_file<-paste(doc_path, 'temp-miRNA-seq.fa', sep='')
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-paste(doc_path, 'temp-miRNA-secondary-struct.txt', sep='')
  secondary_structure<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  secondary_structure<-strsplit(secondary_structure, ' ', fixed=T)[[1]][1]
  return(toString(secondary_structure))
}

getPriMirna<-function(mirna){
  return(unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA']))
}

getPrimarySequence<-function(mirna){
  pri_seqs<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(mirna) | mirna_data$MIRNA==mirna), 'PRI.SEQUENCE'])
  extended_seq<-getExtendedSequence(mirna)
  for(p in pri_seqs){
    if(gregexpr(toString(p), toString(extended_seq))[[1]]!=-1){
      return(toString(p))
    }
  }
}

getPrecursorSequence<-function(primirna){
  primirna<-getPriMirna(primirna)
  pri_seq<-getPrimarySequence(primirna)
  pri_secondary_struct<-getSecondaryStructure(pri_seq)
  seq_5p<-get5PSequence(primirna)
  seq_3p<-get3PSequence(primirna)
  if(nchar(seq_5p)>0){
    start<-gregexpr(seq_5p, pri_seq)[[1]][1]
  }else{
    end<-gregexpr(seq_3p, pri_seq)[[1]]+nchar(seq_3p)-1 
    temp_pri_struct<-substr(pri_secondary_struct, 1, end)
    precurs_brackets_3p<-str_count(temp_pri_struct, '\\)')
    mor_5p_brackets<-str_count(pri_secondary_struct, '\\(')
    precursor_start<-mor_5p_brackets-precurs_brackets_3p+1
    start<-gregexpr('(', pri_secondary_struct, fixed=T)[[1]][precursor_start]
  }
  
  if(nchar(seq_3p)>0){
    end<-gregexpr(seq_3p, pri_seq)[[1]]+nchar(seq_3p)-1
  }else{
    start<-gregexpr(seq_5p, pri_seq)[[1]]
    temp_pri_struct<-substr(pri_secondary_struct, start, nchar(pri_secondary_struct))
    precurs_brackets_5p<-str_count(temp_pri_struct, '\\(')
    end<-gregexpr(')', pri_secondary_struct, fixed=T)[[1]][precurs_brackets_5p]
  }
  
  if(is.na(start) || is.na(end) || length(start)<1 || length(end)<1){
    return('')
  }
  pre_seq<-substr(pri_seq, start, end)
  return(pre_seq)
}

getExtendedSequence<-function(mirna){
  extended_seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(mirna)), 'EXTENDED.SEQUENCE'])
  if(length(extended_seq)>1){
    print(paste('ERROR EXTENDED SEQ', mirna))
    return(NA)
  }
  return(c(toString(extended_seq)))
}

getEditedExtendedSequence <- function(primirna){
  pre_seq <- getPrecursorSequence(primirna)
  if(nchar(pre_seq) < 1) {return('')}
  extended_seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(primirna)), 'EXTENDED.SEQUENCE'])
  pre_start <- gregexpr(pre_seq, extended_seq, fixed=T)[[1]][1]
  pre_end <- pre_start + nchar(pre_seq) - 1
  return(substr(extended_seq, pre_start-30, pre_end+30))
}

getPredictedMatureSequence <- function(motif, pri_seq){
  
  len <- 25
  
  motif_start <- gregexpr(toString(motif), toString(pri_seq), fixed=T)[[1]][1]
  motif_len <- as.integer(nchar(motif)/2)
  
  center <- motif_start + motif_len
  
  mature_start <- center - as.integer(len/2)
  mature_end <- mature_start + len - 1
  
  mature_seq <- substr(pri_seq, mature_start, mature_end)
  return(mature_seq)
}

isUniqueMotif <- function(primirna, motif){
  accession <- mirna_data[which(mirna_data$PRIMIRNA==primirna), 'PRI.ACCESSION']
  sequences <- mirna_data[which(!mirna_data$PRI.ACCESSION%in%accession),'PRI.SEQUENCE']
  for(s in sequences){
    if(gregexpr(toString(motif), toString(s), fixed=T)[[1]][1] != -1){ return(F) }
  }
  return(T)
}

get3PMotif <- function(motif, primirna){
  pri_seq <- getPrimarySequence(primirna)
  if(nchar(pri_seq) < 1){return('')}
  pri_struct <- getSecondaryStructure(pri_seq)
  
  seq_5p <- get5PSequence(primirna)
  seq_start <- gregexpr(seq_5p, pri_seq, fixed=T)[[1]][1]
  seq_end <- seq_start + nchar(seq_5p) - 1
  
  motif_start <- gregexpr(motif, pri_seq, fixed=T)[[1]][1]
  motif_end <- motif_start + nchar(motif) - 1
  motif_struct <- substr(pri_struct, motif_start, motif_end)
  
  if(gregexpr(')', motif_struct, fixed=T)[[1]][1] != -1){
    return('')
  }
  
  start_3p <- seq_end + gregexpr(')', substr(pri_struct, seq_end + 1, nchar(pri_struct)), fixed=T)[[1]][1]
  between_paired <- str_count(substr(pri_struct, motif_end + 1, start_3p), '\\(')
  if(between_paired<1){
    motif_start_3p <- start_3p
  }else{
    motif_start_3p <- start_3p + gregexpr(')', substr(pri_struct, start_3p+1, nchar(pri_struct)), fixed=T)[[1]][between_paired]
  }
  motif_paired <- str_count(motif_struct, '\\(')
  motif_end_3p <- motif_start_3p + gregexpr(')', substr(pri_struct, motif_start_3p, nchar(pri_struct)), fixed=T)[[1]][motif_paired]
  
  motif_3p <- substr(pri_seq, motif_start_3p, motif_end_3p)
  motif_len <- motif_end_3p - motif_start_3p + 1
  motif_center <- motif_start_3p + as.integer(motif_len/2)
  while(isUniqueMotif(primirna, motif_3p) & motif_len > 13){
    motif_len <- motif_len - 1
    motif_center <- motif_start_3p + as.integer(motif_len/2)
    motif_start_3p <- motif_center - as.integer(motif_len/2)
    motif_end_3p <- motif_start_3p + motif_len - 1
    motif_3p <- substr(pri_seq, motif_start_3p, motif_end_3p)
  }
  
  while((!isUniqueMotif(primirna, motif_3p) & motif_len < nchar(seq_5p)) || motif_len < 13){
    motif_len <- motif_len + 1
    motif_center <- motif_start_3p + as.integer(motif_len/2)
    motif_start_3p <- motif_center - as.integer(motif_len/2)
    motif_end_3p <- motif_start_3p + motif_len - 1
    motif_3p <- substr(pri_seq, motif_start_3p, motif_end_3p)
  }
  return(motif_3p)
}

get5PMotif <- function(motif, primirna){
  pri_seq <- getPrimarySequence(primirna)
  if(nchar(pri_seq) < 1){return('')}
  pri_struct <- getSecondaryStructure(pri_seq)
  
  seq_3p <- get3PSequence(primirna)
  seq_start <- gregexpr(seq_3p, pri_seq, fixed=T)[[1]][1]
  
  motif_start <- gregexpr(motif, pri_seq, fixed=T)[[1]][1]
  motif_end <- motif_start + nchar(motif) - 1
  motif_struct <- substr(pri_struct, motif_start, motif_end)
  if(gregexpr('(', motif_struct, fixed=T)[[1]][1] != -1){
    return('')
  }

  
  end_5p <- gregexpr('(', substr(pri_struct, 1, seq_start - 1), fixed=T)[[1]]
  end_5p <- end_5p[length(end_5p)]
  between_paired <- str_count(substr(pri_struct, end_5p, motif_start-1), '\\)')
  motif_end_5p <- gregexpr('(', substr(pri_struct, 1, end_5p), fixed=T)[[1]]
  motif_end_5p <- motif_end_5p[length(motif_end_5p) - between_paired]
  motif_paired <- str_count(motif_struct, '\\)')
  motif_start_5p <- gregexpr('(', substr(pri_struct, 1, motif_end_5p), fixed=T)[[1]]
  motif_start_5p <- motif_start_5p[(length(motif_start_5p) - motif_paired + 1)]
  
  
  
  motif_5p <- substr(pri_seq, motif_start_5p, motif_end_5p)
  motif_len <- motif_end_5p - motif_start_5p + 1
  motif_center <- motif_start_5p + as.integer(motif_len/2)
  while(isUniqueMotif(primirna, motif_5p) & motif_len > 13){
    motif_len <- motif_len - 1
    motif_center <- motif_start_5p + as.integer(motif_len/2)
    motif_start_5p <- motif_center - as.integer(motif_len/2)
    motif_end_5p <- motif_start_5p + motif_len - 1
    motif_5p <- substr(pri_seq, motif_start_5p, motif_end_5p)
  }
  
  while((!isUniqueMotif(primirna, motif_5p)& motif_len < nchar(seq_3p)) || motif_len < 13){
    motif_len <- motif_len + 1
    motif_center <- motif_start_5p + as.integer(motif_len/2)
    motif_start_5p <- motif_center - as.integer(motif_len/2)
    motif_end_5p <- motif_start_5p + motif_len - 1
    motif_5p <- substr(pri_seq, motif_start_5p, motif_end_5p)
  }
  return(motif_5p)
}

unique_primirna <- unique(mirna_data$PRIMIRNA)

motifs_data <- data.frame(PRIMIRNA=NA, MOTIF=NA, STRAND=NA, EXISTS=NA, OTHER.STEM=NA, UNIQUE=NA)
for(p in unique_primirna){

  motif_5p <- toString(mirna_data[which(mirna_data$PRIMIRNA==p & mirna_data$STRAND=='5P'), 'MOTIF.13'])
  motif_3p <- toString(mirna_data[which(mirna_data$PRIMIRNA==p & mirna_data$STRAND=='3P'), 'MOTIF.13'])

  other_5p <- F
  exists_5p <- T
  if(nchar(motif_5p) < 1){
    motif_5p <- get5PMotif(motif_3p, p)
    if(nchar(motif_5p) < 1){other_5p <- T}
    exists_5p <- F
  }
  other_3p<- F
  exists_3p <- T
  if(nchar(motif_3p) < 1){
    motif_3p <- get3PMotif(motif_5p, p)
    if(nchar(motif_3p) < 1){other_3p <- T}
    exists_3p <- F
  }

  unique_5p <- isUniqueMotif(p, motif_5p)
  unique_3p <- isUniqueMotif(p, motif_3p)
  motifs_data<- rbind(motifs_data, data.frame(PRIMIRNA=p, MOTIF=c(motif_5p, motif_3p), EXISTS=c(exists_5p, exists_3p), STRAND=c('5P', '3P'), OTHER.STEM=c(other_5p, other_3p), UNIQUE=c(unique_5p, unique_3p)))

  #write.table(data.frame(SEQUENCE=NA), file=paste0(doc_path, motif_5p, '.', p, '.reads.tsv'), sep='\t', row.names=F)
  #write.table(data.frame(SEQUENCE=NA), file=paste0(doc_path, motif_3p, '.', p, '.reads.tsv'), sep='\t', row.names=F)

}
motifs_data <- motifs_data[-1,]
write.table(motifs_data, file=paste0(doc_path, 'TCGA.motifs.include.paralogs.tsv'), sep='\t', row.names=F)

motifs_data <- read.table(paste0(doc_path, 'TCGA.motifs.include.paralogs.tsv'), sep='\t', header=T, stringsAsFactors = F)
motifs_data$PREDICTED.SEQUENCE <- NA
for(i in 1:nrow(motifs_data)){
  motif <- toString(motifs_data[i, 'MOTIF'])
  primirna <- toString(motifs_data[i, 'PRIMIRNA'])
  pri_seq <- getPrimarySequence(primirna)
  predicted_seq <- getPredictedMatureSequence(motif, pri_seq)
  motifs_data[i, 'PREDICTED.SEQUENCE'] <- predicted_seq
  motifs_data[i, 'PRI.SEQUENCE'] <- pri_seq
  motifs_data[i, 'MOTIF.LEN'] <- nchar(motifs_data[i, 'MOTIF'])
}

write.table(motifs_data, file=paste0(doc_path, 'TCGA.motifs.tsv'), sep='\t', row.names=F)

motifs_data <- read.table(paste0(doc_path, 'TCGA.motifs.tsv'), sep='\t', header=T, stringsAsFactors = F)
motifs_data<-motifs_data[which(motifs_data$EXISTS==F & motifs_data$UNIQUE==T & motifs_data$MOTIF.LEN>12),]
motifs_data <- motifs_data[complete.cases(motifs_data), ]

count <- 0
for( d in list.files(path = base_path, full.names = T, include.dirs = T, all.files = T) ) {

  collapsed_dir <- paste0(d, '/collapsed_fastq/')
  if( file.exists(collapsed_dir) ){
    
    tumor <- strsplit(toString(d), '/', fixed=T)[[1]]
    tumor <- toString(tumor[length(tumor)])
    if(tumor%in%c('ACC', 'BLCA', 'BRCA', 'COAD', 'CESC', 'CHOL', 'DLBC', 'ESCA', 'GBM')){next}
    print(tumor)
    
    max_lines <- 0
    start_time <- as.numeric(Sys.time())
    for( c in list.files(path = collapsed_dir, full.names=T, all.files=T, pattern = '*.collapsed') ){
      count <- count + 1
      if(count < 500){next}
      name <- strsplit(c, '\\.')[[1]][1]
      name <- strsplit(name, '\\/')[[1]]
      name <- toString(name[length(name)])
      
      print(name)
      
      id <- getManifestID(name, tumor)
      
      collapsed_data <- read.table(c, header=F, colClasses=NA, sep='', stringsAsFactors = F, comment.char='', skipNul=T, nrows=10000000, fill=T)
      collapsed_data <- na.omit(collapsed_data)
      colnames(collapsed_data) <- c(id, 'SEQUENCE')
      collapsed_data <- collapsed_data[which(collapsed_data[,id] > 2),]
      
      for(m in 1:nrow(motifs_data)){
        #print(m)
        motif <- toString(motifs_data[m, 'MOTIF'])
        primirna<- toString(motifs_data[m, 'PRIMIRNA'])
        predicted_seq <- toString(motifs_data[m, 'PREDICTED.SEQUENCE'])
        
        motif_file <- paste0(doc_path, motif, '.', primirna, '.reads.tsv')
        motif_data <- data.frame(SEQUENCE=NA, stringsAsFactors = F)
        if(file.exists(motif_file)){
          motif_data <- read.table(motif_file, sep='\t', header=T, stringsAsFactors = F, check.names = F)
        }
        
        if(id %in% colnames(motif_data)){next}
        
        temp_collapsed_data <- collapsed_data[grepl(motif, collapsed_data$SEQUENCE),]
        motif_data <- motif_data[!is.na(colnames(motif_data))]
        motif_data <- motif_data[!grepl('NA', names(motif_data), fixed=T)]
        if(nrow(motif_data) < 1 && !'SEQUENCE'%in%names(motif_data)){
          motif_data <- data.frame(SEQUENCE=NA)
        }else{
          #motif_data <- motif_data[,colSums(is.na(motif_data))<nrow(motif_data)]
        }
        motif_data <- merge(motif_data, temp_collapsed_data, by='SEQUENCE', all=T, na.rm=T)
        
        motif_data <- motif_data[,order(names(motif_data), decreasing=T)]
        motif_data <- motif_data[!is.na(motif_data$SEQUENCE),]
        write.table(motif_data, file=motif_file, sep='\t', row.names=F)
      }
      
      
      collapsed_file <- file(c, open="rb")
      nlines <- 0L
      while (length(chunk <- readBin(collapsed_file, "raw", 65536)) > 0) {
        nlines <- nlines + sum(chunk == as.raw(10L))
      }
      close(collapsed_file)
      
      if(nlines > max_lines){
        max_collapsed_file <- c
        max_lines <- nlines
      }
    }
    
    name <- strsplit(max_collapsed_file, '\\.')[[1]][1]
    name <- strsplit(name, '\\/')[[1]]
    name <- toString(name[length(name)])
    
    print(name)
    
    id <- getManifestID(name, tumor)
    
    collapsed_data <- read.table(max_collapsed_file, header=F, colClasses=NA, sep='', stringsAsFactors = F, comment.char='', skipNul=T, nrows=10000000, fill=T)
    collapsed_data <- na.omit(collapsed_data)
    colnames(collapsed_data) <- c(id, 'SEQUENCE')
    collapsed_data <- collapsed_data[which(collapsed_data[,id] > 2),]
    
    for(m in 1:nrow(motifs_data)){
      #print(m)
      motif <- toString(motifs_data[m, 'MOTIF'])
      primirna<- toString(motifs_data[m, 'PRIMIRNA'])
      predicted_seq <- toString(motifs_data[m, 'PREDICTED.SEQUENCE'])
      
      motif_file <- paste0(doc_path, motif, '.', primirna, '.reads.tsv')
      motif_data <- data.frame(SEQUENCE=NA, stringsAsFactors = F)
      if(file.exists(motif_file)){
        motif_data <- read.table(motif_file, sep='\t', header=T, stringsAsFactors = F, check.names = F)
      }
      
      if(id %in% colnames(motif_data)){next}
      
      temp_collapsed_data <- collapsed_data[grepl(motif, collapsed_data$SEQUENCE),]
      motif_data <- motif_data[!is.na(names(motif_data))]
      if(nrow(motif_data) < 1){
        motif_data[1,'SEQUENCE'] <- NA
      }else{
        motif_data <- motif_data[,colSums(is.na(motif_data))<nrow(motif_data)]
      }
      motif_data <- merge(motif_data, temp_collapsed_data, by='SEQUENCE', all=T, na.rm=T)
      
      motif_data <- motif_data[,order(names(motif_data), decreasing=T)]
      motif_data <- motif_data[!is.na(motif_data$SEQUENCE),]
      write.table(motif_data, file=motif_file, sep='\t', row.names=F)
    }
    
    print((as.numeric(Sys.time()) - start_time)/60)
  }
}

for(m in 1:nrow(motifs_data)){

  motif <- toString(motifs_data[m, 'MOTIF'])
  primirna<- toString(motifs_data[m, 'PRIMIRNA'])
  predicted_seq <- motifs_data[m, 'PREDICTED.SEQUENCE']
  pri_seq <- motifs_data[m, 'PRI.SEQUENCE']

  motif_file <- paste0(doc_path, motif, '.', primirna, '.reads.tsv')
  motif_data <- read.table(motif_file, sep='\t', header=T, stringsAsFactors = F, check.names = F)
  motif_data[,'SUM'] <- rowSums(motif_data[,!colnames(motif_data)%in%c('SUM', 'SEQUENCE', 'DISTANCE.FROM.PREDICTED', 'TEMPLATED'), drop=FALSE], na.rm=T)
  for(i in 1:nrow(motif_data)){
    seq <- toString(motif_data[i, 'SEQUENCE'])
    motif_data[i, 'DISTANCE.FROM.PREDICTED'] <- stringdist(seq, predicted_seq)
    motif_data[i, 'TEMPLATED'] <- grepl(seq, pri_seq)
  }
  motif_data <- motif_data[,order(names(motif_data), decreasing=T)]


  write.table(motif_data, file=motif_file, sep='\t', row.names=F)
}
