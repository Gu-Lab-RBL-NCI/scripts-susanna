
## PURPOSE: get unique motif starting from 5P end of strand
## INPUT: miRBase data                                    miRBase21-master.tsv
## OUTPUT: table containing unique 5p motifs    miRBase21-unique-5p-motifs.tsv

base_path<-'/Users/chens22/Documents/miRBase21/'

mirna_file <- paste0(base_path, 'miRBase21-master.tsv')
mirna_data <- read.table(file=mirna_file, sep='\t', header=T, stringsAsFactors=F)

get5PMotif <- function(sequence, len){
  motif <- substr(sequence, 1, len)
  return(motif)
}

isUniqueMotif <- function(motif, acc, duplimirna){
  sequences <- mirna_data[which(mirna_data$ACCESSION != acc & !mirna_data$MIRNA%in%getDuplimirna(acc) & !mirna_data$MIRNA%in%duplimirna), 'SEQUENCE']
  for(s in sequences){
    if(gregexpr(toString(motif), toString(s), fixed=T) != -1){
      print(s)
      return(FALSE)
    }
  }
  return(TRUE)
}

getMirna <- function(acc){
  primirna <- toString(unique(mirna_data[which(mirna_data$ACCESSION==acc), 'MIRNA']))
  return(primirna)
}

getPrimirna <- function(acc){
  primirna <- toString(unique(mirna_data[which(mirna_data$ACCESSION==acc), 'PRIMIRNA']))
  return(primirna)
}

getSequence <- function(acc){
  seq <- toString(unique(mirna_data[which(mirna_data$ACCESSION==acc), 'SEQUENCE']))
  return(seq)
}

getAccession <- function(mirna){
  seq <- toString(unique(mirna_data[which(mirna_data$MIRNA==mirna), 'ACCESSION']))
  return(seq)
}

isDuplimotif <- function(acc){
  dupli <- toString(unique(mirna_data[which(mirna_data$ACCESSION==acc), 'DUPLIMOTIF']))
  return(dupli)
}

getDuplimirna <- function(acc){
  dupli <- toString(unique(mirna_data[which(mirna_data$ACCESSION==acc), 'DUPLI.ID']))
  return(strsplit(dupli, ',', fixed=T)[[1]])
}

unique_acc <- unique(mirna_data$ACCESSION)
unique_motifs <- data.frame(MIRNA=character(), MOTIF=character(), END = c(), DUPLIMOTIF = c(), stringsAsFactors = F)
for(a in unique_acc){
  mirna <- getMirna(a)
  primirna <- getPrimirna(toString(a))
  seq <- getSequence(toString(a))
  len <- 13
  duplimotif <- F
  motif <- get5PMotif(seq, len)
  while(!isUniqueMotif(motif, a, c()){
    len <- len + 1
    if(len > nchar(seq)){
      duplimotif <- T
      break
    }
    motif <- get5PMotif(seq, len)
  }

  
  unique_motifs <- rbind(unique_motifs, data.frame(MOTIF=motif, MIRNA=mirna, END=len, DUPLIMOTIF=duplimotif, stringsAsFactors = F))
}



for(i in 1:nrow(unique_motifs)){
  
  mirna <- toString(unique_motifs[i, 'MIRNA'])
  motif <- toString(unique_motifs[i, 'MOTIF'])
  
  
  contained_mirna <- c()
  for(x in 1:nrow(unique_motifs)){
    curr_mirna <- toString(unique_motifs[x, 'MIRNA'])
    if(gregexpr(getSequence(getAccession(mirna)), getSequence(getAccession(curr_mirna)), fixed=T)[[1]][1] == 1 & 
       gregexpr(mirna, curr_mirna, fixed=T)[[1]][1] == -1){
        curr_mirna_end <- strsplit(curr_mirna, '-', fixed=T)[[1]]
        curr_mirna_end <- paste(curr_mirna_end[3:length(curr_mirna_end)], collapse='-')
        contained_mirna <- c(contained_mirna, curr_mirna_end)
    }
  }
  
  unique_motifs[i, 'DUPLIMIRNA'] <- paste(contained_mirna, collapse=',')
  
}

unique_motifs <- read.table(paste0(base_path, 'miRBase21-unique-5p-motifs.tsv'), sep='\t', stringsAsFactors = F, header=T)
remove_mirna <- c()

for(i in 1:nrow(unique_motifs)){
  mirna <- unique_motifs[i, 'MIRNA']
  if(mirna%in%remove_mirna) { next }
  duplimirna <- unique_motifs[i, 'DUPLIMIRNA']
  if(length(duplimirna) < 1) { next }
  duplimirna <- strsplit(duplimirna, ',', fixed=T)[[1]]
  motif <- unique_motifs[i, 'MOTIF']
  final_duplimirna <- duplimirna
  for(d in duplimirna){
    curr_mirna <- paste0('hsa-miR-', d)
    curr_duplimirna <- unique_motifs[which(unique_motifs$MIRNA==curr_mirna), 'DUPLIMIRNA']
    if(length(curr_duplimirna) < 1){ next }
    curr_duplimirna <- strsplit(curr_duplimirna, ',', fixed=T)[[1]]

    if(gregexpr(getSequence(getAccession(mirna)), getSequence(getAccession(curr_mirna)), fixed=T)[[1]][1] == 1){
      remove_mirna <- c(curr_mirna, remove_mirna)
      for(c in curr_duplimirna){
        if(gregexpr(getSequence(getAccession(mirna)), getSequence(getAccession(paste0('hsa-miR-', c))), fixed=T)[[1]][1] == 1){
          final_duplimirna <- c(final_duplimirna, c)
          if(gregexpr(c, mirna, fixed=T)[[1]][1] == -1){
            remove_mirna <- c(paste0('hsa-miR-', c), remove_mirna)
          }
        }
      }
    }
  }
  for(f in final_duplimirna){
    if(gregexpr(toString(f), toString(unique_motifs[i, 'MIRNA'])) != -1){
      final_duplimirna <- final_duplimirna[!final_duplimirna == f]
    }
  }
  unique_motifs[i, 'DUPLIMIRNA'] <- paste(unique(final_duplimirna), collapse=',')
  unique_motifs[i, 'DUPLI.COUNT']<- length(final_duplimirna)
}

unique_motifs<-unique_motifs[!(unique_motifs$MIRNA%in%remove_mirna),]

unique_motifs$UNIQUE <- TRUE

for(i in 1:nrow(unique_motifs)){
  mirna <- unique_motifs[i, 'MIRNA']
  motif <- unique_motifs[i, 'MOTIF']
  duplimirna <- unique_motifs[i, 'DUPLIMIRNA']
  duplimirna <- strsplit(duplimirna, ',', fixed=T)[[1]]
  full_duplimirna <- c()
  for(d in duplimirna){
    full_duplimirna <- c(full_duplimirna, paste0('hsa-miR-', d))
  }
  len <- 13
  curr_motif <- get5PMotif(getSequence(getAccession(mirna)), len)
  while(!isUniqueMotif(curr_motif, getAccession(mirna), as.vector(full_duplimirna))){  
      len <- len+1
      curr_motif <- get5PMotif(getSequence(getAccession(mirna)), len)
      if(len >= nchar(getSequence(getAccession(mirna)))){
        unique_motifs[i, 'UNIQUE'] <- FALSE
        break
      }
  }
  unique_motifs[i, 'MOTIF'] <- curr_motif
  unique_motifs[i, 'END'] <- len
}

unique_motifs$DUPLI.COUNT<-NULL

for(i in 1:nrow(unique_motifs)){
  mirna <- unique_motifs[i, 'MIRNA']
  motif <- unique_motifs[i, 'MOTIF']
  duplimirna <- strsplit(unique_motifs[i, 'DUPLIMIRNA'], ',', fixed=T)[[1]]
  edited_duplimirna <- c()
  end <- unique_motifs[i, 'END']
  
  edited_mirna <- mirna
  for(d in duplimirna){
    edited_mirna <- paste0(edited_mirna, '-', d)
    edited_duplimirna <- c(edited_duplimirna, paste0('hsa-miR-', d))
  }
  edited_mirna <- paste0(edited_mirna, '-k', end)
  
  unique_motifs[i, 'EDITED.MIRNA'] <- edited_mirna
  
  unique_motifs[i, 'IS.UNIQUE'] <- F
  if(isUniqueMotif(motif, getAccession(mirna), edited_duplimirna)){
    unique_motifs[i, 'IS.UNIQUE'] <- T
  }
}


fa_data<-c()

for(i in 1:nrow(unique_motifs)){
  mirna <- unique_motifs[i, 'MIRNA']
  edited_mirna <- unique_motifs[i, 'EDITED.MIRNA']
  motif <- unique_motifs[i, 'MOTIF']
  sequence <- getSequence(getAccession(mirna))
  fa_data<-c(fa_data, paste(paste0('>', edited_mirna), toString(motif), sep='\t'), sequence)
}

motif_file<-file(paste0(base_path, 'miRBase21-unique-5p-motifs-collapsed-dupli.fa'))
writeLines(fa_data, motif_file)
close(motif_file)

write.table(unique_motifs, paste0(base_path, 'miRBase21-unique-5p-motifs.tsv'), sep='\t', row.names=F)
