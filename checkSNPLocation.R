
## PURPOSE: locate SNP and determine if structural changes occur when replaced with alternative SNPs
## INPUT: miRNA SNP data from UCSC  mirna_hg38.tsv
##        miRBase data        miRBase21-master.tsv
## OUTPUT: fasta file containing location of structural change if any   UCSC-SNP-2d-struct-change.txt
##         table containing 2d struct of w/ alternate SNPs                   UCSC-SNPs-2d-structs.tsv

library(data.table)
require('stringr')
library(stringdist)

strReverse <- function(x){
  return(sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
}

strReverseComplement <- function(str){
  str <- chartr("ATGC","TACG",str)
  str <- sapply(lapply(strsplit(str, NULL), rev), paste, collapse="")
  return(str)
}

getExtendedSequence<-function(mirna){
  extended_seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(mirna)), 'EXTENDED.SEQUENCE'])
  if(length(extended_seq)>1){
    print(paste('ERROR EXTENDED SEQ', mirna))
    return(NA)
  }
  return(c(toString(extended_seq)))
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
  input_data<-c('>', sequence)
  input_file<-paste(base_path, 'temp-miRNA-seq.fa', sep='')
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-paste(base_path, 'temp-miRNA-secondary-struct.txt', sep='')
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
    start<-gregexpr(seq_5p, pri_seq)[[1]]
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
  
  pre_seq <- substr(pri_seq, start, end)
  
  if(is.na(pre_seq)){
    return('')
  }
  return(pre_seq)
}


base_path <-'/Users/chens22/Documents/'

mirna_snp_file <- paste0(base_path, 'genome/miRNA_hg38.tsv')
mirna_snp_data <- read.table(mirna_snp_file, sep='\t', stringsAsFactors=F, header=T)
mirna_snp_data <- mirna_snp_data[which(mirna_snp_data$type=='miRNA'),]
mirna_snp_data <- mirna_snp_data[which(nchar(mirna_snp_data$snp)>0),]

mirna_file <- paste0(base_path, 'miRBase21/miRBase21-master.tsv')
mirna_data <- read.table(mirna_file, sep='\t', header=T, stringsAsFactors=F)
mirna_data <- mirna_data[which(mirna_data$PRIMIRNA%in%mirna_snp_data$name),]

#mirna_data <- mirna_data[,c('PRI.ACCESSION', 'PRIMIRNA','ACCESSION', 'STRAND','EXTENDED.SEQUENCE', 'X.COORDINATE', 'Y.COORDINATE', 'SEQUENCE', 'CHROMOSOME')]


#reverse complement of extended sequences
# for(i in 1:nrow(mirna_data)){
# 
#   extended_seq <- toString(mirna_data[i, 'EXTENDED.SEQUENCE'])
#   direction <- toString(mirna_data[i, 'DIRECTION'])
#   if(!direction%in%'-'){ next }
# 
#   print(extended_seq)
#   extended_seq <- chartr("ATGC","TACG",extended_seq)
#   extended_seq <- strReverse(c(extended_seq))[1]
#   print(extended_seq)
#   mirna_data[i, 'EXTENDED.SEQUENCE'] <- extended_seq
# 
# }


precursor_data <- data.frame(PRIMIRNA = unique(mirna_data$PRIMIRNA), PRI.ACCESSION = unique(mirna_data$PRI.ACCESSION))
pri_accession <- unique(mirna_data$PRI.ACCESSION)
for(i in 1:nrow(precursor_data)){
  
  a <- precursor_data[i, 'PRI.ACCESSION']
  primirna <- precursor_data[i, 'PRIMIRNA']
  #if(primirna%in%except){ next }
  extended_seq <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==a), 'EXTENDED.SEQUENCE'])
  seq_5p <- mirna_data[which(mirna_data$PRI.ACCESSION==a & mirna_data$STRAND=='5P'), 'SEQUENCE']
  seq_3p <- mirna_data[which(mirna_data$PRI.ACCESSION==a & mirna_data$STRAND=='3P'), 'SEQUENCE']
  direction <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==a), 'DIRECTION'])
  chromosome <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==a), 'CHROMOSOME'])
  both_mature <- T
  
  if(length(seq_5p) < 1 || length(seq_3p) < 1){ 
    both_mature <- F
    pre_seq <- getPrecursorSequence(primirna)
    start_5p <- gregexpr(pre_seq, extended_seq, fixed=T)[[1]][1]
    end_3p <- start_5p + nchar(pre_seq) - 1
  } 
  
  if(direction%in%'-'){
    extended_seq <- strReverseComplement(extended_seq)
    seq_5p <- chartr("ATGC","TACG",seq_5p)
    seq_5p <- strReverse(c(seq_5p))[1]
    seq_3p <- chartr("ATGC","TACG",seq_3p)
    seq_3p <- strReverse(c(seq_3p))[1]
    temp <- seq_5p
    seq_5p <- seq_3p
    seq_3p <- temp
    
  }

  
  
  start <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==a), 'X.COORDINATE'])
  end <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==a), 'Y.COORDINATE'])
  
  if(nchar(toString(seq_5p)) > 0 & nchar(toString(seq_3p)) > 0 & toString(seq_5p) != 'NULL' & toString(seq_3p) != 'NULL'){
    start_5p <- gregexpr(toString(seq_5p), toString(extended_seq), fixed=T)[[1]][1]
    start_3p <-gregexpr(toString(seq_3p), toString(extended_seq), fixed=T)[[1]][1]
    end_3p <- start_3p + nchar(seq_3p) - 1
  }

  
  pre_seq <- substr(extended_seq, start_5p, end_3p)
  pre_start <- start_5p + start - 1
  pre_end <- end_3p + start - 1
  
  pri_start <- start_5p - 13
  pri_end <- end_3p + 13
  pri_seq <- substr(extended_seq, pri_start, pri_end)
  pri_start <- pri_start + start
  pri_end <- pri_end + start
  
  if(nchar(pri_seq)<1){ next }
  precursor_data[i, 'PRE.SEQUENCE'] <- pre_seq
  precursor_data[i, 'PRE.START'] <- pre_start
  precursor_data[i, 'PRE.END'] <- pre_end
  
  precursor_data[i, 'PRI.SEQUENCE'] <- pri_seq
  precursor_data[i, 'PRI.START'] <- pri_start
  precursor_data[i, 'PRI.END'] <- pri_end
  
  precursor_data[i, 'CHROMOSOME'] <- chromosome
  precursor_data[i, 'DIRECTION'] <- direction
  precursor_data[i, 'BOTH.MATURE'] <- both_mature
  
}

precursor_data <- precursor_data[complete.cases(precursor_data),]


snp <- c()

for(s in mirna_snp_data$snp){
  if(nchar(s)>0){
    curr_snp <- strsplit(toString(s), ',', fixed=T)[[1]]
    snp <- c(snp, curr_snp)
  }
}


snp_file <- paste0(base_path, 'genome/SNP_hg38_Common_SNPs150')
#snp_data <- fread(snp_file, sep='\t', header=T, stringsAsFactors=F)
snp_data <- snp_data[which(snp_data$name%in%snp),]

data <- data.frame(SNP.NAME=character(), stringsAsFactors = F)
for(i in 1:nrow(mirna_snp_data)){
  
  mirna <- mirna_snp_data[i, 'name']
  snp <- mirna_snp_data[i, 'snp']
  snp <- strsplit(snp, ',', fixed=T)[[1]]
  
  pri_seq <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'PRI.SEQUENCE']
  pri_start <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'PRI.START']
  pri_end <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'PRI.END']
  chromosome <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'CHROMOSOME']
  direction <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'DIRECTION']
  accession <- unique(precursor_data[which(precursor_data$PRIMIRNA==mirna),'PRI.ACCESSION'])
  
  for(s in snp){
    
    snp_chromosome <- snp_data[which(snp_data$name==s), 'chrom']
    if(toString(chromosome)!=toString(snp_chromosome)){
      print(s)
      next
    }
    
    snp_start <- snp_data[which(snp_data$name==s & snp_data$chrom==chromosome), 'chromStart']
    if(is.na(as.integer(snp_start))){
      snp_start <- snp_data[which(snp_data$name==s), 'chromStart']
    }
    snp_end <- snp_data[which(snp_data$name==s & snp_data$chrom==chromosome), 'chromEnd']
    if(is.na(as.integer(snp_end))){
      snp_end <- snp_data[which(snp_data$name==s), 'chromEnd']
    }
    
    snp_direction <- snp_data[which(snp_data$name==s & snp_data$chrom==chromosome), 'strand']
    
    if(s=="rs2368393"){
     # snp_direction <- '+'
    }
    
    if(length(pri_seq) < 1){ next }
    
    observed <- snp_start >= pri_start && snp_start <= pri_end && snp_end <= pri_end
    
    observed_snp <- ''
    if(observed==T){
      start <- as.integer(snp_start - pri_start) + 3#R starts at 1
      if(start < 1){ start <- 1 }
      end <- as.integer(start + (snp_end-snp_start-1))
      if(snp_direction%in%'-'&direction%in%'-'){
        pri_seq <- strReverseComplement(pri_seq)
        start <- as.integer((nchar(pri_seq)) - start - (snp_end-snp_start-1))
        if(start < 1){ start <- 1 }
        end <- as.integer(start + (snp_end-snp_start-1))
      }
      if(s%in%c('rs2368393', 'rs11020790')){
        pri_seq <- strReverseComplement(pri_seq)
      }
      if(snp_direction%in%'-' & direction%in%'+'){
        start <- as.integer(nchar(pri_seq) - start - (snp_end-snp_start-1)) + 1
        if(start < 1){ start <- 1 }
        end <- as.integer(start + (snp_end-snp_start-1))
        pri_seq <- strReverseComplement(pri_seq)
      }
      observed_snp <- substr(toString(pri_seq), start, end)
    }
    
    if(snp_end==pri_end & snp_start >= snp_end - 1 & (snp_end-snp_start)>0){
      observed_snp <- substr(pri_seq, nchar(pri_seq), nchar(pri_seq))
    }
      
    row <- nrow(data)+1
    data[row, 'SNP.NAME'] <- toString(s)
    data[row, 'SNP'] <- snp_data[which(snp_data$name==s & snp_data$chrom==chromosome), 'observed']
    data[row, 'OBSERVED'] <- observed
    data[row, 'OBSERVED.SNP'] <- toString(observed_snp)
    data[row, 'PRI.DIRECTION'] <- direction
    data[row, 'SNP.DIRECTION'] <- snp_direction
    data[row, 'PRIMIRNA'] <- mirna
    data[row, 'PRI.ACCESSION'] <- accession
    data[row, 'SNP.START'] <- snp_start
    data[row, 'SNP.END'] <- snp_end
    data[row, 'CHROMOSOME'] <- chromosome
    data[row, 'PRI.SEQUENCE'] <- pri_seq
    data[row, 'PRI.START'] <- pri_start
    data[row, 'PRI.END'] <- pri_end
    data[row, 'START'] <- start
    data[row, 'END'] <- end
    data[row, 'BOTH.MATURE'] <- precursor_data[which(precursor_data$PRIMIRNA==mirna),'BOTH.MATURE']
    
    
  }
  
}

data <- data[which(data$OBSERVED==T),]

unmatched1 <- c()
for(i in 1:nrow(data)){
  
  observed_snp <- data[i, 'OBSERVED.SNP']
  snp <- data[i, 'SNP']
  
  if(gregexpr(observed_snp, snp, fixed=T)[[1]][1] == -1 || (nchar(observed_snp) < 1 && gregexpr('-', snp, fixed=T)[[1]][1] == -1)){
    unmatched1 <- c(data[i, 'SNP.NAME'], unmatched1)
  }
}

unmatched <- c()
for(i in 1:nrow(snp_data)){
  
  observed_snp <- snp_data[i, 'refUCSC']
  snp <- snp_data[i, 'observed']
  
  if(gregexpr(observed_snp, snp, fixed=T)[[1]][1] == -1 || (nchar(observed_snp) < 1 && gregexpr('-', snp, fixed=T)[[1]][1] == -1)){
    unmatched <- c(toString(snp_data[i, 'name']), unmatched)
  }
}

unmatched <- intersect(unmatched, data$SNP.NAME)

struct_data <- data.frame(SNP.NAME=character(), SNP=character(), SECONDARY.STRUCTURE=character(), FOUND.UCSC=logical(), stringsAsFactors=F)
for(i in 1:nrow(data)){

  snp_name <- data[i, 'SNP.NAME']
  if(snp_name%in%c('rs781774720', 'rs13447640', 'rs2368392', 'rs2291418', 'rs2241347', 'rs2747232')){next}
  
  snp_start <- data[i, 'START']
  snp_end <- data[i, 'END']
  pri_seq <- data[i, 'PRI.SEQUENCE']
  primirna <- data[i, 'PRIMIRNA']
  observed_snp <- c(data[i, 'OBSERVED.SNP'])
  
  snp <- data[i, 'SNP']
  snp <- strsplit(toString(snp), '/', fixed=T)[[1]]
  snp <- replace(snp, snp=='-', '')
  
  for(s in snp){
    
    curr_seq <- paste0(substr(pri_seq, 1, snp_start-1), s, substr(pri_seq, snp_end + 1, nchar(pri_seq)))
    struct <- getSecondaryStructure(curr_seq)
    
    row <- nrow(struct_data)+1
    struct_data[row, 'SNP.NAME'] <- snp_name
    struct_data[row, 'SNP'] <- toString(s)
    struct_data[row, 'SECONDARY.STRUCTURE'] <- struct
    struct_data[row, 'FOUND.UCSC'] <- s %in% observed_snp
    struct_data[row, 'PRIMIRNA'] <- primirna
    
  }
  
}





struct_data <- unique(struct_data)
for(i in 1:nrow(struct_data)){
  
  snp <- struct_data[i, 'SNP.NAME']

  primirna <- struct_data[i, 'PRIMIRNA']
  observed_struct <- unique(struct_data[which(struct_data$SNP.NAME==snp & struct_data$FOUND.UCSC==T & struct_data$PRIMIRNA==primirna), 'SECONDARY.STRUCTURE'])
  if(nchar(toString(observed_struct)) < 1){
    observed_snp <- toString(snp_data[which(snp_data$name==snp_name), 'refUCSC'])
    observed_struct <- struct_data[which(struct_data$SNP.NAME==snp & struct_data$SNP==observed_snp & struct_data$PRIMIRNA==primirna), 'SECONDARY.STRUCTURE']
    struct_data[which(struct_data$SNP.NAME==snp & struct_data$SNP==observed_snp & struct_data$PRIMIRNA==primirna), 'FOUND.UCSC'] <- T
    
  }
  struct <- struct_data[i, 'SECONDARY.STRUCTURE']
  struct <- gsub(')', 'c', gsub('(', 'b', gsub('.', 'a', struct, fixed=T), fixed=T), fixed=T)
  observed_struct <- gsub(')', 'c', gsub('(', 'b', gsub('.', 'a', observed_struct, fixed=T), fixed=T), fixed=T)
  struct_data[i, 'DISTANCE.FROM.UCSC'] <- stringdist(struct, toString(observed_struct))
  
}
struct_data <- merge(struct_data, data[, c('SNP.NAME', 'PRIMIRNA', 'PRI.ACCESSION','PRI.SEQUENCE', 'PRI.START', 'PRI.END', 'PRI.DIRECTION')], by=c('SNP.NAME', 'PRIMIRNA'), all=F)

write.table(struct_data, file=paste0(base_path, '/genome/UCSC-SNPs-2d-structs.tsv'), sep='\t', row.names=F)

fasta_data <- c()
for(i in 1:nrow(struct_data)){
  accession <-  struct_data[i, 'PRI.ACCESSION']
  snp_name <- struct_data[i, 'SNP.NAME']
  snp <- struct_data[i, 'SNP']
  
  pri_start <- struct_data[i, 'PRI.START']
  pri_end <- struct_data[i, 'PRI.END']
  
  pre_start <- unique(precursor_data[which(precursor_data$PRI.ACCESSION==accession), 'PRE.START'])
  pre_end <- unique(precursor_data[which(precursor_data$PRI.ACCESSION==accession), 'PRE.END'])
  
  seq_5p <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==accession & mirna_data$STRAND=='5P'), 'SEQUENCE'])
  seq_3p <- unique(mirna_data[which(mirna_data$PRI.ACCESSION==accession & mirna_data$STRAND=='3P'), 'SEQUENCE'])
  
  if(nchar(toString(seq_5p))==0){
    seq_5p <- seq_3p
  }
  if(nchar(toString(seq_3p))==0){
    seq_3p <- seq_5p
  }
  
  start <- pre_start - pri_start + 2
  end <- start + (pre_end - pre_start)
  end_5p <- start+nchar(seq_5p)-1
  start_3p <- end - nchar(seq_3p) + 1
  
  snp_start <- data[which(data$PRI.ACCESSION==accession & data$SNP.NAME==snp_name), 'START']
  observed_snp <- unique(struct_data[which(struct_data$SNP.NAME==snp_name & struct_data$FOUND.UCSC==T & struct_data$PRI.ACCESSION==accession), 'SNP'])
  
  snp_len_diff <- nchar(observed_snp) - nchar(snp)
  
  positions <- c(start, end_5p, start_3p, end)
  positions <- positions[order(positions)]
  
  snp_location <- positions[which.min(abs(positions-snp_start))]

  while(snp_start > snp_location & snp_location < positions[4]){
    snp_location <- positions[match(snp_location, positions)+1]
  }
  
  if(snp_location!=positions[4]){
    positions[match(snp_location, positions)] <- positions[match(snp_location, positions)] + snp_len_diff
  }
  
  
  struct <- struct_data[i, 'SECONDARY.STRUCTURE']
  observed_struct <- unique(struct_data[which(struct_data$SNP.NAME==snp_name & struct_data$FOUND.UCSC==T & struct_data$PRI.ACCESSION==accession), 'SECONDARY.STRUCTURE'])
  
  struct <- c(substr(struct, 1, positions[1]-1), substr(struct, positions[1], end_5p), substr(struct, positions[2]+1, positions[3]-1),  substr(struct, positions[3], end), substr(struct, positions[4], nchar(struct)))
  observed_struct <- c(substr(observed_struct, 1, start-1), substr(observed_struct, start, start+nchar(seq_5p)-1), substr(observed_struct, start+nchar(seq_5p), end - nchar(seq_3p) + 1),  substr(observed_struct, end - nchar(seq_3p), end), substr(observed_struct, end, nchar(observed_struct)))

  parts <- c('5p-lowerstem', '5p-arm', 'loop', '3p-arm', '3p-lowerstem')
  
  change <- c()
  for(x in 1:length(struct)){
    curr_struct <- struct[x]
    curr_observed <- observed_struct[x]
    
    dist <- stringdist(curr_struct, curr_observed)
    if(dist > 0){
      change <- c(parts[x], change)
    }
  }
  if(length(change)>0){
    #change <- 'NONE'
    fasta_data <- c(fasta_data, paste(accession, snp_name, ', '))
    fasta_data <- c(fasta_data, change)
    fasta_data <- c(fasta_data, paste(observed_struct, collapse=' '))
    fasta_data <- c(fasta_data, paste(struct, collapse=' '), '\n')
  }
  
  struct_data[i, 'CHANGE'] <- paste(change, collapse=', ')
  
}

file<-file(paste0(base_path, '/genome/UCSC-SNP-2d-struct-change.txt'))
writeLines(fasta_data, file)
close(file)
write.table(struct_data, file=paste0(base_path, '/genome/UCSC-SNPs-2d-structs.tsv'), sep='\t', row.names=F)

#add line under every change 5p-arm, 5p-lowerstem, loop



