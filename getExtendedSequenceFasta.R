
## PURPOSE: create fasta file containing extended sequences of specified miRNA
## INPUT: miRBase data containing UCSC extended sequences  miRBase21-master.tsv
## OUTPUT: fasta file of extended sequences   

require('readtext')

base_path<-'/Users/chens22/Documents/miRNA/structure/'

mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=T)

low_fidelity_mirna<-c('hsa-mir-105-1','hsa-mir-105-2', 'hsa-mir-101-1', 'hsa-mir-101-2', 
                      'hsa-mir-214', 'hsa-mir-30e', 'hsa-mir-455', 'hsa-mir-148a', 'hsa-mir-516b-1', 
                      'hsa-mir-516b-2', 'hsa-mir-526b', 'hsa-mir-10b', 'hsa-mir-29b-2', 'hsa-mir-217', 
                      'hsa-mir-200b', 'hsa-mir-136', 'hsa-mir-216b', 'hsa-mir-29a', 'hsa-mir-758', 
                      'hsa-mir-1266', 'hsa-mir-511', 'hsa-mir-24-2', 'hsa-mir-483', 'hsa-mir-34c', 
                      'hsa-mir-127', 'hsa-mir-370', 'hsa-mir-200a', 'hsa-mir-378a', 'hsa-mir-144', 
                      'hsa-mir-191', 'hsa-mir-532', 'hsa-mir-138-1', 'hsa-mir-138-2', 'hsa-let-7i')
high_fid_mirna<-c('hsa-mir-374a', 'hsa-mir-195', 'hsa-mir-96', 'hsa-mir-182', 'hsa-mir-155', 'hsa-mir-450a-1', 'hsa-mir-450a-2','hsa-mir-493', 'hsa-mir-154', 'hsa-mir-194-1', 'hsa-mir-194-2','hsa-mir-190a', 'hsa-mir-1271', 'hsa-mir-891a', 'hsa-mir-15b', 'hsa-mir-33a', 'hsa-mir-503', 'hsa-mir-885', 'hsa-mir-20a', 'hsa-mir-199a-1', 'hsa-mir-199a-2','hsa-mir-15a')

getPriMirna<-function(mirna){
  return(unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA']))
}

get5PSequence<-function(primirna){
  seq_5p<-unique(mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='5P'), 'SEQUENCE'])
  if(length(seq_5p)<1){
    print(paste0('NO 5P SEQUENCE: ', primirna))
    return('')
  }
  return(toString(seq_5p))
}

get3PSequence<-function(primirna){
  seq_3p<-unique(mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='3P'), 'SEQUENCE'])
  if(length(seq_3p)<1){
    print(paste0('NO 3P SEQUENCE: ', primirna))
    return('')
  }
  return(toString(seq_3p))
}

getSecondaryStructure<-function(sequence, mirna){
  input_data<-c(paste('>', mirna, sep=''), sequence)
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
  pri_seqs<-unique(mirna_data[which(mirna_data$PRIMIRNA==toString(getPriMirna(mirna))), 'PRI.SEQUENCE'])
  extended_seq<-getExtendedSequence(mirna)
  for(p in pri_seqs){
    if(gregexpr(toString(p), toString(extended_seq))[[1]][1]!=-1){
      return(toString(p))
    }
  }
}

getPrecursorSequence<-function(primirna){
  primirna<-toString(getPriMirna(primirna))
  pri_seq<-getPrimarySequence(primirna)
  pri_secondary_struct<-getSecondaryStructure(pri_seq, primirna)
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
  
  pre_seq<-substr(pri_seq, start, end)
  return(pre_seq)
}

getExtendedSequence<-function(primirna){
  seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(primirna)), 'EXTENDED.SEQUENCE'])
  return(toString(seq))
}

fa_data<-c()

for(m in high_fid_mirna){

  extended_seq<-getExtendedSequence(toString(m))
  pre_seq<-getPrecursorSequence(toString(factor(m)))

  pre_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  pre_end<-pre_start+nchar(pre_seq)-1

  pri_seq<-substr(extended_seq, pre_start-30, pre_end+30)
  struct<-getSecondaryStructure(pri_seq, '')
  fa_data<-c(fa_data, paste0('>', m), toString(pri_seq), struct)
}

seq_file<-file(paste0(base_path, 'high-fid-mirna-seq.txt'))
writeLines(fa_data, seq_file)
close(seq_file)

getComplementaryBase<-function(base){
  if(base%in%'A'){
    return('t')
  }
  if(base%in%'G'){
    return('c')
  }
  if(base%in%'C'){
    return('g')
  }
  if(base%in%'T'){
    return('a')
  }
}
