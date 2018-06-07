## PURPOSE: return moR base pairing data for given miRNA
## INPUT: miRBase data  miRBase21-master.tsv
##        extended 2d structures from getConservedExtended2dStruct.R  extended-2d-struct-max-moR-pairing-all-miRBase.tsv
## OUTPUT: moR pairing                  miRNA-5p-moR-pairing-by-fidelity.tsv
##         primiRNA 2d struct data  primary-extended-secondary-structure.tsv

require('stringr')
require('plyr')

base_path<-'/Users/chens22/Documents/structure/'
mirna_file<-'/Users/chens22/Documents/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors=F)
conserved_struct_file<-'/Users/chens22/Documents/structure/extended-2d-struct-max-moR-pairing-all-miRBase.tsv'
conserved_struct_data<-read.table(file=conserved_struct_file, sep='\t', header=T, stringsAsFactors=F)

struct_len<-30


low_fid_mirna<-c('hsa-miR-758-5p', 'hsa-miR-153-5p','hsa-miR-383-5p','hsa-miR-382-5p','hsa-miR-183-5p',
                  'hsa-miR-425-5p','hsa-miR-136-5p','hsa-miR-628-5p','hsa-miR-769-5p','hsa-miR-30e-5p',
                  'hsa-miR-98-5p','hsa-miR-138-5p-1-2','hsa-miR-144-5p','hsa-let-7i-5p', 'hsa-miR-191-5p',
                 'hsa-miR-532-5p')
low_fid_mirna<-c('hsa-mir-3195', 'hsa-mir-4516', 'hsa-mir-7704', 'hsa-mir-29c', 'hsa-mir-452', 
                  'hsa-mir-4286', 'hsa-mir-142', 'hsa-mir-24-1', 'hsa-mir-342', 'hsa-mir-411', 
                  'hsa-mir-183', 'hsa-mir-942', 'hsa-mir-382', 'hsa-mir-654', 'hsa-mir-425', 
                'hsa-mir-299', 'hsa-mir-582', 'hsa-mir-10a', 'hsa-mir-4791', 'hsa-mir-577', 
                'hsa-mir-27b', 'hsa-mir-383', 'hsa-mir-628', 'hsa-mir-378g', 'hsa-mir-769')
low_fid_mirna<-c('hsa-mir-3195', 'hsa-mir-4516', 'hsa-mir-7704', 'hsa-mir-29c', 'hsa-mir-452', 'hsa-mir-4286', 'hsa-mir-142', 'hsa-mir-24-1', 'hsa-mir-342', 'hsa-mir-411', 'hsa-mir-183', 'hsa-mir-942', 'hsa-mir-382', 'hsa-mir-654', 'hsa-mir-425', 'hsa-mir-299', 'hsa-mir-582', 'hsa-mir-10a', 'hsa-mir-4791', 'hsa-mir-577', 'hsa-mir-27b', 'hsa-mir-383', 'hsa-mir-628', 'hsa-mir-378g', 'hsa-mir-769', 'hsa-mir-105-1', 'hsa-mir-105-2', 'hsa-mir-101-1', 'hsa-mir-101-2', 'hsa-mir-214', 'hsa-mir-30e', 'hsa-mir-455', 'hsa-mir-148a', 'hsa-mir-516b-1', 'hsa-mir-516b-2', 'hsa-mir-526b', 'hsa-mir-10b', 'hsa-mir-29b-2', 'hsa-mir-217', 'hsa-mir-200b', 'hsa-mir-136', 'hsa-mir-216b', 'hsa-mir-29a', 'hsa-mir-758', 'hsa-mir-1266', 'hsa-mir-511', 'hsa-mir-24-2', 'hsa-mir-483', 'hsa-mir-34c', 'hsa-mir-127', 'hsa-mir-370', 'hsa-mir-200a', 'hsa-mir-378a')
mid_fid_mirna<-c('hsa-miR-139-5p','hsa-miR-432-5p', 'hsa-miR-9-5p-1-2-3','hsa-let-7g-5p',
                  'hsa-miR-218-5p-1-2','hsa-miR-186-5p','hsa-miR-93-5p','hsa-miR-99a-5p','hsa-let-7d-5p',
                  'hsa-miR-16-5p-1-2','hsa-miR-744-5p','hsa-miR-182-5p','hsa-miR-584-5p',
                  'hsa-miR-181a-5p-1-2','hsa-miR-20a-5p','hsa-miR-885-5p','hsa-miR-15a-5p','hsa-miR-409-5p')
high_fid_mirna<-c('hsa-miR-100-5p','hsa-miR-106b-5p','hsa-miR-129-5p-1-2','hsa-miR-589-5p','hsa-miR-125b-5p-1-2',
                  'hsa-miR-30a-5p','hsa-miR-125a-5p','hsa-miR-134-5p','hsa-miR-99b-5p','hsa-miR-330-5p',
                  'hsa-miR-26a-5p-1-2','hsa-miR-21-5p','hsa-let-7b-5p','hsa-miR-1307-5p','hsa-miR-185-5p',
                  'hsa-miR-361-5p','hsa-miR-379-5p','hsa-let-7f-5p-1-2','hsa-let-7e-5p','hsa-miR-340-5p',
                  'hsa-miR-149-5p','hsa-miR-204-5p','hsa-miR-26b-5p','hsa-miR-30d-5p')
high_fid_mirna<-c('hsa-mir-508', 'hsa-mir-30d', 'hsa-mir-26b', 'hsa-mir-193a', 'hsa-let-7e', 
                  'hsa-let-7b', 'hsa-let-7f-1', 'hsa-let-7f-2', 'hsa-mir-204', 'hsa-mir-450b', 
                  'hsa-mir-361', 'hsa-mir-1287', 'hsa-mir-4662a', 'hsa-mir-379', 'hsa-mir-26a-1', 
                  'hsa-mir-26a-2', 'hsa-mir-134', 'hsa-mir-21', 'hsa-mir-185', 'hsa-mir-211', 
                  'hsa-mir-660', 'hsa-mir-149', 'hsa-mir-1307', 'hsa-mir-340', 'hsa-mir-4774')
high_fid_mirna<-c('hsa-mir-508', 'hsa-mir-30d', 'hsa-mir-26b', 'hsa-mir-193a', 'hsa-let-7e', 'hsa-let-7b', 'hsa-let-7f-1', 'hsa-let-7f-2', 'hsa-mir-204', 'hsa-mir-450b', 'hsa-mir-361', 'hsa-mir-1287', 'hsa-mir-4662a', 'hsa-mir-379', 'hsa-mir-26a-1', 'hsa-mir-26a-2', 'hsa-mir-134', 'hsa-mir-21', 'hsa-mir-185', 'hsa-mir-211', 'hsa-mir-660', 'hsa-mir-149', 'hsa-mir-1307', 'hsa-mir-340', 'hsa-mir-4774', 'hsa-mir-125a', 'hsa-mir-374a', 'hsa-mir-195', 'hsa-mir-96', 'hsa-mir-153-2', 'hsa-mir-99b', 'hsa-mir-125b-1', 'hsa-mir-125b-2', 'hsa-mir-182', 'hsa-mir-155', 'hsa-mir-450a-1', 'hsa-mir-450a-2', 'hsa-mir-493', 'hsa-mir-154', 'hsa-mir-194-1', 'hsa-mir-194-2', 'hsa-mir-330', 'hsa-mir-190a', 'hsa-mir-1271', 'hsa-mir-891a', 'hsa-mir-15b', 'hsa-mir-30a', 'hsa-mir-33a', 'hsa-mir-503', 'hsa-mir-885', 'hsa-mir-106b', 'hsa-mir-20a', 'hsa-mir-199a-1', 'hsa-mir-199a-2', 'hsa-mir-589', 'hsa-mir-15a')

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

getMaximumBracketsExtendedSequence<-function(primirna){
  primirna<-getPriMirna(primirna)
  extended_seq<-getExtendedSequence(primirna)
  if(!primirna%in%conserved_struct_data$PRIMIRNA){
    return('')
  }
  start<-as.numeric(conserved_struct_data[which(conserved_struct_data$PRIMIRNA==primirna),'START']*-1)
  end<-as.numeric(conserved_struct_data[which(conserved_struct_data$PRIMIRNA==primirna),'END'])
  pre_seq<-getPrecursorSequence(primirna)
  if(is.na(pre_seq)){
    print(paste('cannot retrieve precursor seq: ', primirna))
    return('')
  }
  pre_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  pre_end<-pre_start+nchar(pre_seq)-1
  start<-pre_start-start
  end<-pre_end+end
  max_extended_seq<-substr(extended_seq, start, end)
  return(toString(max_extended_seq))
}

getSequence<-function(mirna){
  return(toString(unique(mirna_data[which(mirna_data$MIRNA==mirna), 'SEQUENCE'])))
}

get5PMorSecondaryStructure<-function(primirna, extended_seq){ #get secondary structure of mature miRNA
  primirna<-getPriMirna(primirna)
  mirna<-get5PMirna(primirna)
  seq<-getSequence(mirna)
  pri_seq<-getPrimarySequence(primirna)
  for(p in pri_seq){
    seq_start<-gregexpr(toString(seq), toString(p), fixed=T)[[1]]
    pri_start<-gregexpr(toString(p), extended_seq, fixed=T)[[1]]
    if(length(seq_start)>1){
      print('mature sequence appears in primary sequence multiple times')
      print(mirna)
    }
    if(seq_start>0 & pri_start>0){
      pri_seq<-toString(p)
      break
    }
  }
  pri_seq_start<-gregexpr(pri_seq, extended_seq, fixed=T)[[1]]
  pre_seq<-getPrecursorSequence(primirna)
  pre_seq_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  extended_seq<-substr(extended_seq, pre_seq_start-struct_len, pre_seq_start+nchar(pre_seq)+struct_len-1)
  secondary_struct<-getSecondaryStructure(extended_seq, mirna)#getPriSecondaryStructure(mirna)
  if(nchar(pre_seq)<1){
    print(paste0('NO PRECURSOR:', primirna))
    return('')
  }
  seq_start<-gregexpr(seq, extended_seq, fixed=TRUE)[[1]]
  if(length(seq_start)>1){
    pri_seq_start<-gregexpr(pri_seq, extended_seq, fixed=T)[[1]]
    seq_start_pri<-gregexpr(seq, pri_seq, fixed=T)[[1]][1]
    seq_start<-pri_seq_start+seq_start_pri-1
  }
  bases_len<-0
  mor_secondary_struct<-c()
  temp_mor_secondary_struct<-c()
  i<-0
  #seq_start<-1
  while(seq_start-1-i>=1 & bases_len<struct_len){
    curr_base<-substr(secondary_struct, seq_start-1-i, seq_start-1-i)
    if(curr_base%in%'(' || curr_base%in%'.'){
      mor_secondary_struct<-c(curr_base, mor_secondary_struct)
      temp_mor_secondary_struct<-c(curr_base, temp_mor_secondary_struct)
      bases_len<-bases_len+1
    }else{
      paired_bases<-0
      while(curr_base%in%')'||  curr_base%in%'.'){
        curr_base<-substr(secondary_struct, seq_start-1-i, seq_start-1-i)
        if(curr_base%in%')'){
          paired_bases<-paired_bases+1
          mor_secondary_struct<-c('*', mor_secondary_struct)
        }else{
          mor_secondary_struct<-c('.', mor_secondary_struct)
        }
        temp_mor_secondary_struct<-c(curr_base, temp_mor_secondary_struct)
        i<-i+1
      }
      paired_bases<-paired_bases-1
      mor_secondary_struct<-c('*', mor_secondary_struct)
      temp_mor_secondary_struct<-c(curr_base, temp_mor_secondary_struct)
      while(paired_bases>1 & seq_start-1-i>0){
        if(curr_base%in%'('){
          paired_bases<-paired_bases-1
          mor_secondary_struct<-c('*', mor_secondary_struct)
        }else{
          mor_secondary_struct<-c('.', mor_secondary_struct)
        }
        temp_mor_secondary_struct<-c(curr_base, temp_mor_secondary_struct)
        i<-i+1
        curr_base<-substr(secondary_struct, seq_start-1-i, seq_start-1-i)
      }
      
    }
    i<-i+1
  }
  mor_secondary_struct<-mor_secondary_struct[(length(mor_secondary_struct)-(struct_len-1)):length(mor_secondary_struct)]
  temp_mor_secondary_struct<-temp_mor_secondary_struct[(length(temp_mor_secondary_struct)-(struct_len-1)):length(temp_mor_secondary_struct)]
  mor_seq<-substr(extended_seq, length(temp_mor_secondary_struct)-(struct_len-1), length(temp_mor_secondary_struct))
  return(c(paste(mor_secondary_struct, collapse=''), paste(temp_mor_secondary_struct, collapse=''), mor_seq, extended_seq, secondary_struct))
}

get3PMorSecondaryStructure<-function(primirna, extended_seq){ #get secondary structure of mature miRNA
  primirna<-getPriMirna(primirna)
  mirna<-get3PMirna(primirna)
  seq<-getSequence(mirna)
  pri_seq<-getPrimarySequence(primirna)
  for(p in pri_seq){
    seq_start<-gregexpr(toString(seq), toString(p), fixed=T)[[1]]
    pri_start<-gregexpr(p, extended_seq, fixed=T)[[1]]
    if(seq_start>0 & pri_start>0){
      pri_seq<-p
      break
    }
  }
  pri_seq_start<-gregexpr(pri_seq, extended_seq, fixed=T)[[1]]
  # extended_seq<-substr(extended_seq, pri_seq_start-struct_len, pri_seq_start+nchar(pri_seq)+struct_len-1)
  # secondary_struct<-getSecondaryStructure(extended_seq, mirna)
  pre_seq<-getPrecursorSequence(primirna)
  if(nchar(pre_seq)<1){
    print(paste0('NO PRECURSOR:', primirna))
    return('')
  }
  pre_seq_start<-gregexpr(pre_seq, extended_seq, fixed=TRUE)[[1]]
  extended_seq<-substr(extended_seq, pre_seq_start-struct_len, pre_seq_start+nchar(pre_seq)+struct_len-1)
  secondary_struct<-getSecondaryStructure(extended_seq, mirna)
  seq_start<-gregexpr(seq, extended_seq, fixed=TRUE)[[1]]
  if(length(seq_start)>1){
    pri_seq_start<-gregexpr(pri_seq, extended_seq, fixed=T)[[1]]
    seq_start_pri<-gregexpr(seq, pri_seq, fixed=T)[[1]]
    seq_start<-pri_seq_start+seq_start_pri
  }
  
  
  bases_len<-0
  mor_secondary_struct<-c()
  temp_mor_secondary_struct<-c()
  i<-0
  while(seq_start+nchar(seq)+i<=nchar(secondary_struct) & bases_len<struct_len){
    curr_base<-substr(secondary_struct, seq_start+nchar(seq)+i, seq_start+nchar(seq)+i)
    if(curr_base%in%')' || curr_base%in%'.'){
      mor_secondary_struct<-c(mor_secondary_struct, curr_base)
      temp_mor_secondary_struct<-c(temp_mor_secondary_struct, curr_base)
      bases_len<-bases_len+1
    }else{
      paired_bases<-0
      while(curr_base%in%'('||  curr_base%in%'.'){
        curr_base<-substr(secondary_struct, seq_start+nchar(seq)+i, seq_start+nchar(seq)+i)
        if(curr_base%in%'('){
          paired_bases<-paired_bases+1
          mor_secondary_struct<-c(mor_secondary_struct, '*')
        }else{
          mor_secondary_struct<-c(mor_secondary_struct, '.')
        }
        temp_mor_secondary_struct<-c(temp_mor_secondary_struct, curr_base)
        i<-i+1
      }
      paired_bases<-paired_bases-1
      mor_secondary_struct<-c(mor_secondary_struct, '*')
      temp_mor_secondary_struct<-c(temp_mor_secondary_struct, curr_base)
      while(paired_bases>1 & seq_start-1-i>0){
        if(curr_base%in%')'){
          paired_bases<-paired_bases-1
          mor_secondary_struct<-c(mor_secondary_struct, '*')
        }else{
          mor_secondary_struct<-c(mor_secondary_struct, '.')
        }
        
        temp_mor_secondary_struct<-c(temp_mor_secondary_struct, curr_base)
        i<-i+1
        curr_base<-substr(secondary_struct, seq_start+nchar(seq)+i, seq_start+nchar(seq)+i)
      }
    }
    i<-i+1
  }
  mor_secondary_struct<-mor_secondary_struct[1:struct_len]
  temp_mor_secondary_struct<-temp_mor_secondary_struct[1:struct_len]
  mor_seq<-substr(extended_seq, 1, struct_len)

  return(c(paste(mor_secondary_struct, collapse=''), paste(temp_mor_secondary_struct, collapse=''), mor_seq, extended_seq, secondary_struct))
}

get5PMorSecondaryStructureList<-function(primirna_list){
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, EXTENDED.SEQUENCE=NA, SECONDARY.STRUCTURE=NA, stringsAsFactors = F)
  for(p in primirna_list){
    m<-get5PMirna(p)
    print(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      extended_sequences<-getExtendedSequence(p)
      mor_struct<-get5PMorSecondaryStructure(p, toString(extended_sequences))
      if(!any(is.na(mor_struct))){
        curr_data[, 'MOR.EDITED']<-mor_struct[1]
        curr_data[, 'MOR']<-mor_struct[2]
        curr_data[, 'MOR.SEQUENCE']<-mor_struct[3]
        curr_data[, 'EXTENDED.SEQUENCE']<-mor_struct[4]
        curr_data[, 'SECONDARY.STRUCTURE']<-mor_struct[5]
        for(i in 1:nchar(mor_struct[1])){
          curr_data[, toString(-(nchar(mor_struct[1])-i+1))]<-substr(mor_struct[1], i, i)
        }
        data<-rbind.fill(data, curr_data)
        curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
      }
    }else{
      next
    }
  }
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

get3PMorSecondaryStructureList<-function(primirna_list){
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, EXTENDED.SEQUENCE=NA, SECONDARY.STRUCTURE=NA, stringsAsFactors = F)
  for(p in primirna_list){
    print(p)
    m<-get3PMirna(p)
    extended_sequences<-c()
    extended_sequences<-getExtendedSequence(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      extended_sequences<-getExtendedSequence(p)
      #edit extended sequence to sequence with maximum brackets in 2d structure
      mor_struct<-get3PMorSecondaryStructure(p, toString(extended_sequences))
      if(!any(is.na(mor_struct))){
        curr_data[, 'MOR.EDITED']<-mor_struct[1]
        curr_data[, 'MOR']<-mor_struct[2]
        curr_data[, 'MOR.SEQUENCE']<-mor_struct[3]
        curr_data[, 'EXTENDED.SEQUENCE']<-mor_struct[4]
        curr_data[, 'SECONDARY.STRUCTURE']<-mor_struct[5]
        for(i in 1:nchar(mor_struct[1])){
          if(i!=0){
            curr_data[, toString(i)]<-substr(mor_struct[1], i, i) 
          }
        }
        if(nrow(curr_data)>1){print(extended_sequences)}
        data<-rbind.fill(data, curr_data)
        curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
        
      }
    }else{
      next
    }
  }
  #data<-data[complete.cases(data),]
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

get5PMorSecondaryStructureListFromFasta<-function(fasta){
  
  primirna_list<-fasta$PRIMIRNA
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, stringsAsFactors = F)
  
  for(p in primirna_list){
    m<-get3PMirna(p)
    extended_sequences<-c()
    extended_sequences<-getExtendedSequence(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      #edit extended sequence to sequence with maximum brackets in 2d structure
      secondary_struct<-toString(fasta[which(fasta$PRIMIRNA==p), 'SECONDARY.STRUCTURE'])
      mor_struct<-substr(secondary_struct, 1, 11)
      sequence<-toString(fasta[which(fasta$PRIMIRNA==p), 'SEQUENCE'])
      mor_seq<-substr(sequence, 1, 11)
      curr_data[, 'MOR.EDITED']<-mor_struct
      curr_data[, 'MOR']<-mor_struct
      curr_data[, 'MOR.SEQUENCE']<-mor_seq
      for(i in 1:nchar(mor_struct)){
        if(i!=0){
          curr_data[, toString(i)]<-substr(mor_struct, i, i) 
        }
      }
      data<-rbind.fill(data, curr_data)
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
    }else{
      next
    }
  }
  #data<-data[complete.cases(data),]
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

get3PMorSecondaryStructureListFromFasta<-function(fasta){
  
  primirna_list<-fasta$PRIMIRNA
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, stringsAsFactors = F)
  
  for(p in primirna_list){
    m<-get3PMirna(p)
    extended_sequences<-c()
    extended_sequences<-getExtendedSequence(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      #edit extended sequence to sequence with maximum brackets in 2d structure
      secondary_struct<-toString(fasta[which(fasta$PRIMIRNA==p), 'SECONDARY.STRUCTURE'])
      mor_struct<-substr(secondary_struct, nchar(secondary_struct)-13, nchar(secondary_struct))
      sequence<-toString(fasta[which(fasta$PRIMIRNA==p), 'SEQUENCE'])
      mor_seq<-substr(sequence, nchar(secondary_struct)-13, nchar(secondary_struct))
      curr_data[, 'MOR.EDITED']<-mor_struct
      curr_data[, 'MOR']<-mor_struct
      curr_data[, 'MOR.SEQUENCE']<-mor_seq
      for(i in 1:nchar(mor_struct)){
        if(i!=0){
          curr_data[, toString(i)]<-substr(mor_struct, i, i) 
        }
      }
      data<-rbind.fill(data, curr_data)
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
    }else{
      next
    }
  }
  #data<-data[complete.cases(data),]
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

getMaxBracket5PMorSecondaryStructureList<-function(primirna_list){
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MAX.BRACKET=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, stringsAsFactors = F)
  for(p in primirna_list){
    m<-get5PMirna(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      extended_sequences<-getMaximumBracketsExtendedSequence(p)
      if(nchar(extended_sequences)>0){
        curr_data[, 'MAX.BRACKET']<-T
      }else{
        next
        curr_data[, 'MAX.BRACKET']<-F
        extended_sequences<-getExtendedSequence(p)
        pre_seq<-getPrecursorSequence(p)
        if(is.na(pre_seq)){
          print(paste('cannot retrieve precursor seq: ', p))
        }else{
          pre_start<-gregexpr(pre_seq, extended_sequences, fixed=T)[[1]]
          start<-pre_start-30
          end<-pre_start+nchar(pre_seq)-1+30
          extended_sequences<-substr(extended_sequences, start, end)
        }
      }
      mor_struct<-get5PMorSecondaryStructure(p, toString(extended_sequences))
      if(!any(is.na(mor_struct))){
        curr_data[, 'MOR.EDITED']<-mor_struct[1]
        curr_data[, 'MOR']<-mor_struct[2]
        curr_data[, 'MOR.SEQUENCE']<-mor_struct[3]
        for(i in 1:nchar(mor_struct[1])){
          curr_data[, toString(-(nchar(mor_struct[1])-i+1))]<-substr(mor_struct[1], i, i)
        }
        data<-rbind.fill(data, curr_data)
        curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
      }
    }else{
      next
    }
  }
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

getMaxBracket3PMorSecondaryStructureList<-function(primirna_list){
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, MAX.BRACKET=NA, stringsAsFactors = F)
  for(p in primirna_list){
    print(p)
    m<-get3PMirna(p)
    extended_sequences<-c()
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      extended_sequences<-getMaximumBracketsExtendedSequence(p)
      if(nchar(extended_sequences)>0){
        curr_data[, 'MAX.BRACKET']<-T
      }else{
        next
        curr_data[, 'MAX.BRACKET']<-F
        extended_sequences<-getExtendedSequence(p)
        pre_seq<-getPrecursorSequence(p)
        if(is.na(pre_seq)){
          print(paste('cannot retrieve precursor seq: ', p))
        }else{
          pre_start<-gregexpr(pre_seq, extended_sequences, fixed=T)[[1]]
          start<-pre_start-30
          end<-pre_start+nchar(pre_seq)-1+30
          extended_sequences<-substr(extended_sequences, start, end)
        }
      }
      #edit extended sequence to sequence with maximum brackets in 2d structure
      mor_struct<-get3PMorSecondaryStructure(p, toString(extended_sequences))
      if(!any(is.na(mor_struct))){
        curr_data[, 'MOR.EDITED']<-mor_struct[1]
        curr_data[, 'MOR']<-mor_struct[2]
        curr_data[, 'MOR.SEQUENCE']<-mor_struct[3]
        for(i in 1:nchar(mor_struct[1])){
          if(i!=0){
            curr_data[, toString(i)]<-substr(mor_struct[1], i, i) 
          }
        }
        if(nrow(curr_data)>1){print(extended_sequences)}
        data<-rbind.fill(data, curr_data)
        curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p))
        
      }
    }else{
      next
    }
  }
  #data<-data[complete.cases(data),]
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

get5PBracketSD<-function(mor_struct_list){
  number_of_base<-ncol(mor_struct_list)-6
  data<-data.frame()
  for(i in 1:400){
    positions<-round(runif(25, 1, nrow(mor_struct_list)))
    temp_data<-mor_struct_list[positions,]
    for(x in 1:number_of_base){
      curr_pos<-(number_of_base-x+1)*-1
      curr_bases<-toString(paste(temp_data[,toString(curr_pos)], collapse=''))
      bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
      data[i, toString(curr_pos)]<-bracket/nchar(curr_bases)
    }
  }
  mor_file<-paste(base_path, 'miRBase21-5p-moR-subset.tsv', sep='')
  write.table(data, file=mor_file, sep='\t', row.names = F)
  sd<-apply(data, 2, sd)
  means<-colMeans(data)
  max<-data
  min<-data
  for(i in 1:struct_len){
    curr_col<-as.vector(data[,toString(-i)])
    curr_col<-sort(curr_col)
    #max[,toString(-i)]<-NULL
    max[1,toString(-i)]<-curr_col[length(curr_col)-9]
    #min[,toString(-i)]<-NULL
    min[1,toString(-i)]<-curr_col[10]
  }
  max<-max[1,]
  min<-min[1,]
  #max<-means+(sd*2)#sapply(data, max, na.rm = TRUE)
  #min<-means-(sd*2)#sapply(data, min, na.rm = TRUE)
  result<-data.frame()
  result<-rbind(means, sd, max, min)
  result<-cbind(result, BASE=c('BRACKET.SAMPLE', 'SD.SAMPLE','MAX.SAMPLE', 'MIN.SAMPLE'))
  return(as.data.frame(result))
}

get3PBracketSD<-function(mor_struct_list){
  number_of_base<-ncol(mor_struct_list)-6
  data<-data.frame()
  for(i in 1:400){
    positions<-round(runif(25, 1, nrow(mor_struct_list)))
    temp_data<-mor_struct_list[positions,]
    for(x in 1:number_of_base){
      curr_bases<-toString(paste(temp_data[,toString(x)], collapse=''))
      bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
      data[i, toString(x)]<-bracket/nchar(curr_bases)
    }
  }
  mor_file<-paste(base_path, 'miRBase21-3p-moR-subset.tsv', sep='')
  write.table(data, file=mor_file, sep='\t', row.names = F)
  means<-colMeans(data)
  max<-data
  min<-data
  for(i in 1:struct_len){
    curr_col<-as.vector(data[,toString(i)])
    curr_col<-sort(curr_col)
    max[1,toString(i)]<-curr_col[length(curr_col)-9]
    min[1,toString(i)]<-curr_col[10]
  }
  max<-max[1,]
  min<-min[1,]
  sd<-apply(data, 2, sd)
  #max<-means+(sd*2)#sapply(data, max, na.rm = TRUE)
  #min<-means-(sd*2)#sapply(data, min, na.rm = TRUE)
  result<-data.frame()
  result<-rbind(means, sd, max, min)
  result<-cbind(result, BASE=c('BRACKET.SAMPLE', 'SD.SAMPLE', 'MAX.SAMPLE', 'MIN.SAMPLE'))
  return(as.data.frame(result))
}

get5PBasePairRatioPerPosition<-function(mor_struct_list){ #dataframe with miRNA and mor secondary structs
  data<-data.frame(BASE=c('DOT', 'OTHER.STEM', 'BRACKET'), stringsAsFactors = F)
  number_of_base<-ncol(mor_struct_list) - grep("-30", colnames(mor_struct_list), fixed=T)[1]
  temp_mor_struct_list<-mor_struct_list
  for(i in 1:number_of_base){
    curr_pos<-(number_of_base-i+1)*-1
    curr_bases<-toString(paste(mor_struct_list[,toString(curr_pos)], collapse=''))
    dot<-str_count(curr_bases, '\\.')
    bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
    asterik<-str_count(curr_bases, '\\*')
    data[which(data$BASE=='DOT'), toString(curr_pos)]<-dot/nchar(curr_bases)
    data[which(data$BASE=='BRACKET'), toString(curr_pos)]<-bracket/nchar(curr_bases)
    data[which(data$BASE=='OTHER.STEM'), toString(curr_pos)]<-asterik/nchar(curr_bases)
    for(x in 1:nrow(temp_mor_struct_list)){
      if(temp_mor_struct_list[x, toString(curr_pos)]%in%'*'){
        temp_mor_struct_list[x,]<-''
      }
    }
    curr_bases<-toString(paste(temp_mor_struct_list[,toString(curr_pos)], collapse=''))
    bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
    data[which(data$BASE=='BRACKET.NO.ASTERIK'), toString(curr_pos)]<-bracket/nchar(curr_bases)
  }
  #data[,'SD']<-as.data.frame(transform(data, SD=apply(data, 1, sd, na.rm = TRUE)))[,'SD']
  return(data)
}

get3PBasePairRatioPerPosition<-function(mor_struct_list){ #dataframe with miRNA and mor secondary structs
  data<-data.frame(BASE=c('DOT', 'OTHER.STEM', 'BRACKET'), stringsAsFactors = F)
  number_of_base<-ncol(mor_struct_list) - grep("1", colnames(mor_struct_list))[1]
  temp_mor_struct_list<-mor_struct_list
  for(i in 1:number_of_base){
    curr_bases<-toString(paste(mor_struct_list[,toString(i)], collapse=''))
    dot<-str_count(curr_bases, '\\.')
    bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
    asterik<-str_count(curr_bases, '\\*')
    data[which(data$BASE=='DOT'), toString(i)]<-dot/nchar(curr_bases)
    data[which(data$BASE=='OTHER.STEM'), toString(i)]<-asterik/nchar(curr_bases)
    data[which(data$BASE=='BRACKET'), toString(i)]<-bracket/nchar(curr_bases)
    for(x in 1:nrow(temp_mor_struct_list)){
      if(temp_mor_struct_list[x, i]%in%'*'){
        temp_mor_struct_list[x,]<-''
      }
    }
    curr_bases<-toString(paste(temp_mor_struct_list[,toString(i)], collapse=''))
    bracket<-str_count(curr_bases, '\\)')+str_count(curr_bases, '\\(')
    data[which(data$BASE=='BRACKET.NO.ASTERIK'), toString(i)]<-bracket/nchar(curr_bases)
  }
  #data[,'SD']<-transform(data, SD=apply(data, 1, sd, na.rm = TRUE))['SD']
  return(data)
}

get5PBasePairData<-function(high_fid_mirna, mid_fid_mirna, low_fid_mirna){
  high_fid_mor_struct<-get5PMorSecondaryStructureList(high_fid_mirna)
  high_fid_base_pair_ratio<-get5PBasePairRatioPerPosition(high_fid_mor_struct)
  for(i in 1:nrow(high_fid_base_pair_ratio)){
    high_fid_base_pair_ratio[i,1]<-toString(paste('HIGH.FID.', toString(high_fid_base_pair_ratio[i,1]), sep=''))
  }
  
  mid_fid_mor_struct<-get5PMorSecondaryStructureList(mid_fid_mirna)
  mid_fid_base_pair_ratio<-get5PBasePairRatioPerPosition(mid_fid_mor_struct)
  for(i in 1:nrow(mid_fid_base_pair_ratio)){
    mid_fid_base_pair_ratio[i,1]<-toString(paste('MID.FID.', toString(mid_fid_base_pair_ratio[i,1]), sep=''))
  }
  
  low_fid_mor_struct<-get5PMorSecondaryStructureList(low_fid_mirna)
  low_fid_base_pair_ratio<-get5PBasePairRatioPerPosition(low_fid_mor_struct)
  for(i in 1:nrow(low_fid_base_pair_ratio)){
    low_fid_base_pair_ratio[i,1]<-toString(paste('LOW.FID.', toString(low_fid_base_pair_ratio[i,1]), sep=''))
  }
  base_pair_ratio<-rbind(low_fid_base_pair_ratio, mid_fid_base_pair_ratio, high_fid_base_pair_ratio)
  
  file<-paste(base_path, 'miRNA-5p-moR-pairing-by-fidelity.tsv', sep='')

  write.table(base_pair_ratio, file=file, sep='\t', row.names=F)
}

get3PBasePairData<-function(high_fid_mirna, mid_fid_mirna, low_fid_mirna){
  high_fid_mor_struct<-get3PMorSecondaryStructureList(high_fid_mirna)
  high_fid_base_pair_ratio<-get3PBasePairRatioPerPosition(high_fid_mor_struct)
  for(i in 1:nrow(high_fid_base_pair_ratio)){
    high_fid_base_pair_ratio[i,1]<-toString(paste('HIGH.FID.', toString(high_fid_base_pair_ratio[i,1]), sep=''))
  }

  mid_fid_mor_struct<-get3PMorSecondaryStructureList(mid_fid_mirna)
  mid_fid_base_pair_ratio<-get3PBasePairRatioPerPosition(mid_fid_mor_struct)
  for(i in 1:nrow(mid_fid_base_pair_ratio)){
    # mid_fid_base_pair_ratio[i,1]<-toString(paste('MID.FID.', toString(mid_fid_base_pair_ratio[i,1]), sep=''))
  }
  
  low_fid_mor_struct<-get3PMorSecondaryStructureList(low_fid_mirna)
  low_fid_base_pair_ratio<-get3PBasePairRatioPerPosition(low_fid_mor_struct)
  for(i in 1:nrow(low_fid_base_pair_ratio)){
    low_fid_base_pair_ratio[i,1]<-toString(paste('LOW.FID.2', toString(low_fid_base_pair_ratio[i,1]), sep=''))
  }
  
  base_pair_ratio<-rbind(low_fid_base_pair_ratio, mid_fid_base_pair_ratio, high_fid_base_pair_ratio)
  
  file<-paste(base_path, 'miRNA-3p-moR-pairing-by-fidelity.tsv', sep='')
  
  write.table(base_pair_ratio, file=file, sep='\t', row.names=F)

}

get5PMorStructFile<-function(mor_struct_list){
  result<-c()
  for(i in 1:nrow(mor_struct_list)){
    mirna<-mor_struct_list[i,'MIRNA']
    mor_seq<-mor_struct_list[i,'MOR.SEQUENCE']
    mor_seq<-substr(mor_seq, 18, 30)
    mor_struct<-toString(mor_struct_list[i, 'MOR'])
    mor_struct<-substr(mor_struct, 18, 30)
    result<-c(result, c(paste('>', mirna, sep=''), mor_seq, mor_struct))
  }
  return(result)
}

get3PMorStructFile<-function(mor_struct_list){
  result<-c()
  for(i in 1:nrow(mor_struct_list)){
    mirna<-mor_struct_list[i,'MIRNA']
    mor_seq<-mor_struct_list[i,'MOR.SEQUENCE']
    mor_seq<-substr(mor_seq, 1, 11)
    mor_struct<-toString(mor_struct_list[i, 'MOR'])
    mor_struct<-substr(mor_struct, 1, 11)
    result<-c(result, c(paste('>', mirna, sep=''), mor_seq, mor_struct))
  }
  return(result)
}

getRefinedJunctionSecondaryStructure<-function(primirna){
  primirna<-getPriMirna(primirna)
  extended_seq<-toString(getExtendedSequence(primirna))
  pre_seq<-toString(getPrecursorSequence(primirna))
  pre_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  start<-conserved_struct_data[which(conserved_struct_data$PRIMIRNA==primirna),'START']*-1
  end<-conserved_struct_data[which(conserved_struct_data$PRIMIRNA==primirna),'END']
  start<-pre_start-start
  end<-pre_start+nchar(pre_seq)+end-1
  
  edited_extended_seq<-substr(extended_seq, start, end)
  edited_extended_struct<-toString(getSecondaryStructure(edited_extended_seq, primirna))
  edited_pre_start<-gregexpr(pre_seq, toString(edited_extended_seq), fixed=T)[[1]]
  junction_start<-edited_pre_start-11
  junction_end<-edited_pre_start+nchar(pre_seq)+13-1
  junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
  n_left_bracket<-str_count(junction_struct, '\\(')
  n_right_bracket<-str_count(toString(junction_struct), '\\)')
  
  if(n_left_bracket != n_right_bracket){
    while(n_left_bracket<n_right_bracket){
      junction_start<-junction_start-1
      junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
      n_left_bracket<-str_count(junction_struct, '\\(')
    }
    while(n_right_bracket<n_left_bracket){
      junction_end<-junction_end+1
      junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
      n_right_bracket<-str_count(junction_struct, '\\)')
    }
  }
  junction_seq<-substr(edited_extended_seq, junction_start, junction_end)
  #remove 5p bulges
  junction_pre_start<-gregexpr(pre_seq, toString(junction_seq), fixed=T)[[1]]
  mor_5p<-substr(junction_struct, 1, junction_pre_start-1)
  bulge_end<-gregexpr(')', mor_5p, fixed=T)[[1]]
  if(bulge_end!=-1){
    bulge_end<-bulge_end[length(bulge_end)]
    junction_start<-bulge_end+1
  }
  
  #remove 3p bulges
  junction_pre_end<-junction_pre_start+nchar(pre_seq)
  mor_3p<-substr(junction_struct, junction_pre_end, nchar(junction_struct))
  bulge_start<-gregexpr('(', mor_3p, fixed=T)[[1]][1]
  if(bulge_start!=-1){
    junction_end<-junction_pre_end+bulge_start-1-1
  }
  
  junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
  junction_seq<-substr(edited_extended_seq, junction_start, junction_end)
  
  n_left_bracket<-str_count(junction_struct, '\\(')
  n_right_bracket<-str_count(toString(junction_struct), '\\)')
  
  if(n_left_bracket != n_right_bracket){
    while(n_left_bracket<n_right_bracket){
      junction_start<-junction_start-1
      junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
      n_left_bracket<-str_count(junction_struct, '\\(')
    }
    while(n_right_bracket<n_left_bracket){
      junction_end<-junction_end+1
      junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
      n_right_bracket<-str_count(junction_struct, '\\)')
    }
  }
  junction_struct<-toString(substr(edited_extended_struct, junction_start, junction_end))
  junction_seq<-substr(edited_extended_seq, junction_start, junction_end)
  
  #remove unpaired bases at end
  junction_start<-gregexpr('(', junction_struct, fixed=T)[[1]][1]
  junction_end<-gregexpr(')', junction_struct, fixed=T)[[1]]
  junction_end<-junction_end[length(junction_end)]

  junction_struct<-toString(substr(junction_struct, junction_start, junction_end))
  junction_seq<-substr(junction_seq, junction_start, junction_end)
  
  junction_pre_start<-gregexpr(pre_seq, junction_seq, fixed=T)[[1]]
  added_start<-junction_pre_start-1
  added_end<-nchar(junction_seq)-added_start-nchar(pre_seq)
  return(c(junction_seq, junction_struct, -added_start, added_end))
}

getRefinedJunctionSecondaryStructureData<-function(mirna_list, posi_file){
  result<-c()
  posi_data<-data.frame(PRIMIRNA=NA, ADDED.START=NA, ADDED.END=NA)
  for(m in unique(mirna_list)){
    primirna<-getPriMirna(m)
    for(p in primirna){
      secondary_struct<-getRefinedJunctionSecondaryStructure(p)
      posi_data<-rbind(posi_data, data.frame(PRIMIRNA=p, ADDED.START=secondary_struct[3], ADDED.END=secondary_struct[4]))
      p<-gsub('-', 'X', p)
      print(secondary_struct[3:4])
      #p<-paste(p, paste(secondary_struct[3:4], collapse='X'), sep='X')
      result<-c(result, paste('>', p, sep=''))
      result<-c(result, secondary_struct[1:2])
    }
  }
  write.table(posi_data, file=posi_file, sep='\t', row.names=F)
  return(unique(result))
}

getEditedExtendedSecondaryStructureData<-function(mirna_list){
  result<-c()
  for(m in unique(mirna_list)){
    primirna<-getPriMirna(m)
    for(p in primirna){
      secondary_struct<-getEditedExtendedSecondaryStructure(p)
      p<-gsub('-', '0', p)
      result<-c(result, paste('>', p, sep=''))
      result<-c(result, secondary_struct)
    }
  }
  return(unique(result))
}

getPrimarySecondaryStructure<-function(primirna){
  pri_seq<-getPrimarySequence(primirna)
  pri_secondary_struct<-getSecondaryStructure(pri_seq, primirna)
  return(pri_secondary_struct)
}

#get primary secondary structure without unpaired bases at the 5p/3p ends
getEditedPrimarySecondaryStructure<-function(primirna){
  pri_seq<-getPrimarySequence(p)
  pri_secondary_struct<-getSecondaryStructure(pri_seq, primirna)
  left_bracket_pos<-gregexpr('(', pri_secondary_struct, fixed=T)[[1]]
  first_left_bracket_pos<-left_bracket_pos[1]
  right_bracket_pos<-gregexpr(')', pri_secondary_struct, fixed=T)[[1]]
  last_right_bracket_pos<-right_bracket_pos[length(right_bracket_pos)]
  pri_secondary_struct<-substr(pri_secondary_struct, first_left_bracket_pos, last_right_bracket_pos)
  return(pri_secondary_struct)
}

getPrecursorSecondaryStructure<-function(primirna){
  seq_5p<-get5PSequence(primirna)
  seq_3p<-get3PSequence(primirna)
  pri_seq<-getPrimarySequence(primirna)
  pri_secondary_struct<-getSecondaryStructure(pri_seq, primirna)
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
  if(length(start)==0 || length(end)==0){
    return(NA)
  }
  pre_secondary_struct<-substr(pri_secondary_struct, start, end)
  return(pre_secondary_struct)
}

hasJunctionMirbase<-function(primirna){
  pri_seq<-getPrimarySequence(primirna)
  seq_5p<-get5PSequence(primirna)
  seq_3p<-get3PSequence(primirna)
  
  start_5p<-gregexpr(seq_5p, pri_seq, fixed=T)[[1]]
  end_3p<-gregexpr(seq_3p, pri_seq, fixed=T)[[1]]+nchar(seq_3p)
  
  if(start_5p>11 & end_3p<(nchar(pri_seq)-12)){
    return(TRUE)
  }
  return(FALSE)
}
  
get5PJunction<-function(primirna){
  extended_seq<-getExtendedSequence(primirna)
  mirna_5p<-get5PMirna(primirna)
  struct_5p<-get5PMorSecondaryStructure(primirna, extended_seq)[1]
  junction<-gregexpr('*', struct_5p, fixed=T)[[1]]
  junction<-junction[length(junction)]
  if(junction!=-1){
    junction<-nchar(struct_5p)-junction
  }else{
    struct_3p<-get3PMorSecondaryStructure(primirna, getExtendedSequence(primirna))[1]
    junction_3p<-gregexpr('*', struct_3p, fixed=T)[[1]]
    if(junction_3p[1]!=-1){
      junction_3p<-junction[1]-1
      struct_3p<-substr(struct_3p, 1, junction)
      padding<-str_count(struct_3p, '\\)')
      junction<-gregexpr('(', struct_5p, fixed=T)[[1]][padding+1]+1
    }else{
      junction<-NA
    }
  }
  return(-junction)
}

get3PJunction<-function(primirna){
  extended_seq<-getExtendedSequence(primirna)
  mirna_3p<-get3PMirna(primirna)
  struct_3p<-get3PMorSecondaryStructure(primirna, extended_seq)[1]
  junction<-gregexpr('*', struct_3p, fixed=T)[[1]][1]
  if(junction!=-1){
    junction<-junction-1
  }else{
    struct_5p<-get5PMorSecondaryStructure(primirna, getExtendedSequence(primirna))[1]
    junction_5p<-gregexpr('*', struct_5p, fixed=T)[[1]]
    if(junction_5p[1]!=-1){
      junction_5p<-junction[length(junction)]+1
      struct_5p<-substr(struct_5p, junction, nchar(struct_5p))
      padding<-str_count(struct_5p, '\\(')
      junction<-gregexpr(')', struct_3p, fixed=T)[[1]][padding+1]-1
    }else{
      junction<-NA
    }
    
  }
  return(junction)
}

getBrackets5PJunctionToPrecursor<-function(primirna){
  struct_5p<-get5PMorSecondaryStructure(primirna, getExtendedSequence(primirna))
  junction<-gregexpr('*', struct_5p, fixed=T)[[1]]
  if(junction[1]==-1){return(NA)}
  junction<-junction[length(junction)]+1
  mor_to_pre<-substr(struct_5p, junction, nchar(struct_5p))
  return(str_count(mor_to_pre, '\\('))
}

getBrackets3PJunctionToPrecursor<-function(primirna){
  struct_3p<-get3PMorSecondaryStructure(primirna, getExtendedSequence(primirna))
  junction<-gregexpr('*', struct_3p, fixed=T)[[1]]
  if(junction[1]==-1){return(NA)}
  junction<-junction[1]-1
  mor_to_pre<-substr(struct_3p, 1, junction)
  return(str_count(mor_to_pre, '\\)'))
}

get5PMorBrackets<-function(pre_struct, pri_struct){
  pre_start<-gregexpr(pre_struct, pri_struct, fixed=T)[[1]]
  mor_5p<-substr(pri_struct, 1, pre_start-1)
  num_brackets<-str_count(mor_5p, '\\(')
  return(num_brackets)
}

get3PMorBrackets<-function(pre_struct, pri_struct){
  pre_start<-gregexpr(pre_struct, pri_struct, fixed=T)[[1]]
  mor_3p<-substr(pri_struct, pre_start+nchar(pre_struct), nchar(pri_struct))
  num_brackets<-str_count(mor_3p, '\\)')
  return(num_brackets)
}

writeFile<-function(data, filename){
  file<-file(filename)
  writeLines(data, file)
  close(file)
}

getFastaData<-function(fasta_file){
  fasta_data<-read.table(file=fasta_file)
  data<-data.frame(PRIMIRNA=NA, SEQUENCE=NA, SECONDARY.STRUCTURE=NA)
  for(i in seq(from=1, to=nrow(fasta_data), by=3)){
    primirna<-toString(fasta_data[i,])
    if(gregexpr('corr', toString(primirna), fixed=T)[[1]]!=-1){
      primirna<-substr(primirna, 2, nchar(primirna))
      primirna<-substr(primirna, 1, nchar(primirna)-4)
      sequence<-toString(fasta_data[i+1,])
      struct<-toString(fasta_data[i+2,])
      data<-rbind(data, data.frame(PRIMIRNA=primirna,
                                   SEQUENCE = sequence,
                                   SECONDARY.STRUCTURE=struct))
    }
  }
  data<-data[-1,]
  return(data)
}


#compare secondary structure from miRBase21 with UCSC extended sequences
data<-data.frame(PRIMIRNA=NA, 
                 PRIMARY.SECONDARY.STRUCTURE=NA, 
                 EXTENDED.SECONDARY.STRUCTURE=NA, 
                 PRIMARY.IN.EXTENDED=NA,
                 PRE.IN.EXT=NA,
                 HAS.JUNCTION.MIRBASE=NA, 
                 EXTENDED.JUNCTION.5P=NA,
                 EXTENDED.JUNCTION.3P=NA,
                 PAD.BRACKETS.5P=NA,
                 PAD.BRACKETS.3P=NA,
                 FIDELITY=NA)
for(m in c(low_fid_mirna, high_fid_mirna)){
  primirna<-getPriMirna(m)
  for(p in primirna){
    extended_struct<-getEditedExtendedSecondaryStructure(p)[2]
    primary_struct<-getEditedPrimarySecondaryStructure(p)
    pri_in_extended<-F
    if(gregexpr(primary_struct, extended_struct, fixed=T)[[1]]!=-1){
      pri_in_extended<-T
    }
    pre_struct<-getPrecursorSecondaryStructure(p)
    pre_in_pri<-F
    if(gregexpr(pre_struct, primary_struct, fixed=T)[[1]]!=-1){
      pre_in_pri<-T
    }
    pre_in_ext<-F
    if(gregexpr(pre_struct, extended_struct, fixed=T)[[1]]!=-1){
      pre_in_ext<-T
    }
    
    has_junction<-hasJunctionMirbase(p)
    junction_5p<-get5PJunction(p)
    junction_3p<-get3PJunction(p)
    pad_brackets_5p<-getBrackets5PMorToPrecursor(p)
    pad_brackets_3p<-getBrackets3PMorToPrecursor(p)
    
    fidelity='HIGH'
    if(m %in% low_fid_mirna){fidelity='LOW'}
    
    data<-rbind(data, data.frame(PRIMIRNA=p, 
                                 PRIMARY.SECONDARY.STRUCTURE=primary_struct, 
                                 EXTENDED.SECONDARY.STRUCTURE=extended_struct,
                                 PRE.IN.EXT=pre_in_ext,
                                 PRIMARY.IN.EXTENDED=pri_in_extended,
                                 HAS.JUNCTION.MIRBASE=has_junction,
                                 EXTENDED.JUNCTION.5P=junction_5p,
                                 EXTENDED.JUNCTION.3P=junction_3p,
                                 PAD.BRACKETS.5P=pad_brackets_5p,
                                 PAD.BRACKETS.3P=pad_brackets_3p,
                                 FIDELITY=fidelity))
  }
}

data<-data[-1,]

write.table(data, file=paste(base_path, 'primary-extended-secondary-structure.tsv', sep=''), sep='\t', row.names=F)



get5PBasePairData(high_fid_mirna, mid_fid_mirna, low_fid_mirna)
get3PBasePairData(high_fid_mirna, mid_fid_mirna, low_fid_mirna)

#GET MOR WITH SECONDARY STRUCTURE FROM EDITED FASTA [-11, 13]
fasta_file<-paste0(base_path, 'low-fid-mirna-edited-seq.txt')
fasta_data<-getFastaData(fasta_file)


#GET MOR WITH SECONDARY STRUCTURE
# unique_primirna<-unique(mirna_data$PRIMIRNA)
# unique_primirna<-unique_primirna[!is.na(unique_primirna)]


# mor_struct_3p<-get3PMorSecondaryStructureList(unique_primirna)
# mor_struct_file<-paste0(base_path, 'miRBase21-moR-struct-3p-final.tsv')
# write.table(mor_struct_3p, file=mor_struct_file, sep='\t', row.names=F)
# mor_struct_3p[,'0']<-NULL
# base_pair_ratio_3p<-get3PBasePairRatioPerPosition(mor_struct_3p)
# #base_pair_ratio_3p<-rbind(base_pair_ratio_3p, get3PBracketSD(mor_struct_3p))
# 
# file<-paste(base_path, 'miRBase21-moR-base-pairing-3p-final.tsv', sep='')
# #
# write.table(base_pair_ratio_3p, file=file, sep='\t', row.names=F)
#
#
# mor_struct_5p<-get5PMorSecondaryStructureList(unique_primirna)
# mor_struct_file<-paste0(base_path, 'miRBase21-moR-struct-5p-final.tsv')
# write.table(mor_struct_5p, file=mor_struct_file, sep='\t', row.names=F)
# base_pair_ratio_5p<-get5PBasePairRatioPerPosition(mor_struct_5p)
#base_pair_ratio_5p<-rbind(base_pair_ratio_5p, get5PBracketSD(mor_struct_5p))
#
# file<-paste(base_path, 'miRBase21-moR-base-pairing-5p-final.tsv', sep='')
#
# write.table(base_pair_ratio_5p, file=file, sep='\t', row.names=F)

#GETTING FA FILE FOR PRIMARY SEQUENCES
# low_posi_file<-paste0(base_path, 'low fidelity tertiary structure [-11, 13] (max brackets in 2d)/low-fidelity-secondary-structure-position.tsv')
# low_mor_struct<-getRefinedJunctionSecondaryStructureData(low_fid_mirna, low_posi_file)
# filename<-paste(base_path, 'low-fidelity-secondary-struct.txt', sep='')
# writeFile(low_mor_struct, filename)
# 
# high_posi_file<-paste0(base_path, 'high fidelity tertiary structure [-11, 13] (max brackets in 2d)/high-fidelity-secondary-structure-position.tsv')
# high_mor_struct<-getRefinedJunctionSecondaryStructureData(high_fid_mirna, high_posi_file)
# filename<-paste(base_path, 'high-fidelity-secondary-struct.txt', sep='')
# writeFile(high_mor_struct, filename)
