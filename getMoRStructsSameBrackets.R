
## PURPOSE: get 5P and 3P moR 2d structs with the same number of pairing
## INPUT: miRBase data                                                          miRBase21-master.tsv
##        conserved extended primiRNA 2d structs  extended-2d-struct-max-moR-pairing-all-miRBase.tsv
## OUTPUT: table containing moR with even brackets      low-fidelity-moR-struct-3p-even-brackets.tsv
require('stringr')
require('plyr')

base_path<-'/Users/chens22/Documents/miRNA/structure/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors=F)
conserved_struct_file<-'/Users/chens22/Documents/miRNA/structure/extended-2d-struct-max-moR-pairing-all-miRBase.tsv'
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
  
  pre_seq<-substr(pri_seq, start, end)
  if(is.na(pre_seq)){
    return('')
  }
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
  extended_seq<-substr(extended_seq, pri_seq_start-struct_len, pri_seq_start+nchar(pri_seq)+struct_len-1)
  secondary_struct<-getSecondaryStructure(extended_seq, mirna)#getPriSecondaryStructure(mirna)
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
  print(paste(mor_secondary_struct, collapse=''))
  mor_secondary_struct<-mor_secondary_struct[(length(mor_secondary_struct)-(struct_len-1)):length(mor_secondary_struct)]
  temp_mor_secondary_struct<-temp_mor_secondary_struct[(length(temp_mor_secondary_struct)-(struct_len-1)):length(temp_mor_secondary_struct)]
  extended_seq<-substr(extended_seq, length(temp_mor_secondary_struct)-(struct_len-1), length(temp_mor_secondary_struct))
  return(c(paste(mor_secondary_struct, collapse=''), paste(temp_mor_secondary_struct, collapse=''), extended_seq))
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
  extended_seq<-substr(extended_seq, pri_seq_start-struct_len, pri_seq_start+nchar(pri_seq)+struct_len-1)
  secondary_struct<-getSecondaryStructure(extended_seq, mirna)
  pre_seq<-getPrecursorSequence(primirna)
  if(nchar(pre_seq)<1){
    print(paste0('NO PRECURSOR:', primirna))
    return('')
  }
  seq_start<-gregexpr(pre_seq, extended_seq, fixed=TRUE)[[1]]
  if(length(seq_start)>1){
    pri_seq_start<-gregexpr(pri_seq, extended_seq, fixed=T)[[1]]
    seq_start_pri<-gregexpr(seq, pri_seq, fixed=T)[[1]]
    seq_start<-pri_seq_start+seq_start_pri
  }
  
  
  bases_len<-0
  mor_secondary_struct<-c()
  temp_mor_secondary_struct<-c()
  i<-0
  while(seq_start+nchar(pre_seq)+i<=nchar(secondary_struct) & bases_len<struct_len){
    curr_base<-substr(secondary_struct, seq_start+nchar(pre_seq)+i, seq_start+nchar(pre_seq)+i)
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
  extended_seq<-substr(extended_seq, 1, struct_len)
  
  return(c(paste(mor_secondary_struct, collapse=''), paste(temp_mor_secondary_struct, collapse=''), extended_seq))
}


get5PMorSecondaryStructureList<-function(primirna_list){
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, BRACKETS=NA, stringsAsFactors = F)
  for(p in primirna_list){
    m<-get5PMirna(p)
    if(nchar(m)>0){
      curr_data<-data.frame(MIRNA=toString(m), PRIMIRNA=toString(p), stringsAsFactors=F)
      extended_sequences<-getExtendedSequence(p)
      mor_struct<-get5PMorSecondaryStructure(p, toString(extended_sequences))
      if(!any(is.na(mor_struct))){
        curr_data[, 'MOR.EDITED']<-mor_struct[1]
        curr_data[, 'MOR']<-mor_struct[2]
        curr_data[, 'MOR.SEQUENCE']<-mor_struct[3]
        curr_data[, 'BRACKETS']<-str_count(mor_struct[1], pattern='\\(')
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
  data<-data.frame(MIRNA=NA, PRIMIRNA=NA, MOR.EDITED=NA, MOR.SEQUENCE=NA, BRACKETS=NA, stringsAsFactors = F)
  count<-1
  for(p in primirna_list){
    count-count+1
    print(p)
    print(count)
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
        curr_data[, 'BRACKETS']<-str_count(mor_struct[1], pattern='\\)')
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
      print(paste0('NO 3P PRIMIRNA:', p))
      next
    }
  }
  #data<-data[complete.cases(data),]
  data<-data[-1,]
  data<-unique(data)
  return(data)
}

get5PBasePairRatioPerPosition<-function(mor_struct_list){ #dataframe with miRNA and mor secondary structs
  data<-data.frame(BASE=c('DOT', 'OTHER.STEM', 'BRACKET'), stringsAsFactors = F)
  number_of_base<-ncol(mor_struct_list)-grep('-30', colnames(mor_struct_list))[1]-1
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
  number_of_base<-ncol(mor_struct_list)-grep('1', colnames(mor_struct_list))[1]-1
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


unique_primirna<-unique(mirna_data$PRIMIRNA)

mor_struct_3p<-get3PMorSecondaryStructureList(low_fid_mirna)
mor_struct_5p<-get5PMorSecondaryStructureList(low_fid_mirna)

unique_primirna<-intersect(mor_struct_3p$PRIMIRNA, mor_struct_5p$PRIMIRNA)

for(p in unique_primirna){
  brackets_5p<-mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), 'BRACKETS']
  brackets_3p<-mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), 'BRACKETS']
  diff<-brackets_5p-brackets_3p
  if(brackets_5p>brackets_3p){
    struct_5p<-mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), 'MOR.EDITED']
    bracket_pos<-gregexpr('(', toString(struct_5p), fixed=T)[[1]]
    
    for(i in 1:diff){
      curr_pos<-(struct_len-bracket_pos[i]+1)*-1
      mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), toString(curr_pos)]<-'*'
    }
    
    pair_col_start<-grep('-30', colnames(mor_struct_5p))
    new_struct_5p<-paste(mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), pair_col_start:ncol(mor_struct_5p)], collapse='')
    mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), 'MOR.EDITED']<-new_struct_5p
    mor_struct_5p[which(mor_struct_5p$PRIMIRNA==p), 'BRACKETS']<-str_count(toString(new_struct_5p), '\\(')
  }else if(brackets_3p>brackets_5p){
    diff<-diff*-1
    print(p)
    print(diff)
    struct_3p<-mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), 'MOR.EDITED']
    print(struct_3p)
    bracket_pos<-gregexpr(')', toString(struct_3p), fixed=T)[[1]]
    for(i in seq(from=0, to=diff-1, by=1)){
      curr_pos<-bracket_pos[length(bracket_pos)-i]
      mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), toString(curr_pos)]<-'*'
    }
    
    pair_col_start<-grep('1', colnames(mor_struct_3p))[1]
    new_struct_3p<-paste(mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), pair_col_start:ncol(mor_struct_3p)], collapse='')
    mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), 'MOR.EDITED']<-new_struct_3p
    mor_struct_3p[which(mor_struct_3p$PRIMIRNA==p), 'BRACKETS']<-str_count(toString(new_struct_3p), '\\)')
    print(new_struct_3p)
  }
}


mor_struct_file_3p<-paste0(base_path, 'low-fidelity-moR-struct-3p-even-brackets.tsv')
write.table(mor_struct_3p, file=mor_struct_file_3p, sep='\t', row.names=F)
mor_struct_3p[,'0']<-NULL
base_pair_ratio_3p<-get3PBasePairRatioPerPosition(mor_struct_3p)

file<-paste(base_path, 'low-fidelity-moR-base-pairing-3p-even-brackets.tsv', sep='')

write.table(base_pair_ratio_3p, file=file, sep='\t', row.names=F)


mor_struct_file_5p<-paste0(base_path, 'low-fidelity-moR-struct-5p-even-brackets.tsv')
write.table(mor_struct_5p, file=mor_struct_file_5p, sep='\t', row.names=F)
mor_struct_5p[,'0']<-NULL
base_pair_ratio_5p<-get5PBasePairRatioPerPosition(mor_struct_5p)

file<-paste(base_path, 'low-fidelity-moR-base-pairing-5p-even-brackets.tsv', sep='')
#
write.table(base_pair_ratio_5p, file=file, sep='\t', row.names=F)
