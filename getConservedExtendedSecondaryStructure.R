
## PURPOSE: get secondary structure of primiRNA sequence with the most moR pairing
## INPUT: miRBase data  miRBase21-master.tsv
## OUTPUT: 2d struct with max moR pairing   extended-2d-struct-max-moR-pairing-all-miRBase.tsv

require('plyr')
require('stringr')

base_path<-'/Users/chens22/Documents/miRNA/structure/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors=F)
conserved_struct_file<-paste0(base_path, 'extended-2d-struct-max-moR-pairing-all-miRBase.tsv')

low_fid_mirna<-c('hsa-miR-758-5p', 'hsa-miR-153-5p','hsa-miR-383-5p','hsa-miR-382-5p','hsa-miR-183-5p',
                 'hsa-miR-425-5p','hsa-miR-136-5p','hsa-miR-628-5p','hsa-miR-769-5p','hsa-miR-30e-5p',
                 'hsa-miR-98-5p','hsa-miR-138-5p-1-2','hsa-miR-144-5p','hsa-let-7i-5p', 'hsa-miR-191-5p',
                 'hsa-miR-532-5p')
mid_fid_mirna<-c('hsa-miR-139-5p','hsa-miR-432-5p', 'hsa-miR-9-5p-1-2-3','hsa-let-7g-5p',
                 'hsa-miR-218-5p-1-2','hsa-miR-186-5p','hsa-miR-93-5p','hsa-miR-99a-5p','hsa-let-7d-5p',
                 'hsa-miR-16-5p-1-2','hsa-miR-744-5p','hsa-miR-182-5p','hsa-miR-584-5p',
                 'hsa-miR-181a-5p-1-2','hsa-miR-20a-5p','hsa-miR-885-5p','hsa-miR-15a-5p','hsa-miR-409-5p')
high_fid_mirna<-c('hsa-miR-100-5p','hsa-miR-106b-5p','hsa-miR-129-5p-1-2','hsa-miR-589-5p','hsa-miR-125b-5p-1-2',
                  'hsa-miR-30a-5p','hsa-miR-125a-5p','hsa-miR-134-5p','hsa-miR-99b-5p','hsa-miR-330-5p',
                  'hsa-miR-26a-5p-1-2','hsa-miR-21-5p','hsa-let-7b-5p','hsa-miR-1307-5p','hsa-miR-185-5p',
                  'hsa-miR-361-5p','hsa-miR-379-5p','hsa-let-7f-5p-1-2','hsa-let-7e-5p','hsa-miR-340-5p',
                  'hsa-miR-149-5p','hsa-miR-204-5p','hsa-miR-26b-5p','hsa-miR-30d-5p')


getPriMirna<-function(mirna){
  return(unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA']))
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
  return(pre_seq)
}

getExtendedSequence<-function(mirna){
  extended_seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(mirna)), 'EXTENDED.SEQUENCE'])
  if(length(extended_seq)>1){
    print(paste('ERROR EXTENDED SEQ', mirna))
    return(NA)
  }
  return(toString(extended_seq))
}

getSecondaryStructure<-function(sequence, mirna){
  input_data<-c(paste('>', mirna, sep=''), sequence)
  input_file<-paste(base_path, 'temp-miRNA-conserved-seq.fa', sep='')
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-paste(base_path, 'temp-miRNA-secondary-struct.txt', sep='')
  secondary_structure<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  secondary_structure<-strsplit(secondary_structure, ' ', fixed=T)[[1]][1]
  return(secondary_structure)
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

unique_mirna<-unique(mirna_data$MIRNA)
conserved_structs<-read.table(file=conserved_struct_file, sep='\t', header=T, stringsAsFactors = F)
conserved_structs<-conserved_structs[!(conserved_structs$END>30&!is.na(conserved_structs$END)),]
count<-0
#unique_mirna<-c('hsa-mir-378d-1')
for(m in unique_mirna){
  count<-count+1
  if(count%%100==0){print(count/length(unique_mirna))}
  primirna<-getPriMirna(m)
  for(p in primirna){
    if(!p%in%conserved_structs$PRIMIRNA){
      pre_struct<-getPrecursorSecondaryStructure(p)
      if(is.na(pre_struct)||nchar(pre_struct)<1){
        print(paste('cannot retrieve precursor seq: ', p))
        conserved_structs<-rbind(conserved_structs, data.frame(PRIMIRNA=p,
                                                               PRECURSOR.STRUCT=NA,
                                                               SECONDARY.STRUCT=NA,
                                                               PRE.IN.PRI=F,
                                                               START=NA,
                                                               END=NA))
        next
      }
      extended_seq<-getExtendedSequence(p)
      pre_seq<-getPrecursorSequence(p)
      if(nchar(pre_seq)>0){
        start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
        end<-start+nchar(pre_seq)-1
        start<-start-11
        end<-end+13
        curr_seq<-substr(extended_seq, start, end)
        curr_struct<-getSecondaryStructure(curr_seq, p)
        added_start<-gregexpr(pre_struct, curr_struct, fixed=T)[[1]]-1
        added_end<-nchar(curr_struct)-(added_start+nchar(pre_struct))
        max_brackets_5p<-get5PMorBrackets(pre_struct, curr_struct)
        max_brackets_3p<-get3PMorBrackets(pre_struct, curr_struct)
        if(gregexpr(toString(pre_struct), toString(curr_struct), fixed=T)[[1]][1]!=-1){
          while(gregexpr(pre_struct, getSecondaryStructure(substr(extended_seq, start, end), p), fixed=T)[[1]][1]!=-1 &
                ((added_start<30 & added_end<=30)) &
                (get5PMorBrackets(pre_struct, getSecondaryStructure(substr(extended_seq, start, end), p))>=max_brackets_5p) &
                (get3PMorBrackets(pre_struct, getSecondaryStructure(substr(extended_seq, start, end), p))>=max_brackets_3p)
          ){
            curr_seq<-substr(extended_seq, start, end)
            curr_struct<-getSecondaryStructure(curr_seq, p)
            added_start<-gregexpr(pre_struct, curr_struct, fixed=T)[[1]]-1
            added_end<-nchar(curr_struct)-(added_start+nchar(pre_struct))
            max_brackets_5p<-get5PMorBrackets(pre_struct, curr_struct)
            max_brackets_3p<-get3PMorBrackets(pre_struct, curr_struct)
            if(nchar(curr_struct)%%2==0){
              if(added_start<30){
                start<-start-1
              }
              if(gregexpr(pre_struct, getSecondaryStructure(substr(extended_seq, start, end), p), fixed=T)[[1]][1]==-1 & added_end<30){
                start<-start+1
                end<-end+1
              }
            }else if(added_start<30 & added_end==30){
              start<-start-1
            }else{
              if(added_end<30){
                end<-end+1
              }
              if(gregexpr(pre_struct, getSecondaryStructure(substr(extended_seq, start, end), p), fixed=T)[[1]][1]==-1 & added_start<30){
                end<-end-1
                start<-start-1
              }
            }
            added_start<-gregexpr(pre_struct, curr_struct, fixed=T)[[1]]-1
            added_end<-nchar(curr_struct)-(added_start+nchar(pre_struct))
            # if(added_start>30){
            #   start<-start+1
            #   curr_seq<-substr(extended_seq, start, end)
            #   curr_struct<-getSecondaryStructure(curr_seq, p)
            # }
            # if(added_end>30){
            #   end<-end-1
            #   curr_seq<-substr(extended_seq, start, end)
            #   curr_struct<-getSecondaryStructure(curr_seq, p)
            # }
          }
          start<-gregexpr(pre_struct, curr_struct, fixed=T)[[1]]
          if(start==-1){pre_in_pri<-F}
          else{pre_in_pri<-T}
          end<-nchar(curr_struct)-(start+nchar(pre_struct)-1)
          start<-(start-1)*-1
          conserved_structs<-rbind(conserved_structs, data.frame(PRIMIRNA=p,
                                                                 PRECURSOR.STRUCT=pre_struct,
                                                                 SECONDARY.STRUCT=curr_struct,
                                                                 PRE.IN.PRI=pre_in_pri,
                                                                 START=start,
                                                                 END=end))
        }else{
          conserved_structs<-rbind(conserved_structs, data.frame(PRIMIRNA=p,
                                                                 PRECURSOR.STRUCT=pre_struct,
                                                                 SECONDARY.STRUCT=curr_struct,
                                                                 PRE.IN.PRI=FALSE,
                                                                 START=NA,
                                                                 END=NA))
        }
      
      }else{
        conserved_structs<-rbind(conserved_structs, data.frame(PRIMIRNA=p,
                                                               PRECURSOR.STRUCT=NA,
                                                               SECONDARY.STRUCT=NA,
                                                               PRE.IN.PRI=FALSE,
                                                               START=NA,
                                                               END=NA))
        write.table(conserved_structs, file=paste0(base_path, 'extended-2d-struct-max-moR-pairing-all-miRBase.tsv'), sep='\t', row.names=F)
      }
    }
    write.table(conserved_structs, file=paste0(base_path, 'extended-2d-struct-max-moR-pairing-all-miRBase.tsv'), sep='\t', row.names=F)
  }
  conserved_structs<-read.table(file=conserved_struct_file, sep='\t', header=T)
}
conserved_structs<-unique(conserved_structs)
#conserved_structs<-conserved_structs[-1,]
write.table(conserved_structs, file=paste0(base_path, 'extended-2d-struct-max-moR-pairing-all-miRBase.tsv'), sep='\t', row.names=F)
