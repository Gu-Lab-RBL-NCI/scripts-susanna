
## PURPOSE: determine if precursor 2d struct is found in primiRNA 2d struct (generated by RNAfold)
## INPUT: miRBase data      miRBase21-master.tsv
##        isomiR summary data   tumor.isomir.tsv
## OUTPUT: table containing fidelity & RPM of mature strands and if pre struct in pri struct    precursor-vs-primary-2d-struct-miRBase.tsv


require('stringr')
require('plyr')
require('stringdist')

base_path<-'/Users/chens22/Documents/miRNA/'
mirna_file<- paste0(base_path, 'miRBase21-master.tsv')
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors=F)
summary_ext<-'.isomir.tsv'

tumors <- c('ACC')


getPriMirna<-function(mirna){
  return(unique(mirna_data[which(mirna_data$MIRNA==mirna | mirna_data$PRIMIRNA==mirna), 'PRIMIRNA']))
}

get3PMirna<-function(primirna){
  pri_mirna<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(primirna)), 'PRIMIRNA'])
  mirna_3p<-toString(unique(mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='3P'), 'MIRNA']))
  if(nchar(mirna_3p)>0){
    return(toString(mirna_3p))
  }else{
    print(paste0('no 3p mirna: ', primirna))
    return('')
  }
}

get5PMirna<-function(primirna){
  pri_mirna<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(primirna)), 'PRIMIRNA'])
  mirna_5p<-toString(unique(mirna_data[which(mirna_data$PRIMIRNA==primirna & mirna_data$STRAND=='5P'), 'MIRNA']))
  if(nchar(mirna_5p)>0){
    return(toString(mirna_5p))
  }else{
    print(paste0('no 5p mirna: ', primirna))
    return('')
  }
}

getExtendedSequence<-function(primirna){
  extended_seq<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(primirna)), 'EXTENDED.SEQUENCE'])
  if(length(extended_seq)>1){
    print(paste('ERROR EXTENDED SEQ', mirna))
    return(NA)
  }
  pre_seq<-getPrecursorSequence(primirna)
  if(nchar(pre_seq)<1){
    return('')
  }
  pre_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  pre_end<-pre_start+nchar(pre_seq)-1+struct_len
  pre_start<-pre_start-struct_len
  extended_seq<-substr(extended_seq, pre_start, pre_end)
  return(toString(extended_seq))
}

getPrimarySequence<-function(mirna){
  pri_seqs<-unique(mirna_data[which(mirna_data$PRIMIRNA==getPriMirna(mirna)), 'PRI.SEQUENCE'])
  return(toString(pri_seqs))
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

  if(length(start)<1 || length(end)<1){
    return('')
  }
  pre_seq<-substr(pri_seq, start, end)
  return(toString(pre_seq))
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
  input_file<-paste(base_path, '/structure/temp-miRNA-seq.fa', sep='')
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-paste(base_path, '/structure/temp-miRNA-secondary-struct.txt', sep='')
  secondary_structure<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  secondary_structure<-strsplit(secondary_structure, ' ', fixed=T)[[1]][1]
  return(toString(secondary_structure))
}

getEnergy<-function(sequence, mirna){
  input_data<-c(paste('>', mirna, sep=''), sequence)
  input_file<-paste(base_path, '/structure/temp-miRNA-seq.fa', sep='')
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-paste(base_path, '/structure/temp-miRNA-secondary-struct.txt', sep='')
  output_data<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  energy<-as.numeric(strsplit(output_data, ' ', fixed=T)[[1]][4])
  energy<-gsub(')', '', toString(energy), fixed=T)
  energy<-gsub('(', '', toString(energy), fixed=T)
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][2]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][3]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  return(toString(energy))
}

getExtendedSecondaryStructure<-function(primirna, struct_len){
  extended_seq<-getExtendedSequence(primirna)
  extended_struct<-getSecondaryStructure(extended_seq, primirna)
  return(toString(extended_struct))
}

getPrecursorSecondaryStructure<-function(primirna){
  # pri_seq<-getPrimarySequence(primirna)
  # pre_seq<-getPrecursorSequence(primirna)
  # if(nchar(pre_seq)<1){
  #   return('')
  # }
  # pre_start<-gregexpr(toString(pre_seq), toString(pri_seq), fixed=T)[[1]]
  # pre_end<-pre_start+nchar(pre_seq)-1+13
  # pre_start<-pre_start-11
  # pri_seq<-substr(pri_seq, pre_start, pre_end)
  # pri_secondary_struct<-getSecondaryStructure(pri_seq, primirna)
  # 
  # seq_5p<-get5PSequence(primirna)
  # seq_3p<-get3PSequence(primirna)
  # if(nchar(seq_5p)>0){
  #   start<-gregexpr(seq_5p, pri_seq)[[1]] 
  # }else{
  #   end<-gregexpr(seq_3p, pri_seq)[[1]]+nchar(seq_3p)-1 
  #   temp_pri_struct<-substr(pri_secondary_struct, 1, end)
  #   precurs_brackets_3p<-str_count(temp_pri_struct, '\\)')
  #   mor_5p_brackets<-str_count(pri_secondary_struct, '\\(')
  #   precursor_start<-mor_5p_brackets-precurs_brackets_3p+1
  #   start<-gregexpr('(', pri_secondary_struct, fixed=T)[[1]][precursor_start]
  # }
  # if(nchar(seq_3p)>0){
  #   end<-gregexpr(seq_3p, pri_seq)[[1]]+nchar(seq_3p)-1 
  # }else{
  #   start<-gregexpr(seq_5p, pri_seq)[[1]]
  #   temp_pri_struct<-substr(pri_secondary_struct, start, nchar(pri_secondary_struct))
  #   precurs_brackets_5p<-str_count(temp_pri_struct, '\\(')
  #   end<-gregexpr(')', pri_secondary_struct, fixed=T)[[1]][precurs_brackets_5p]
  # }
  # 
  # if(length(start)==0 || length(end)==0){
  #   return('')
  # }
  # pre_secondary_struct<-substr(pri_secondary_struct, start, end)
  # 
  pre_seq<-getPrecursorSequence(primirna)
  pre_secondary_struct<-getSecondaryStructure(pre_seq, primirna)
  
  #take out unpaired bases at the end
  start<-gregexpr('(', pre_secondary_struct, fixed=T)[[1]][1]
  end<-gregexpr(')', pre_secondary_struct, fixed=T)[[1]]
  end<-end[length(end)]
  pre_secondary_struct<-substr(pre_secondary_struct, start, end)
  return(toString(pre_secondary_struct))
}

getExtendedPrecursorSecondaryStructure<-function(primirna, struct_len){
  pre_seq<-getPrecursorSequence(primirna)
  pre_secondary_struct<-getSecondaryStructure(getPrecursorSequence(primirna), primirna)
  
  #take out unpaired bases at the end
  start<-gregexpr('(', pre_secondary_struct, fixed=T)[[1]][1]
  end<-gregexpr(')', pre_secondary_struct, fixed=T)[[1]]
  end<-end[length(end)]
  pre_seq<-substr(pre_seq, start, end)
  
  #get pre struct in extended struct
  extended_seq<-getExtendedSequence(primirna)
  pre_start<-gregexpr(pre_seq, extended_seq, fixed=T)[[1]]
  pre_end<-pre_start+nchar(pre_seq)-1
  extended_struct<-getExtendedSecondaryStructure(primirna, struct_len)
  extended_pre_struct<-substr(extended_struct, pre_start, pre_end)
  return(extended_pre_struct)
}

getFidelity3P<-function(primirna, tumor){
  mirna<-get3PMirna(primirna)
  if(nchar(mirna)<1){
    return(NA)
  }
  fid_data<-summary_data[which(summary_data$MIRNA==mirna,
                                  summary_data$TUMOR==tumor),]
  if(nrow(fid_data)<1){
    return(NA)
  }
  fidelity<-mean(fid_data$FIDELITY.5P, na.rm=T)
  return(fidelity)
}

getFidelity5P<-function(primirna, tumor){
  mirna<-get5PMirna(primirna)
  if(nchar(mirna)<1){
    return(NA)
  }
  fid_data<-summary_data[which(summary_data$MIRNA==mirna,
                                  summary_data$TUMOR==tumor),]
  if(nrow(fid_data)<1){
    return(NA)
  }
  fidelity<-mean(fid_data$FIDELITY.5P, na.rm=T)
  return(fidelity)
}

get5PRPM<-function(primirna, tumor){
  mirna<-get5PMirna(primirna)
  if(nchar(mirna)<1){return(NA)}
  fid_data<-summary_data[which(summary_data$MIRNA==mirna,
                                   summary_data$TUMOR==tumor),]
  if(nrow(fid_data)<1){
    return(NA)
  }
  fid_data[,'RPM']<-(fid_data$TOTAL.READS/fid_data$TOTAL.READS.IN.SAMPLE)*1000000
  return(mean(fid_data$RPM, na.rm=T))
}

get3PRPM<-function(primirna, tumor){
  mirna<-get3PMirna(primirna)
  if(nchar(mirna)<1){return(NA)}
  fid_data<-summary_data[which(summary_data$MIRNA==mirna,
                                  summary_data$TUMOR==tumor),]
  if(nrow(fid_data)<1){
    return(NA)
  }
  fid_data[,'RPM']<-(fid_data$TOTAL.READS/fid_data$TOTAL.READS.IN.SAMPLE)*1000000
  return(mean(fid_data$RPM, na.rm=T))
}



for(t in tumors){

  summary_file<-paste0('/Users/chens22/Documents/miRNA/', t, '/summary_files/', t, summary_ext)
  summary_data<-read.table(file=summary_file, sep='\t', header=T)
  summary_data[,'TUMOR']<-t


  unique_primirna<-unique(mirna_data$PRIMIRNA)
  unique_primirna<-unique_primirna[!is.na(unique_primirna)]
  pre_in_pri_data<-data.frame(PRIMIRNA=NA,
                               PRE.STRUCT.ENERGY=NA,
                               PRECURSOR.STRUCT=NA,
                               PRIMARY.PRECURSOR.STRUCT=NA,
                               PRE.DIST=NA,
                               PRIMARY.STRUCT=NA,
                               PRE.IN.PRI=NA,
                               FIDELITY.5P=NA,
                               RPM.5P=NA,
                               FIDELITY.3P=NA,
                               RPM.3P=NA,
                               stringsAsFactors = F)

  for(p in unique_primirna){
    pre_struct<-toString(getPrecursorSecondaryStructure(p))
    pre_struct_energy<-getEnergy(getPrecursorSequence(p), p)
    pri_struct<-toString(getExtendedSecondaryStructure(p, 30))
    #pre_in_pri<-F
    pre_in_pri<-gregexpr(pre_struct, pri_struct, fixed=T)[[1]]!=-1
    pri_pre_struct<-getExtendedPrecursorSecondaryStructure(p, 30)
    pre_struct_dist<-stringdist(toString(pre_struct), toString(pri_pre_struct))
    
    fidelity_5p<-getFidelity5P(p, t)
    rpm_5p<-get5PRPM(p, t)

    fidelity_3p<-getFidelity3P(p, t)
    rpm_3p<-get3PRPM(p, t)
    
    pre_in_pri_data<-rbind.fill(pre_in_pri_data, data.frame(PRIMIRNA=toString(p),
                                                            PRE.STRUCT.ENERGY=toString(pre_struct_energy),
                                                            PRECURSOR.STRUCT=pre_struct,
                                                            PRIMARY.PRECURSOR.STRUCT=pri_pre_struct,
                                                            PRE.DIST=pre_struct_dist,
                                                            PRIMARY.STRUCT=pri_struct,
                                                            PRE.IN.PRI=pre_in_pri, 
                                                            FIDELITY.5P=fidelity_5p,
                                                            RPM.5P=rpm_5p,
                                                            FIDELITY.3P=fidelity_3p,
                                                            RPM.3P=rpm_3p,
                                                            stringsAsFactors = F))
  }

  pre_in_pri_data<-pre_in_pri_data[-1,]
  pre_in_pri_data<-unique(pre_in_pri_data)
  file<-paste0(base_path, '/structure/precursor-vs-primary-2d-struct-miRBase.tsv')
  write.table(pre_in_pri_data, file=file, sep='\t', row.names=F)


}
