
## PURPOSE: edit TCGA sequences to reflect only templated bases
## INPUT: miRBase data                   miRBase21-master.tsv
##        isomiR expression data  tumor.isomir.expression.tsv

mapTemplatedBases<-function(mirbase_data, input_data){
  
  mapped_bases<-input_data[,c('MIRNA', 'SEQUENCE')]

  for(i in 1:nrow(mapped_bases)){
    
    mirna<-toString(mapped_bases[i, 'MIRNA'])
    consensus<-toString(unique(mirbase_data[which(mirbase_data$MIRNA==mirna), 'SEQUENCE']))
    motif<-toString(unique(mirbase_data[which(mirbase_data$MIRNA==mirna), 'MOTIF.13']))
    sequence<-toString(mapped_bases[i, 'SEQUENCE'])
    pri_sequence<-mirbase_data[which(mirbase_data$MIRNA==mirna), 'PRI.SEQUENCE'] #primary sequences of miRNA
    
    curr_mapped_base<-strsplit(sequence, '')[[1]]
    #get start of consensus
    motif_start<-gregexpr(motif, consensus, fixed=TRUE)[[1]][1]
    if(motif_start==-1){ 
      print(sequence) 
      next 
    }
    consensus_start<-gregexpr(motif, sequence, fixed=TRUE)[[1]]-motif_start
    mapped_bases[i,'CONSENSUS.START']<-consensus_start
    
    
    
    threshold<-3
    mapped_5p<-c()
    if(consensus_start>0){
      mor_bases<-substr(sequence, 1, consensus_start)
      mapped_5p<-strsplit(mor_bases, '')[[1]]
      if(nchar(mor_bases)>=threshold){
        for(x in 1:(nchar(mor_bases)-threshold+1)){
          temp_mor_bases<-substr(mor_bases, nchar(mor_bases)-threshold-x+2, nchar(mor_bases)-x+1)
          templated<-FALSE
          for(p in pri_sequence){
            pri_templated_start<-gregexpr(consensus, toString(p), fixed=TRUE)[[1]][1]
            pri_base<-substr(toString(p), pri_templated_start-threshold-x+1, pri_templated_start-x)
            if(gregexpr(temp_mor_bases, pri_base, fixed=TRUE)[[1]][1]==1 || gregexpr(pri_base, temp_mor_bases, fixed=TRUE)[[1]][1]==1){
              templated<-TRUE
              break
            }
          }

          if(templated){
            mapped_5p[(length(mapped_5p)-nchar(temp_mor_bases)-x+2):(length(mapped_5p)-x+1)]<-'.'
          }else{
            checked_seq<-mapped_5p[length(mapped_5p)-nchar(temp_mor_bases)-x+2:length(mapped_5p)-x+1]
            if(!'.'%in%checked_seq){
              mapped_5p[(nchar(mor_bases)-nchar(temp_mor_bases)-x+2):(nchar(mor_bases)-x+1)]<-strsplit(temp_mor_bases, '')[[1]]
            }
          }
        }
      }
    }
    
    if(length(mapped_5p)>0){
      for(x in 1:length(mapped_5p)){
        curr_base<-mapped_5p[x]
        mapped_bases[i, x+2]<-curr_base
        colnames(mapped_bases)[x+2]<-toString((length(mapped_5p)*-1)+(x-1))
      } 
    }

    sequence<-substr(sequence, consensus_start+1, nchar(sequence))
    mapped_consensus<-c()
  
    for(x in 1:nchar(sequence)){
      consensus_base<-substr(consensus, x, x)
      sequence_base<-substr(sequence, x, x)
      if(consensus_base%in%sequence_base){
        mapped_consensus[x]<-'.'
      }else{
        mapped_3p<-strsplit(substr(sequence, x, nchar(sequence)), '')[[1]]
        consensus_end<-x
        break
      }
    }
    
    templated_seq<-paste(mapped_consensus, collapse='')
    #threshold for marking templated nucleotides
    threshold<-3
    if(length(mapped_3p)>=threshold){

      templated_seq<-toString(substr(sequence, 1, length(curr_mapped_base))) #templated bases found so far in sequence
      for(x in 1:(length(mapped_3p)-threshold+1)){
        temp_mapped_3p<-paste(mapped_3p[x:(x+threshold-1)], collapse='')
        templated<-FALSE
        for(p in pri_sequence){
          pri_templated_end<-gregexpr(templated_seq, toString(p), fixed=TRUE)[[1]][1]+nchar(templated_seq) #end of templated bases in primary sequence
          pri_base<-substr(toString(p), x+pri_templated_end, x+pri_templated_end+threshold-1)
          if(gregexpr(temp_mapped_3p, pri_base, fixed=TRUE)[[1]][1]==1 || gregexpr(pri_base, temp_mapped_3p, fixed=TRUE)[[1]][1]==1){
            templated<-TRUE
            break
          }
        }
        if(templated){
          mapped_3p[x:(x+nchar(temp_mapped_3p)-1)]<-'.'
        }else{
          checked_seq<-mapped_3p[x:(x+threshold-1)]
          if(!'.'%in%checked_seq){
            mapped_3p[x:(x+nchar(temp_mapped_3p)-1)]<-strsplit(temp_mapped_3p, '')[[1]]
          }
        }
      }
    }
    mapped_consensus<-c(mapped_consensus, mapped_3p)
    #curr_mapped_base<-c(curr_mapped_base, strsplit(substr(sequence, length(curr_mapped_base)+1, nchar(sequence)), '')[[1]])

    if(length(mapped_consensus)<nchar(consensus)){
      mapped_consensus[(length(mapped_consensus)+1):nchar(consensus)]<-'*'
    }

    #add mapped bases to dataframe
    for(x in 1:length(mapped_consensus)){
      curr_base<-mapped_consensus[x]
      mapped_bases[i, x+2+length(mapped_5p)]<-curr_base
      colnames(mapped_bases)[x+2+length(mapped_5p)]<-toString(x)
    }
  }

  mapped_bases[is.na(mapped_bases)]<-''
  mapped_bases$MIRNA<-NULL
  mapped_bases$PRI.SEQUENCE<-NULL
  mapped_bases<-merge(input_data, mapped_bases, by='SEQUENCE', all=TRUE, na.rm=TRUE)
  return(unique(mapped_bases))
}


