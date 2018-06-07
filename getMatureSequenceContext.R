
## PURPOSE: get IUPAC sequence of moR and loop of miRNA sequence
## INPUT: miRBase data  miRBAse21-master.tsv
## OUTPUT: context sequences of miRNA   miRBase21-sequence-context.tsv

getContextBase<-function(bases){
  if(length(bases)==2){
    if(all(c('A', 'G')%in%bases)){ return('W') }
    if(all(c('A', 'T')%in%bases)){ return('W') }
    if(all(c('C', 'T')%in%bases)){ return('Y') }
    if(all(c('C', 'G')%in%bases)){ return('S') }
    if(all(c('T', 'G')%in%bases)){ return('K') }
    if(all(c('A', 'C')%in%bases)){ return('M') }
  }
  if(length(bases)==3){
    if(all(c('A', 'T', 'C')%in%bases)){ return('H') }
    if(all(c('A', 'T', 'G')%in%bases)){ return('D') }
    if(all(c('C', 'T', 'G')%in%bases)){ return('B') }
    if(all(c('A', 'G', 'C')%in%bases)){ return('V') }
  }
  return(bases[1])
}

mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-master.tsv'


mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE)

para_mirna<-unique(mirna_data[which(mirna_data$PARALOGS>0), c('MIRNA', 'PARALOGS','STRAND', 'SEQUENCE', 'PRI.SEQUENCE')])

data<-cbind(para_mirna, data.frame(CONTEXT.SEQUENCE=NA))

for(i in 1:nrow(data)){
  mirna<-toString(data[i, 'MIRNA'])
  sequence<-toString(data[i, 'SEQUENCE'])
  strand<-toString(data[i, 'STRAND'])
  curr_pri_sequence<-toString(data[i, 'PRI.SEQUENCE'])
  pri_mirna<-unique(mirna_data[which(mirna_data$MIRNA==mirna), 'PRIMIRNA'])
  pri_data<-unique(mirna_data[which(mirna_data$PRIMIRNA%in%pri_mirna), c('EXTENDED.SEQUENCE', 'STRAND', 'SEQUENCE')])
  pri_seq<-unique(pri_data$EXTENDED.SEQUENCE)
  
  if(strand%in%'5P'){
    
    #get minimum 5p mor length
    # min_mor_len<-Inf
    # for(p in pri_seq){
    #   start_5p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
    #   mor_len<-start_5p-1
    #   if(start_5p<min_mor_len){ min_mor_len<-mor_len }
    # }
    min_mor_len<-5
    
    #get 5p moR sequence
    pri_mor<-c()
    for(p in pri_seq){
      start_5p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
      if(length(start_5p)>1){
        start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
        start_5p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
      }
      curr_pri_5p<-substr(p, start_5p-min_mor_len, start_5p-1)
      pri_mor<-c(pri_mor, curr_pri_5p)
    }
    if(nchar(pri_mor[1])>0){
      data[i, 'MOR']<-paste(pri_mor, collapse=',')
    }
    
    #get minimum loop length
    min_loop_len<-Inf
    pri_3p<-pri_data[which(pri_data$STRAND=='3P'),]
    if(nrow(pri_3p)>0){
      for(x in 1:nrow(pri_3p)){
        p<-pri_3p[x, 'EXTENDED.SEQUENCE']
        seq_3p<-pri_3p[x, 'SEQUENCE']
        start_5p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
        if(length(start_5p)>1){
          start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          start_5p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
        }
        end_5p<-start_5p+nchar(sequence)
        start_3p<-gregexpr(seq_3p, p, fixed=TRUE)[[1]]
        if(length(start_3p)>1){
          start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          start_3p<-start_pri+gregexpr(seq_3p, curr_pri_sequence, fixed=TRUE)[[1]]
        }
        loop_len<-start_3p-end_5p-1
        if(loop_len<min_loop_len){ min_loop_len<-loop_len }
      }
    }else{
      # for(p in pri_seq){
        # start_5p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
        # if(length(start_5p)>1){
          # start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          # start_5p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
        # }
        # end_5p<-start_5p+nchar(sequence)
        # loop_len<-nchar(p)-end_5p
        # if(loop_len<min_loop_len){ min_loop_len<-loop_len }
      # }
      min_loop_len<-10
    }
    
    #get loop sequence
    pri_loop<-c()
    for(p in pri_seq){
      start_5p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
      if(length(start_5p)>1){
        start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
        start_5p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
      }
      end_5p<-start_5p+nchar(sequence)
      curr_pri_loop<-substr(p, end_5p, end_5p+min_loop_len)
      pri_loop<-c(pri_loop, curr_pri_loop)
    }
    if(nchar(pri_loop[1])>0){
      data[i, 'LOOP']<-paste(pri_loop, collapse=',')
    }

  }
  
  if(strand%in%'3P'){
    
    #get minimum 3p mor length
    # min_mor_len<-Inf
    # for(p in pri_seq){
    #   start_3p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
    #   end_3p<-start_3p+nchar(sequence)
    #   mor_len<-nchar(p)-end_3p
    #   if(mor_len<min_mor_len){ min_mor_len<-mor_len }
    # }
    min_mor_len<-5
    
    #get 3p moR sequence
    pri_mor<-c()
    for(p in pri_seq){
      start_3p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
      if(length(start_3p)>1){
        start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
        start_5p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
      }
      end_3p<-start_3p+nchar(sequence)
      curr_pri_3p<-substr(p, end_3p, end_3p+min_mor_len)
      pri_mor<-c(pri_mor, curr_pri_3p)
    }
    if(nchar(pri_mor[1])>0){
      data[i, 'MOR']<-paste(pri_mor, collapse=',')
    }
    
    #get minimum loop length
    min_loop_len<-Inf
    pri_5p<-pri_data[which(pri_data$STRAND=='5P'),]
    if(nrow(pri_5p)>0){
      for(x in 1:nrow(pri_5p)){
        p<-pri_5p[x, 'EXTENDED.SEQUENCE']
        seq_5p<-pri_5p[x, 'SEQUENCE']
        start_5p<-gregexpr(seq_5p, p, fixed=TRUE)[[1]]
        if(length(start_5p)>1){
          start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          start_5p<-start_pri+gregexpr(seq_5p, curr_pri_sequence, fixed=TRUE)[[1]]
        }
        end_5p<-start_5p+nchar(sequence)
        start_3p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
        if(length(start_3p)>1){
          start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          start_3p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
        }
        loop_len<-start_3p-end_5p-1
        if(loop_len<min_loop_len){ min_loop_len<-loop_len }
      }
    }else{
      # for(p in pri_seq){
        # start_3p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
        # if(length(start_3p)>1){
          # start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
          # start_3p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
        # }
        # loop_len<-start_3p-1
        # if(loop_len<min_loop_len){ min_loop_len<-loop_len }
      # }
      min_loop_len<-10
    }
    
    #get loop sequence
    pri_loop<-c()
    for(p in pri_seq){
      start_3p<-gregexpr(sequence, p, fixed=TRUE)[[1]]
      if(length(start_3p)>1){
        start_pri<-gregexpr(curr_pri_sequence, p, fixed=TRUE)[[1]]
        start_3p<-start_pri+gregexpr(sequence, curr_pri_sequence, fixed=TRUE)[[1]]
      }
      curr_pri_loop<-substr(p, start_3p-min_loop_len, start_3p-1)
      pri_loop<-c(pri_loop, curr_pri_loop)
    }
    if(nchar(pri_loop[1])>0){
      data[i, 'LOOP']<-paste(pri_loop, collapse=',')
    }
    
  }
  if(length(min_loop_len)>1){
    print(mirna)
  }
  
  #get context sequence
  context_loop<-c()
  for(x in 1:min_loop_len){
    bases<-c()
    for(l in pri_loop){
      curr_base<-substr(l, x, x)
      bases<-c(bases, curr_base)
    }
    bases<-unique(bases)
    context_base<-getContextBase(bases)
    context_loop<-c(context_loop, context_base)
  }
  context_loop<-paste(context_loop, collapse='')
  
  #get context sequence
  context_mor<-c()
  for(x in 1:nchar(pri_mor[1])){
    bases<-c()
    for(m in pri_mor){
      curr_base<-substr(m, x, x)
      bases<-c(bases, curr_base)
    }
    bases<-unique(bases)
    context_base<-getContextBase(bases)
    context_mor<-c(context_mor, context_base)
  }
  context_mor<-paste(context_mor, collapse='')
  
  if(strand=='3P'){
    context_seq<-paste(context_loop, sequence, context_mor, sep='')
  }else if(strand=='5P'){
    context_seq<-paste(context_mor, sequence, context_loop, sep='')
  }
  
  data[i, 'CONTEXT.SEQUENCE']<-context_seq
}

data$PRI.SEQUENCE<-NULL
data<-unique(data)

result_file<-'/Users/chens22/Documents/miRNA/miRBase21/miRBase21-sequence-context.tsv'
write.table(data, file=result_file, sep='\t', row.names=FALSE)
