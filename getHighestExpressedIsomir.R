
## PURPOSE: get highest expressed isomiR across specified tumors
## INPUT: miRBase data  miRBase21-master.tsv
##        isomiR expression data  tumor.isomir.expression.tsv
## OUTPUT: table containing highest expressed isomiR from TCGA data   most-expressed-TCGA-consensus.tsv

require('stringdist')
require('plyr')


base_path<-'/Users/chens22/Documents/miRNA/'
mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21-master.tsv'
edited_seq_file<-'/Users/chens22/Documents/miRNA/TCGA-sequence.tsv'
summary_ext<-'.isomir.tsv'
expression_ext<-'.isomir.expression.tsv'

tumors<-c('PAAD','PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM') 

tcga_file<-paste(base_path, 'most-expressed-TCGA-consensus.tsv', sep='')
TCGA_data <- data.frame()



for(t in tumors){
  start_time<-Sys.time()
    
  expression_file<-paste(base_path, t, '/summary_files/', t, expression_ext, sep='')
  expression_data<-read.table(file=expression_file, header=TRUE, sep='\t')
  
  expression_data<-expression_data[, c('MIRNA', 'SEQUENCE', 'CPM', 'RATIO')]
  expression_data$RATIO<-as.numeric(gsub('%', '', expression_data$RATIO))
  expression_data$RATIO<-expression_data$RATIO/100
  expression_data<-expression_data[complete.cases(expression_data),]
  expression_data<-cbind(expression_data, data.frame(INSTANCES=1))
  
  expression_data<-expression_data[order(expression_data$MIRNA, decreasing=TRUE),]
  expression_data<-expression_data[order(expression_data$CPM, decreasing=TRUE),]
  unique_mirna<-as.vector(unique(expression_data$MIRNA))
  temp_exp_data<-data.frame(MIRNA=character(), SEQUENCE=character(), CPM=numeric(), RATIO=numeric())
  i<-1
  while(length(unique_mirna)>0){
    mirna<-toString(expression_data[i, 'MIRNA'])
    seq<-expression_data[i, 'SEQUENCE']
    cpm<-expression_data[i, 'CPM']
    ratio<-expression_data[i, 'RATIO']
    
    if(mirna%in%unique_mirna){
      
      curr_data<-expression_data[i,]
      
      i<-i+1
      next_cpm<-expression_data[i, 'CPM']
      next_ratio<-expression_data[i, 'RATIO']
      next_mirna<-expression_data[i, 'MIRNA']
      
      while(next_cpm>=cpm & next_mirna%in%mirna){
        
        next_seq<-expression_data[i, 'MIRNA']
        if(next_seq%in%curr_data){
          curr_data[which(curr_data$SEQUENCE==next_seq), 'RATIO']<-curr_data[which(curr_data$SEQUENCE==next_seq), 'RATIO']+next_ratio
          curr_data[which(curr_data$SEQUENCE==next_seq), 'CPM']<-curr_data[which(curr_data$SEQUENCE==next_seq), 'CPM']+next_cpm
          curr_data[which(curr_data$SEQUENCE==next_seq), 'INSTANCES']<-curr_data[which(curr_data$SEQUENCE==next_seq), 'INSTANCES']+1
        }else{
          new_row<-expression_data[i,]
          curr_data<-rbind(curr_data, new_row)
        }
        
        i<-i+1
        next_cpm<-expression_data[i, 'CPM']
        next_mirna<-expression_data[i, 'MIRNA']
      }
      
      temp_exp_data<-rbind(temp_exp_data, curr_data)
      unique_mirna<-unique_mirna[-which(mirna%in%unique_mirna)]
      
    }else{
      i<-i+1
    }
  }
  expression_data<-temp_exp_data
  # unique_mirna<-unique(expression_data$MIRNA)
  # for(m in unique_mirna){
  #   temp_data<-expression_data[which(expression_data$MIRNA==m),]
  #   if(nrow(temp_data)>0){
  #     most_expressed<-temp_data[which(temp_data$CPM>=max(temp_data$CPM)),]
  #     expression_data<-expression_data[-which(expression_data$MIRNA==m),]
  #     expression_data<-rbind(expression_data, most_expressed)
  #   }
  #   
  # }
  
  # unique_mirna<-unique(expression_data$MIRNA)
  # count<-0
  # for(m in unique_mirna){
  #   count<-count+1
  #   if(count%%200==0){print(count/length(unique_mirna))}
  #   temp_exp_data<-expression_data[which(expression_data$MIRNA==m),]
  #   expression_data<-expression_data[-which(expression_data$MIRNA==m),]
  #   max_cpm<-max(temp_exp_data$CPM)
  #   temp_exp_data<-temp_exp_data[which(temp_exp_data$CPM>=max_cpm/3),]
  #   expression_data<-rbind(expression_data, temp_exp_data)
  #   
  # }
    
  unique_sequences<-unique(expression_data$SEQUENCE)
  unique_sequences<-unique_sequences[!is.na(unique_sequences)]
  if(length(unique_sequences)>100000){
    for(s in unique_sequences[1:50000]){
     temp_exp_data<-expression_data[which(expression_data$SEQUENCE==s),]
     if(nrow(temp_exp_data)>1){
       sample<-temp_exp_data[1, 'SAMPLE']
       mirna<-temp_exp_data[1, 'MIRNA']
       cpm<-sum(as.numeric(temp_exp_data$CPM), na.rm=TRUE)
       ratio<-sum(as.numeric(temp_exp_data$RATIO), na.rm=TRUE)
       instances<-sum(as.numeric(temp_exp_data$TEMP.INSTANCES), na.rm=TRUE)
       expression_data<-expression_data[-which(expression_data$SEQUENCE==s),]
       expression_data<-rbind.fill(expression_data, data.frame(SAMPLE=sample, MIRNA=mirna, SEQUENCE=s, CPM=cpm, RATIO=ratio, INSTANCES=instances))
     }
    }
  }
       
    
  names(expression_data)[names(expression_data)=='CPM']<-'TEMP.CPM'
  names(expression_data)[names(expression_data)=='RATIO']<-'TEMP.RATIO'
  names(expression_data)[names(expression_data)=='INSTANCES']<-'TEMP.INSTANCES'
  
  # unique_samples<-expression_data
  # for(s in unique_samples){
  #   temp_exp_data<-expression_data[which(expression_data$SAMPLE==s),]
  #   temp_exp_data$SAMPLE<-NULL
  TCGA_data<-merge(TCGA_data, expression_data, by=c('SEQUENCE', 'MIRNA'), all=TRUE, na.rm=TRUE)
  TCGA_data$CPM<-rowSums(TCGA_data[,c('CPM', 'TEMP.CPM')], na.rm=TRUE)
  TCGA_data$TEMP.CPM<-NULL
  TCGA_data$RATIO<-rowSums(TCGA_data[,c('RATIO', 'TEMP.RATIO')], na.rm=TRUE)
  TCGA_data$TEMP.RATIO<-NULL
  TCGA_data$INSTANCES<-rowSums(TCGA_data[,c('INSTANCES', 'TEMP.INSTANCES')], na.rm=TRUE)
  TCGA_data$TEMP.INSTANCES<-NULL
  # }
  
  if(nrow(TCGA_data)>100000){
    unique_mirna<-unique(TCGA_data$MIRNA)
    for(m in unique_mirna){
      temp_data<-TCGA_data[which(TCGA_data$MIRNA==m),]
      if(nrow(temp_data)>0){
        most_expressed<-temp_data[which(temp_data$CPM>=max(temp_data$CPM)),]
        TCGA_data<-TCGA_data[-which(TCGA_data$MIRNA==m),]
        TCGA_data<-rbind(TCGA_data, most_expressed)
      }
    }
  }

  #get elapsed time
  print(t)
  print(Sys.time()-start_time)
  
  write.table(TCGA_data, file=tcga_file, sep='\t', row.names=FALSE)
}

if(file.exists(tcga_file)){
  TCGA_data<-read.table(file=tcga_file, header=TRUE, sep='\t')
}

TCGA_data<-TCGA_data[!is.na(TCGA_data$SEQUENCE),]
TCGA_data<-TCGA_data[!is.na(TCGA_data$MIRNA),]

unique_sequences<-unique(expression_data$SEQUENCE)
unique_sequences<-unique_sequences[!is.na(unique_sequences)]
for(s in unique_sequences){
  temp_exp_data<-expression_data[which(expression_data$SEQUENCE==s),]
  temp_exp_data<-unique(temp_exp_data)
  if(nrow(temp_exp_data)>1){
    mirna<-temp_exp_data[1, 'MIRNA']
    cpm<-sum(as.numeric(temp_exp_data$CPM), na.rm=TRUE)
    ratio<-sum(as.numeric(temp_exp_data$RATIO), na.rm=TRUE)
    instances<-sum(as.numeric(temp_exp_data$INSTANCES), na.rm=TRUE)
    expression_data<-expression_data[-which(expression_data$SEQUENCE==s),]
    expression_data<-rbind.fill(expression_data, data.frame(MIRNA=mirna, SEQUENCE=s, CPM=cpm, RATIO=ratio, INSTANCES=instances))
  }
}

unique_mirna<-unique(TCGA_data$MIRNA)
for(m in unique_mirna){

  temp_data<-TCGA_data[which(TCGA_data$MIRNA==m),]
  temp_data<-temp_data[which(temp_data$CPM>=max(temp_data$CPM, na.rm=TRUE)),]
  temp_data<-unique(temp_data)

  if(nrow(temp_data)>1){
    sequence<-paste(temp_data$SEQUENCE, collapse=',')
    cpm<-sum(temp_data$CPM, na.rm=TRUE)
    ratio<-sum(temp_data$CPM, na.rm=TRUE)
    instances<-sum(temp_data$INSTANCES, na.rm=TRUE)
    temp_data<-temp_data[1,]
    temp_data[, 'SEQUENCE']<-as.factor(sequence)
    temp_data[, 'MULTIPLE.SEQ']<-TRUE
    temp_data[,'CPM']<-cpm
    temp_data[,'RATIO']<-ratio
    temp_data[,'INSTANCES']<-instances
  }else{
    temp_data[, 'MULTIPLE.SEQ']<-FALSE
  }

  TCGA_data<-TCGA_data[-which(TCGA_data$MIRNA==m),]
  TCGA_data<-rbind(TCGA_data, temp_data)

}

TCGA_data$RATIO<-TCGA_data$RATIO/(TCGA_data$INSTANCES)
TCGA_data$CPM<-TCGA_data$CPM/(TCGA_data$INSTANCES)
TCGA_data$INSTANCES<-NULL

mirbase_data<-read.table(file=mirna_file, header=TRUE, sep='\t')

mirbase_data<-mirbase_data[,c('MIRNA', 'ACCESSION', 'MOTIF.13', 'SEQUENCE')]
names(mirbase_data)[names(mirbase_data)=='SEQUENCE']<-'MIRBASE21.SEQUENCE'

data<-merge(TCGA_data, mirbase_data, by='MIRNA', all=TRUE, na.rm=TRUE)
data<-data[!is.na(data$SEQUENCE),]
data<-unique(data)
names(data)[names(data)=='SEQUENCE']<-'TCGA.SEQUENCE'
data<-cbind(data, data.frame(TCGA.TEMPLATED=NA))

for(i in 1:nrow(data)){
  multiple_seq<- data[i, 'MULTIPLE.SEQ']
  if(!multiple_seq){
    mirna<-data[i, 'MIRNA']
    mirbase_seq<-as.character(data[i, 'MIRBASE21.SEQUENCE'])
    tcga_seq<-as.character(data[i, 'TCGA.SEQUENCE'])

    data[i, 'IDENTICAL.SEQ']<-identical(mirbase_seq,tcga_seq)

    pri_seq<-mirbase_data[which(mirbase_data$MIRNA==mirna),'MIRBASE21.SEQUENCE']
    for(s in pri_seq){
      if(grepl(tcga_seq, s)){
        data[i, 'TCGA.TEMPLATED']<-TRUE
        break
      }else{
        data[i, 'TCGA.TEMPLATED']<-FALSE
      }
    }
    data[i, 'DISTANCE.MIRBASE21.TCGA']<-stringdist(mirbase_seq, tcga_seq)
  }

}

filename<-paste(base_path, 'updated-TCGA-consensus.tsv', sep='')
write.table(data, file=filename, sep='\t', row.names=FALSE)

consensus_data<-c()

unique_mirna<-unique(data$MIRNA)

for(m in unique_mirna){

  temp_data<-data[which(data$MIRNA==m),]
  templated<-temp_data[,'TCGA.TEMPLATED']
  identical<-temp_data[,'IDENTICAL.SEQ']

  mirbase_seq<-temp_data[,'MIRBASE21.SEQUENCE']
  tcga_seq<-temp_data[,'TCGA.SEQUENCE']

  if(identical%in%'TRUE'){ seq<-mirbase_seq }
  else if(templated%in%'TRUE'){ seq<-tcga_seq }
  else if(templated%in%'FALSE'){ seq<-mirbase_seq }

  motif<-temp_data[,'MOTIF.13']

  first_line<-paste(toString(m), toString(motif), sep='\t')
  second_line<-toString(seq)

  consensus_data<-c(consensus_data, first_line, second_line)
}

consensus_file<-paste(base_path, 'updated-TCGA-consensus.fa', sep='')
file<-file(consensus_file)
writeLines(consensus_data, file)
close(file)
