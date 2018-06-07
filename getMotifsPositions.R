
## PURPOSE: get position of miRNA motifs in primiRNA
## INPUT: miRBase data  miRBase21-master.tsv
## OUTPUT: table containing motif start positions

mirna_file<-'/Users/chens22/Documents/miRNA/miRBase21-master.tsv'
data<-read.table(file=mirna_file, header=TRUE, sep='\t', stringsAsFactors = FALSE)
data<-cbind(data, data.frame(UG.5P.START=NA, CNNC.START=NA, UGUG.START=NA, UG.LOOP.START=NA, HAS.5P=NA, HAS.3P=NA))
#data<-mirna_data[,c('MIRNA', 'SEQUENCE', 'EXTENDED.SEQUENCE', )]
unique_acc<-unique(data$PRI.ACCESSION)

for(a in unique_acc){
  
  seq_5p<-toString(data[which(data$PRI.ACCESSION==a & data$STRAND=='5P'), 'SEQUENCE'])
  seq_3p<-toString(data[which(data$PRI.ACCESSION==a & data$STRAND=='3P'), 'SEQUENCE'])
  extended_seq<-toString(unique(data[which(data$PRI.ACCESSION==a), 'EXTENDED.SEQUENCE']))
  pri_seq<-toString(unique(data[which(data$PRI.ACCESSION==a), 'PRI.SEQUENCE']))
  
  if(nchar(seq_5p)<1){
    data[which(data$PRI.ACCESSION==a), 'HAS.5P']<-FALSE
    has_5p<-FALSE
  }else{
    data[which(data$PRI.ACCESSION==a), 'HAS.5P']<-TRUE
    has_5p<-TRUE
  }
  
  if(nchar(seq_3p)<1){
    data[which(data$PRI.ACCESSION==a), 'HAS.3P']<-FALSE
    has_3p<-FALSE
  }else{
    data[which(data$PRI.ACCESSION==a), 'HAS.3P']<-TRUE
    has_3p<-TRUE
  }
  
  #find position of UG motif before 5p strand
  if(has_5p){
    start<-gregexpr(seq_5p, extended_seq, fixed=TRUE)[[1]][1]-20
    end<-gregexpr(seq_5p, extended_seq, fixed=TRUE)[[1]][1]-9
    sub_seq_5p<-substr(extended_seq, start, end)
    motif_pos<-gregexpr('TG', sub_seq_5p)[[1]]
    motif_pos<-gregexpr('TG', sub_seq_5p)[[1]][length(motif_pos)]
    if(motif_pos!=-1){
      data[which(data$PRI.ACCESSION==a), 'UG.5P.START']<-((12-motif_pos+9))*-1
    }
  }
  
  #find position of CNNC motif after 3p strand
  if(has_3p){
    start<-gregexpr(seq_3p, extended_seq, fixed=TRUE)[[1]][1]+nchar(seq_3p)+9
    end<-gregexpr(seq_3p, extended_seq, fixed=TRUE)[[1]][1]+nchar(seq_3p)+20
    sub_seq_3p<-substr(extended_seq, start, end)
    for(i in 1:nchar(sub_seq_3p)-4){
      temp_seq<-substr(sub_seq_3p, i, i+4)
      first_char<-substr(temp_seq, 1, 1)
      last_char<-substr(temp_seq, 4, 4)
      if(first_char%in%'C' && last_char%in%'C'){
        data[which(data$PRI.ACCESSION==a), 'CNNC.START']<-i+9
      }
    }
  }
    
  #find position of UG and UGUG motif in the loop (between 5p and 3p strand)
  if(has_5p && has_3p){
    start<-gregexpr(seq_5p, pri_seq, fixed=TRUE)[[1]][1]+nchar(seq_5p)
    end<-gregexpr(seq_3p, pri_seq, fixed=TRUE)[[1]][1]
    sub_seq_loop<-substr(pri_seq, start, end)
    motif_pos<-gregexpr('TG', sub_seq_loop)[[1]][1]
    if(motif_pos!=-1){
      data[which(data$PRI.ACCESSION==a), 'UG.LOOP.START']<-motif_pos
    }
    motif_pos<-gregexpr('TGTG', sub_seq_loop)[[1]][1]
    if(motif_pos!=-1){
      data[which(data$PRI.ACCESSION==a), 'UGUG.START']<-motif_pos
    }
  }
}

write.table(data, file=mirna_file, sep='\t', row.names=FALSE)





