## PURPOSE: format miRBase data
## INPUT: miRBase summary file                          miRBase21-summary.csv
##        miRBase results file (motif, paralogs, etc.)  miRBase21-results.tsv
##        UCSC extended sequences           miRNA-extended-sequences-UCSC.tsv
## OUTPUT: formatted miRBase data                         miRBase21-master.tsv


summary_file<-'/Users/chens22/Documents/miRNA/RNA Project/miRBase21-summary.csv'
result_file<-'/Users/chens22/Dropbox/miRNA Analysis/miRBase21-results.tsv'
extended_seq_file<-'/Users/chens22/Documents/miRNA/miRNA-extended-sequences-UCSC.tsv'

summary_data<-read.csv(file=summary_file, header=TRUE)

mirna_data<-data.frame(PRIMIRNA=NA, 
                       PRI.ACCESSION=NA, 
                       PRI.SEQUENCE=NA,
                       ACCESSION=NA,
                       MIRNA=NA,
                       SEQUENCE=NA)

for(i in 1:nrow(summary_data)){
  if(nchar(toString(summary_data[i, 'Mature1_Acc']))>1){
    mature1<-summary_data[i, c('PRIMIRNA', 'ACCESSION', 
                               'SEQUENCE', 'Mature1_Acc',
                               'Mature1_ID', 'Mature1_Seq')]
    names(mature1)[names(mature1) == 'ACCESSION'] <- 'PRI.ACCESSION'
    names(mature1)[names(mature1) == 'SEQUENCE'] <- 'PRI.SEQUENCE'
    names(mature1)[names(mature1) == 'Mature1_Acc'] <- 'ACCESSION'
    names(mature1)[names(mature1) == 'Mature1_ID'] <- 'MIRNA'
    names(mature1)[names(mature1) == 'Mature1_Seq'] <- 'SEQUENCE'
    
    mirna_data<-rbind(mirna_data, mature1)
  }
  
  if(nchar(toString(summary_data[i, 'Mature2_Acc']))>1){
    mature2<-summary_data[i, c('PRIMIRNA', 'ACCESSION', 
                             'SEQUENCE', 'Mature2_Acc',
                             'Mature2_ID', 'Mature2_Seq')]
    names(mature2)[names(mature2) == 'ACCESSION'] <- 'PRI.ACCESSION'
    names(mature2)[names(mature2) == 'SEQUENCE'] <- 'PRI.SEQUENCE'
    names(mature2)[names(mature2) == 'Mature2_Acc'] <- 'ACCESSION'
    names(mature2)[names(mature2) == 'Mature2_ID'] <- 'MIRNA'
    names(mature2)[names(mature2) == 'Mature2_Seq'] <- 'SEQUENCE'
    
    mirna_data<-rbind(mirna_data, mature2)
  }
}


result_data<-read.table(file=result_file, header=TRUE, sep='\t')

result_data<-result_data[,c('MATURE.ID', 'PARALOGS', 'MOTIF.13', 'N.MOTIF',
                              'NON.N', 'DUPLIMOTIF', 'DUPLI.ID', 'UNIQUE.MOTIF',
                              'MOTIF.LEN', 'SEED', 'FAMILY', 'STRAND')]
names(result_data)[names(result_data) == 'MATURE.ID'] <- 'MIRNA'

mirna_data<-merge(mirna_data, result_data, by='MIRNA', all=TRUE, na.rm=TRUE)

#recreate mirna ids with paralogs

for(i in 1:nrow(mirna_data)){
  
  paralogs<-mirna_data[i, 'PARALOGS']
  
  if(!is.na(paralogs) & paralogs>1){
    mirna_id<-mirna_data[i, 'MIRNA']
    mirna_id<-paste(mirna_id, paste(1:paralogs, collapse='-'), sep='-')
    mirna_data[i, 'MIRNA']<-mirna_id
  }
  
}


extended_data<-read.table(file=extended_seq_file, sep='\t', header=TRUE)
extended_data$PRIMIRNA<-NULL
names(extended_data)[names(extended_data) == 'ACCESSION'] <- 'PRI.ACCESSION'
names(extended_data)[names(extended_data) == 'STRAND'] <- 'DIRECTION'

mirna_data<-merge(mirna_data, extended_data, by='PRI.ACCESSION', all=TRUE, na.rm=TRUE)

mirna_data$SEQUENCE<-gsub('U', 'T', mirna_data$SEQUENCE)
mirna_data$PRI.SEQUENCE<-gsub('U', 'T', mirna_data$PRI.SEQUENCE)

write.table(mirna_data, file='/Users/chens22/Documents/miRNA/miRBase21-master.tsv', row.names=FALSE, sep='\t')

