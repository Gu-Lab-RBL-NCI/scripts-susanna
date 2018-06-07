
library("AnnotationDbi")
library("org.Hs.eg.db")
library(stringr)
library(data.table)

## PURPOSE: FIND DIFFERENCE BETWEEN FIDELITY OF SAMPLES WITH HIGH AND LOW LEVELS OF CERTAIN GENE EXPRESSION
## INPUT: formatted manifest file from CGC                 all-tumor-manifest.csv
##       file containing average RPM of all miRNA in tumor    tumor.fidelity.tsv
##       normalized FPKM-UQ file                    tumor.normalized.FPKM-UQ.tsv
##       normalized FPKM-UQ manifest file     tumor-gene-expression-manifest.csv
##       isomiR summary data                                    tumor.isomir.tsv

## OUTPUT: p value of effective number of samples with highest and lowest gene expression   tumor.effective.gene.expression.all.miRNA.p.p_limit.tsv

base_path<-'/Users/chens22/Documents/'
setwd(base_path)

manifest_file<-'all-tumor-manifest.csv'
manifest_data<-read.csv(file=manifest_file, header=T, stringsAsFactors=F)

summary_ext<-'.isomir.tsv'
expression_ext<-'.isomir.expression.tsv'
cleavage_ext<-'.isomir.cleavage.tsv'
fidelity_ext<-'.fidelity.tsv'
gene_ext <- '.normalized.FPKM-UQ.tsv'


num_samples <- 25
p_limit <- 0.05


tumors <- c('UCEC')


for(tumor in tumors){
  
  
  #get 50 miRNA with highest rpm
  fidelity_file <- paste0(base_path, tumor, '/fidelity/', tumor, '.fidelity.tsv')
  fidelity_data <- read.table(fidelity_file, sep='\t', header=T)

  fidelity_data <- fidelity_data[order(-fidelity_data[,'RPM.AVERAGE']),]
  mirna <- fidelity_data[1:50, 'MIRNA']

  gene_file <- paste0('/', tumor, '/expression/', tumor, gene_ext)
  gene_data <- read.table(gene_file, sep='\t', header=T, stringsAsFactors = F)

  gene_expression <- gene_data
  gene_expression[,c('id', 'name', 'sample_id')] <- NULL
  gene_expression <- transpose(gene_expression)
  gene_expression$EXPRESSION <- rowMeans(data.matrix(gene_expression), na.rm=T)
  gene_expression$GENE <- colnames(gene_data)[!colnames(gene_data)%in%c('id', 'name', 'sample_id')]
  gene_expression <- gene_expression[order(-gene_expression$EXPRESSION),]
  genes <- gene_expression[1:as.integer(nrow(gene_expression)/4), 'GENE']
  
  genes <- sort(genes)
  genes <- rev(genes)
  

  gene_manifest_file <- paste0('/', tumor, '/expression/', tumor, '-gene-expression-manifest.csv')
  gene_manifest_data <- read.table(gene_manifest_file, sep=',', header=T, stringsAsFactors = F)
  if(!'sample_id'%in%colnames(gene_data)){
    gene_data <- merge(gene_data, gene_manifest_data[c('name', 'sample_id')], by='name',all=T)
    
  }
  
  summary_file<-paste0('/', tumor, '/summary_files/', tumor, '.isomir.tsv')
  summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors = F)
  
  
  file <- paste0('/', tumor,'/expression/', tumor, '.effective.gene.expression.all.miRNA.p.', p_limit, '.tsv')
  
  if(file.exists(file)){
    data <- read.table(file=file, sep='\t', header=T, stringsAsFactors=F) #read precalculated output file
  }else{
    data <- data.frame(GENE=character(), stringsAsFactors=F)
  }
  genes <- genes[!genes%in%c('id', 'name', 'sample_id')]

  for(g in genes){
    
    if(g%in%data$GENE){ next } #prevent repeat calculations

    row <- nrow(data)+1
    
    # get mean top and bottom gene expression for gene
    data[row, 'GENE'] <- as.character(g)
    gene_data <- gene_data[order(-gene_data[,g]),]
    gene_expression <- mean(gene_data[1:num_samples, g])
    data[row, 'TOP.GENE.EXPRESSION'] <- gene_expression
    
    gene_data <- gene_data[order(gene_data[,g]),]
    gene_expression <- mean(gene_data[1:num_samples, g])
    data[row, 'BOTTOM.GENE.EXPRESSION'] <- gene_expression
    
    for(m in mirna){
      # acquire effective number of samples in the top gene expression
      m <- gsub('.', '-', m, fixed=T)
      gene_data <- gene_data[order(-gene_data[,g]),]
      sample_ids <- gene_data[1:num_samples, 'sample_id']
      samples <- manifest_data[which(manifest_data$SAMPLE.ID%in%sample_ids), 'NAME']
      temp_summ_data <- summary_data[which(summary_data$SAMPLE%in%samples & summary_data$MIRNA==m),]
      
      top_effective <- temp_summ_data[, 'EFFECTIVE.NUMBER']
      
      # acquire effective number of samples in the bottom gene expression
      m <- gsub('.', '-', m, fixed=T)
      gene_data <- gene_data[order(gene_data[,g]),]
      sample_ids <- gene_data[1:num_samples, 'sample_id']
      samples <- manifest_data[which(manifest_data$SAMPLE.ID%in%sample_ids), 'NAME']
      temp_summ_data <- summary_data[which(summary_data$SAMPLE%in%samples & summary_data$MIRNA==m),]
      
      bottom_effective <- temp_summ_data[, 'EFFECTIVE.NUMBER']
      
      # calculate p value of both sets of effective numbers if possible
      m <- gsub('-', '.', m, fixed=T)
      if((length(unique(bottom_effective)) <= 1 || length(unique(top_effective)) <= 1)){
        data[row, m] <- NA
      }else{
        data[row, m] <- t.test(top_effective, bottom_effective)$p.value
      }
  
      data[row, paste0(m, '.FOLD.CHANGE')] <- mean(top_effective)/mean(bottom_effective)
      
      
    }

    write.table(data, file=file, quote=FALSE, sep='\t', row.names = FALSE)
  }
  
  
}

# PREVIOUSLY USED PROGRAM FOR FORMATTING DATA TABLE
# 
# tumors <- c('LAML')
# for(tumor in tumors){
#   data <- data.frame(GENE=character())
#   row_len <- Inf
#   for(m in mirna){
#     data_file <- paste0('/', tumor,'/expression/', tumor, '.', m, '.effective.gene.expression.tsv')
#     curr_data <- read.table(file=data_file, sep='\t', header=T, stringsAsFactors=F)
# 
#     curr_data$FOLD.CHANGE <-  curr_data$TOP.EFFECTIVE.NUMBER/curr_data$BOTTOM.EFFECTIVE.NUMBER
#     data <- merge(data, curr_data[,c('GENE', 'EFFECTIVE.P', 'FOLD.CHANGE')], by='GENE', all=T, na.rm=T)
#     colnames(data)[colnames(data)=='EFFECTIVE.P'] <- m
#     colnames(data)[colnames(data)=='FOLD.CHANGE'] <- paste0(m, '.FOLD.CHANGE')
#   }
#   
#   
#   data$SYMBOL <- mapIds(org.Hs.eg.db, keys=str_split_fixed(data$GENE, '\\.', 2)[,1], column="SYMBOL", keytype="ENSEMBL" ,multiVals="first")
#   write.table(data, file= paste0('/', tumor,'/expression/', tumor, '.effective.gene.expression.tsv'), sep='\t', row.names=F)
#   
# 
#   for(i in 1:nrow(data)){
#     small_p <- c(data[i,mirna[1]] <= p_limit, data[i,mirna[2]] <= p_limit, data[i,mirna[3]] <= p_limit)
#     data[i, 'SMALL.P'] <- length(small_p[small_p==T])
#   }
#   
#   data$DIRECTION <- (data[,paste0(mirna[1], '.FOLD.CHANGE')] > 1 &
#                        data[,paste0(mirna[2], '.FOLD.CHANGE')] > 1 &
#                        data[,paste0(mirna[3], '.FOLD.CHANGE')] > 1) |
#     (data[,paste0(mirna[1], '.FOLD.CHANGE')] < 1 &
#        data[,paste0(mirna[2], '.FOLD.CHANGE')] < 1 &
#        data[,paste0(mirna[3], '.FOLD.CHANGE')] < 1)
#   
#   
#   file <- paste0('/', tumor,'/expression/', tumor, '.effective.gene.expression.p.', p_limit, '.tsv')
#   write.table(data, file=file, sep='\t', row.names=F)
#   
#   
# }
# 
# 
# tumors <- c('BLCA', 'PRAD', 'UCEC', 'READ', 'LAML', 'OV', 'LUSC', 'KIRC')
# data <- data.frame(GENE = character(), SYMBOL=character())
# for(t in tumors){
#   print(t)
#   file <- paste0('/', t,'/expression/', t, '.effective.gene.expression.p.', p_limit, '.tsv')
#   curr_data <- read.table(file=file, sep='\t', header=T, stringsAsFactors=F)
#   curr_data <- unique(curr_data)
#   curr_data <- curr_data[which(curr_data$SMALL.P>0),]
#   for(m in mirna){
#     abbv <- strsplit(m, '-', fixed=T)[[1]][3]
#     if(F%in%grepl(toString(m), colnames(curr_data))){    m<- gsub('-', '.', m, fixed=T)}
# 
#     colnames(curr_data)[colnames(curr_data)==paste0(m, '.FOLD.CHANGE')] <-  paste0(t, '.', abbv,'.FC')
#     colnames(curr_data)[colnames(curr_data)==m] <-  paste0(t, '.', abbv, '.P')
#                                  
#   }
#   colnames(curr_data)[colnames(curr_data)=='SMALL.P'] <-  paste0(t, '.SMALL.P')
#   colnames(curr_data)[colnames(curr_data)=='DIRECTION'] <-  paste0(t, '.DIRECTION')
#   
# 
#   if(nrow(data) < 1){
#     data <- merge(data, curr_data, all=T, na.rm=T, by=c('GENE', 'SYMBOL'))
#   }else{
#     data <- merge(data, curr_data, all=F,na.rm=T, by=c('GENE', 'SYMBOL'))
#   }
# 
#   
# }
# data <- read.table(file=paste0('MECP2.', paste(tumors, collapse='.'), '.effective.gene.expression.tsv'), sep='\t', header=T, stringsAsFactors=F)
# data$TOTAL.SMALL.P <- rowSums(data[ , grepl( ".SMALL.P" , names( data ) ) ])
# write.table(data, file=paste0(base_path, paste(tumors, collapse='.'), '.effective.gene.expression.tsv'), sep='\t', row.names = F)
# data <- read.table(file=paste0(base_path, paste(tumors, collapse='.'), '.effective.gene.expression.tsv'), sep='\t', header=T, stringsAsFactors=F)
# data <- data[complete.cases(data[,colnames(data)!='SYMBOL']),]
# 
# 
