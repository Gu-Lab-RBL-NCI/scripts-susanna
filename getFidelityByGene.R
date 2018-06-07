
## PURPOSE: get average fidelity of samples with top expression of specified gene
## INPUT: isomiR summary data     tumor.isomir.tsv
##        isomiR expression data    tumor.isomir.expression.tsv
##        miRBase data              miRBase21-master.tsv
##        manifest data             all-tumor-manifest.csv
##        gene expression data      tumor.normalized.FPKM-UQ.tsv
##        gene expression manifest data   tumor-gene-expression-manifest.csv
## OUTPUT: table with 
library(ggplot2)

mirna_file<-'/Users/chens22/Documents/miRBase21/miRBase21-master.tsv'
mirna_data<-read.table(file=mirna_file, sep='\t', header=TRUE, stringsAsFactors = F)
manifest_file<-'/Users/chens22/Documents/all-tumor-manifest.csv'
manifest_data<-read.csv(file=manifest_file, header=T, stringsAsFactors=F)
base_path<-'/Users/chens22/Documents/'
summary_ext<-'.isomir.tsv'
expression_ext<-'.isomir.expression.tsv'
cleavage_ext<-'.isomir.cleavage.tsv'
fidelity_ext<-'.fidelity.tsv'
gene_ext <- '.normalized.FPKM-UQ.tsv'


num_samples <- 25
p_limit <- 0.05
mirna <- c('hsa-miR-29a-3p', 'hsa-miR-30e-5p', 'hsa-miR-151a-3p', 'hsa-miR-101-3p-1-2', 'hsa-miR-183-5p', 'hsa-miR-1307-3p', 'hsa-let-7i-5p', 'hsa-miR-182-5p', 'hsa-let-7g-5p', 'hsa-miR-181a-2-3p', 'hsa-miR-532-5p', 'hsa-miR-22-3p', 'hsa-miR-181a-5p-1-2', 'hsa-miR-92a-3p-1-2', 'hsa-miR-143-3p', 'hsa-miR-100-5p', 'hsa-miR-125a-5p', 'hsa-miR-125b-5p-1-2', 'hsa-miR-148b-3p', 'hsa-miR-148a-3p', 'hsa-miR-25-3p', 'hsa-let-7b-5p', 'hsa-let-7e-5p', 'hsa-miR-21-5p', 'hsa-miR-99b-5p', 'hsa-miR-361-5p', 'hsa-miR-26a-5p-1-2', 'hsa-let-7f-5p-1-2', 'hsa-miR-30d-5p', 'hsa-miR-30a-5p', 'hsa-miR-24-3p-1-2', 'hsa-miR-191-5p', 'hsa-miR-93-5p', 'hsa-miR-28-3p', 'hsa-miR-10b-5p')
tumors <- rev(c('BLCA', 'PRAD', 'UCEC', 'READ', 'LAML', 'OV', 'KIRC', 'LUSC'))
genes <- c('ENSG00000135486.16')#c('ENSG00000169057.18') #ENSG00000142544.6


for(g in genes){

  gene <-mapIds(org.Hs.eg.db, keys=c(strsplit(g, '.', fixed=T)[[1]][1]), column="SYMBOL", keytype="ENSEMBL" ,multiVals="first")

  file <- paste0(base_path, gene, '.', paste(tumors, collapse='.'), '.all.miRNA.fidelity.tsv')#paste0('/Users/chens22/Documents/', gene, '.', paste(tumors, collapse='.'), '.all.miRNA.effective.tsv')#
  
  if(file.exists(file)){
    data <- read.table(file=file, sep='\t', header=T, stringsAsFactors=F)
  }else{
    data <- data.frame(MIRNA=mirna)
  }
  
  
  
  for(tumor in tumors){
    
    if(paste0(tumor, '.P')%in%colnames(data)){ next }
    
    gene_file <- paste0('/Users/chens22/Documents/', tumor, '/expression/', tumor, gene_ext)
    gene_data <- read.table(gene_file, sep='\t', header=T, stringsAsFactors = F)
    #gene_data <- gene_data[,c('name', gene)]
   # if(!'ENSG00000105676.12'%in%genes){ next }
    
    
    
    gene_manifest_file <- paste0('/Users/chens22/Documents/', tumor, '/expression/', tumor, '-gene-expression-manifest.csv')
    gene_manifest_data <- read.table(gene_manifest_file, sep=',', header=T, stringsAsFactors = F)
    gene_data <- merge(gene_data, gene_manifest_data[c('name', 'sample_id')], by='name',all=T)
    
    
    summary_file<-paste0(base_path, tumor, '/summary_files/', tumor, summary_ext)
    summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors = F)
    
    if(!'EFFECTIVE.NUMBER'%in%colnames(summary_data)){
      summary_file<-paste0(base_path, tumor, '/summary_files/', tumor, '.isomir.effective.tsv')
      summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors = F)
    }
    
    
    genes <- genes[genes!='id']
    if(g%in%data$GENE){ next }
    
    for(m in mirna){
      gene_data <- gene_data[order(-gene_data[,g]),]
      sample_ids <- gene_data[1:num_samples, 'sample_id']
      samples <- manifest_data[which(manifest_data$SAMPLE.ID%in%sample_ids), 'NAME']
      temp_summ_data <- summary_data[which(summary_data$SAMPLE%in%samples & summary_data$MIRNA==m),]
      
      top_effective <- temp_summ_data[, 'FIDELITY.5P'] #temp_summ_data[, 'EFFECTIVE.NUMBER']
      top_fpkm <- mean(gene_data[1:num_samples, g])
      
      gene_data <- gene_data[order(gene_data[,g]),]
      sample_ids <- gene_data[1:num_samples, 'sample_id']
      samples <- manifest_data[which(manifest_data$SAMPLE.ID%in%sample_ids), 'NAME']
      temp_summ_data <- summary_data[which(summary_data$SAMPLE%in%samples & summary_data$MIRNA==m),]
      
      bottom_effective <- temp_summ_data[, 'FIDELITY.5P'] #temp_summ_data[, 'EFFECTIVE.NUMBER']
      bottom_fpkm <- mean(gene_data[1:num_samples, g])
      #m <- gsub('-', '.', m, fixed=T)
      data[which(data$MIRNA==m), paste0(tumor, '.P')] <- t.test(top_effective, bottom_effective)$p.value
      data[which(data$MIRNA==m), paste0(tumor, '.FOLD.CHANGE')] <- top_fpkm/bottom_fpkm#mean(top_effective)/mean(bottom_effective)
      data[which(data$MIRNA==m), paste0(tumor, '.TOP.FPKM')] <- top_fpkm

      data[which(data$MIRNA==m), paste0(tumor, '.BOTTOM.FPKM')] <- bottom_fpkm      
      write.table(data, file=file, quote=FALSE, sep='\t', row.names = FALSE)
      
      
      
      
      
    }
  }
  
  

}
####compare fidelity 5p with effective values#####

fidelity_file <- paste0(base_path, 'CTU1.', paste(tumors, collapse='.'), '.all.miRNA.fidelity.tsv')#paste0('/Users/chens22/Documents/', gene, '.', paste(tumors, collapse='.'), '.all.miRNA.effective.tsv')#
fidelity_data <- read.table(file=fidelity_file, sep='\t', header=T, stringsAsFactors=F)

effective_file <- paste0(base_path, 'CTU1.', paste(tumors, collapse='.'), '.all.miRNA.effective.tsv')#paste0('/Users/chens22/Documents/', gene, '.', paste(tumors, collapse='.'), '.all.miRNA.effective.tsv')#
effective_data <- read.table(file=effective_file, sep='\t', header=T, stringsAsFactors=F)



tumor <- 'OV'

fid_data <- fidelity_data[, c('MIRNA', paste0(tumor, '.P'))]
#fid_data[,paste0(tumor, '.P')] <- fid_data[,paste0(tumor, '.P')]/max(fid_data[,paste0(tumor, '.P')])
colnames(fid_data)[colnames(fid_data)==paste0(tumor, '.P')] <- 'FIDELITY.P'

eff_data <- effective_data[, c('MIRNA', paste0(tumor, '.P'))]
#eff_data[,paste0(tumor, '.P')] <- eff_data[,paste0(tumor, '.P')]/max(eff_data[,paste0(tumor, '.P')])
colnames(eff_data)[colnames(eff_data)==paste0(tumor, '.P')] <- 'EFFECTIVE.P'

data <- merge(fid_data, eff_data, by='MIRNA')
data <- data[order(-data[,'EFFECTIVE.P']),]
data$MIRNA <- factor(data$MIRNA)
p <- ggplot(data, aes(MIRNA), las=2) 
  p+
  geom_point(aes(y=FIDELITY.P)) + 
  geom_path(aes(y=FIDELITY.P), group=1)+geom_point(aes(y=EFFECTIVE.P)) + 
  geom_path(aes(y=EFFECTIVE.P), group=1)+scale_x_discrete(limits=data$MIRNA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

summary_file<-paste0(base_path, tumor, '/summary_files/', tumor, summary_ext)
summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors = F)
summary_data <- summary_data[which(summary_data$MIRNA%in%data$MIRNA),]
p <- ggplot(summary_data, aes(MIRNA), las=2) 
p +
  geom_point(aes(y=FIDELITY.5P)) + 
  #geom_path(aes(y=FIDELITY.5P), group=1) + 
  geom_point(aes(y=EFFECTIVE.NUMBER)) + 
  #geom_path(aes(y=EFFECTIVE.NUMBER), group=1)+
  scale_x_discrete(limits=data$MIRNA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_color_manual(values=c("#CC6666", "#9999CC"))
