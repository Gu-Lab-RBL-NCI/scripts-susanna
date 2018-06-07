## PURPOSE: normalize miRNA expression files
## INPUT: miRNA expression files
## OUTPUT: table containing normalized miRNA expression data
require(dplyr)
require(gdata)

install.packages("clusterSim")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("org.Hs.eg.db")
biocLite("AnnotationDbi")

library(clusterSim)
library(limma)
library(edgeR)
library("AnnotationDbi")
library("org.Hs.eg.db")
#install.packages("gtools")
library(gtools)
# 
# ## https://www.bioconductor.org/help/workflows/RNAseq123/

directory <- list.files("/Users/chens22/Documents/", full.names = T)

for(f in directory){
  
  if((!dir.exists(paste0(f, '/miRNA_expression/')))){ next }
  tumor <- strsplit(f, '/')[[1]]
  tumor <- tumor[length(tumor)]

  
  setwd(paste0(f))
  disease <- c('')
  if(tumor %in% c('TARGET')){ next }#disease <- list.files(getwd())}
  else{}
  
  
  for(d in disease){
    
    if(file.exists(paste0(getwd(), '/', tumor, '.normalized.miRNA.tsv'))){ next }
    print(tumor)
   
    datapatients <- read.csv(paste0(getwd(), '/',tumor, '-miRNA-expression-manifest.csv'), header = TRUE)
    setwd(paste0(f, '/miRNA_expression/', d, '/'))
    fileNames <- Sys.glob("*.mirna.quantification.txt")
    temp_file_names <- gsub('.txt', '', fileNames, fixed=T)
    
    
    missing <- datapatients[ !(datapatients$name %in% temp_file_names), ]
    write.table(missing, paste0(tumor, ".missing.miRNA.expression.tsv"), sep="\t", append = FALSE, row.names = TRUE)
    datapatients$mirna <- regexpr(".mirna.quantification.txt", datapatients$name)
    #datapatients <- datapatients[which(datapatients$mirna>-1),]
    
    
    
    ## Importing data
    read.delim(fileNames[1], nrow=5, header = FALSE)
    x <- readDGE(fileNames, columns=c(1,2)) #FPKM-UQ
    
    cpm <- cpm(x)
    lcpm <- cpm(x, log=TRUE)
    
    ## Trimming data with low expression
    samplenum <- length(fileNames)
    table(rowSums(x$counts==0)==samplenum)
    
    keep.exprs <- rowSums(cpm>1)>=samplenum/3
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    dim(x)
    
    ## Normalization
    x <- calcNormFactors(x, method = "TMM")
    x$samples$norm.factors
    lcpm <- cpm(x, log=TRUE)
    cpm <- cpm(x)
    normalizedmiRNA <- as.data.frame(cpm)
    #plotMDS(lcpm)
    
    nameid <- datapatients[,c('name','id')]
    nameid$name <- as.character(nameid$name)
    nameid$name <- gsub('.txt', '', nameid$name, fixed=T)
    tnormalizedmiRNA <- as.data.frame(t(normalizedmiRNA))
    tnormalizedmiRNA$name <- rownames(tnormalizedmiRNA)
    tnormalizedmiRNA$name <- as.character(tnormalizedmiRNA$name)
    tnormalizedmiRNA <- merge(tnormalizedmiRNA, nameid, by="name")
    
    setwd(f)
    write.table(tnormalizedmiRNA, paste0(tumor, '.normalized.miRNA.tsv'), sep="\t", append = FALSE, row.names = F)

  }
}
