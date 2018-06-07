## PURPOSE: normalize gene expression files
## INPUT: gene expression files
## OUTPUT: table containing normalized gene data
# # require(dplyr)
# # require(gdata)
# # 
# # install.packages("clusterSim")
# # source("https://bioconductor.org/biocLite.R")
# # biocLite("edgeR")
# # biocLite("org.Hs.eg.db")
# # biocLite("AnnotationDbi")
# # 
# # library(clusterSim)
# # library(limma)
# # library(edgeR)
# # library("AnnotationDbi")
# # library("org.Hs.eg.db")
# # #install.packages("gtools")
# # library(gtools)

## https://www.bioconductor.org/help/workflows/RNAseq123/



directory <- list.files("/Users/chens22/Documents/", full.names = T)

for(f in directory){
  if((!dir.exists(paste0(f, '/gene_expression/'))) || length(list.files(paste0(f, '/gene_expression/'))) < 1){ next }
  tumor <- strsplit(f, '/')[[1]]
  tumor <- tumor[length(tumor)]
  
  setwd(paste0(f))
  disease <- c()
  if(tumor %in% c('TARGET')){ next }#disease <- list.files(getwd())}
  #else{next}
  disease <- c('')
  for(d in disease){
    #setwd(paste0(f, '/gene_expression/', d, '/'))
    if(file.exists(paste0(getwd(), '/', tumor, '.normalized.FKQM-UQ.tsv')) || file.exists(paste0(getwd(), '/', tumor, '.normalized.FPKM-UQ.tsv'))){ next }
    print(tumor)
    
  
    ## Importing data
    setwd(paste0(f, '/gene_expression/', d, '/'))
    fileNames <- Sys.glob("*.FPKM-UQ.txt.gz")
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
    normalizedFPKM <- as.data.frame(cpm)
    #plotMDS(lcpm)
    
    setwd(paste0(f))
    datapatients <- read.csv(paste0(getwd(), '/', tumor,'-gene-expression-manifest.csv'), header = TRUE)
    
    nameid <- datapatients[,c('name','id')]
    nameid$name <- as.character(nameid$name)
    nameid$name <- substr(nameid$name, 1, 48)
    class(nameid$name)
    tnormalizedFPKM <- as.data.frame(t(normalizedFPKM))
    tnormalizedFPKM$name <- rownames(tnormalizedFPKM)
    tnormalizedFPKM <- merge(tnormalizedFPKM, nameid, by="name")
    
    
    write.table(tnormalizedFPKM, paste0(getwd(), '/', tumor, '.normalized.FPKM-UQ.tsv'), sep="\t", append = FALSE, row.names = F)
    
  }
  
}

