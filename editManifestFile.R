## PURPOSE: format manifest file
## INPUT: manifest data from CGC	tumor.manifest.tsv
##		  file containing manifest data for all tumors 		all-tumor-manifest.tsv
## OUTPUT: formatted table containing manifest data for all tumors 	all-tumor-manifest.tsv

require(data.table)
require(plyr)

base_path <- '/Volumes/2TB (MAC)/Susanna/'#'/media/user/2TB (MAC)/Susanna/'

tumor <- 'KIRC'
manifest_file <- paste0(base_path, tumor, '/', tumor, '.manifest.csv')
manifest_data<- read.csv(manifest_file, header=T, stringsAsFactors=F)
manifest_data <- data.frame(lapply(manifest_data, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
for(i in 1:ncol(manifest_data)){
  manifest_data[,i] <- gsub(' ', '.', manifest_data[,i], fixed=T)
}

colnames(manifest_data) <- toupper(colnames(manifest_data))
colnames(manifest_data) <- gsub('_', '.', colnames(manifest_data), fixed=T)



all_manifest_file <- paste0(base_path, 'all-tumor-manifest.csv')
all_manifest_data <- read.table(file=all_manifest_file, sep=',', stringsAsFactors=F, header=T)

manifest_data <- merge(manifest_data, all_manifest_data[,c('INVESTIGATION', 'DISEASE.ABBV')], by='INVESTIGATION')

all_manifest_data <- rbind.fill(all_manifest_data, manifest_data)
all_manifest_data <- unique(all_manifest_data)

for(i in 1:nrow(all_manifest_data)){
  name <- all_manifest_data[i, 'NAME']
  new_name <- strsplit(toString(name), '.', fixed=T)[[1]][1]
  all_manifest_data[which(all_manifest_data$NAME==name), 'NAME'] <- new_name
}

write.table(all_manifest_data, file=all_manifest_file, sep=',', row.names=F)
