
## PURPOSE: calculate new EFFECTIVE.NUMBER column in tumor.isomir.tsv summary files
## INPUT:	isomiR summary file 				tumor.isomir.tsv
##			isomiR expression file 	tumor.isomir.expression.tsv
## OUTPUT: adds new EFFECTIVE.NUMBER column to tumor.isomir.tsv

base_path <- '/Volumes/2TB (MAC)/Susanna/'

expression_ext<-'.isomir.expression.tsv'
summary_ext<-'.isomir.tsv'

#calculate effective number for set of reads
getEffectiveNumber<-function(reads){
  reads_freq<-reads/sum(reads)
  effective<-1/sum(reads_freq^2)
  return(effective)
 }


for(f in list.files(base_path, full.names=T)){ #iterate through file system


	tumor <- strsplit(toString(f), '/', fixed=T)[[1]]
	tumor <- tumor[length(tumor)]
	

  	f <- paste0(base_path, tumor)

	summary_file <- paste0(base_path, tumor, '/summary_files/', tumor, summary_ext)
	summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors=F)
	#delete duplicate columns
	colnames(summary_data) <- gsub('.x', '', colnames(summary_data), fixed=T) 
	colnames(summary_data) <- gsub('.y', '', colnames(summary_data), fixed=T)
	summary_data <- summary_data[, !duplicated(colnames(summary_data))]
	#edit sample IDs
	for(s in unique(summary_data$SAMPLE)){
	  new_sample <- strsplit(toString(s), '.', fixed=T)[[1]][1]
	  summary_data[which(summary_data$SAMPLE==s), 'SAMPLE'] <- new_sample
	}
	write.table(summary_data, file=summary_file, sep='\t', row.names=F)
	
	#skip if effective number is already calculated in summary file
	if('EFFECTIVE.NUMBER' %in% colnames(summary_data)){ 
	  next 
	}else{ 
	  summary_data$EFFECTIVE.NUMBER <- NA 
	}

	
	expression_file <- paste0(f, '/summary_files/', tumor, expression_ext)
	expression_data <- fread(expression_file, sep='\t', header=T, stringsAsFactors=F)
	#edit sample IDs
	for(s in unique(expression_data$SAMPLE)){
		new_sample <- strsplit(toString(s), '.', fixed=T)[[1]][1]
		expression_data[which(expression_data$SAMPLE==s), 'SAMPLE'] <- new_sample
	}
	expression_data <- expression_data[,c('MIRNA', 'SAMPLE', 'READS', 'SEQUENCE')]
	expression_data$SEED<-substr(expression_data$SEQUENCE, 2,8)
  
	
	effective_data <- as.data.frame(unique(summary_data[, c('SAMPLE', 'MIRNA')]))
	for(s in unique(expression_data$SAMPLE)){

		curr_summ_data <- summary_data[which(summary_data$SAMPLE==toString(s)),c('MIRNA', 'SAMPLE')]
		curr_exp_data <- expression_data[which(expression_data$SAMPLE==toString(s)),]

		seeds_data<-aggregate(curr_exp_data$READS, by=list(MIRNA=curr_exp_data$MIRNA, SEED=curr_exp_data$SEED), sum)
		colnames(seeds_data)[which(names(seeds_data)=='x')]<-'READS'
		seeds_data<-seeds_data[order(-seeds_data$READS),] 
		seeds_data$SAMPLE<-s
		
		mirnas <- summary_data[which(summary_data$SAMPLE==s & is.na(summary_data$EFFECTIVE.NUMBER)==T), 'MIRNA']
		for(m in mirnas){

			curr_eff_num <- getEffectiveNumber(seeds_data[which(seeds_data$MIRNA==m), 'READS'])
			summary_data[which(summary_data$SAMPLE==s & summary_data$MIRNA==m), 'EFFECTIVE.NUMBER'] <- as.numeric(curr_eff_num)
		}
		expression_data <- expression_data[-which(expression_data$SAMPLE==s),]
		write.table(summary_data, file=summary_file, sep='\t', row.names=F)
	}

	if('EFFECTIVE.NUMBER' %in% colnames(summary_data)){
	  summary_data$EFFECTIVE.NUMBER <- NULL
	}

	new_summary_data <- merge(summary_data, effective_data, by=c('SAMPLE', 'MIRNA'))
	summary_data <- read.table(summary_file, sep='\t', header=T, stringsAsFactors=F)
	summary_data <- merge(summary_data, new_summary_data, by=c('SAMPLE', 'MIRNA'))
	write.table(summary_data, file=summary_file, sep='\t', row.names=F)
}
