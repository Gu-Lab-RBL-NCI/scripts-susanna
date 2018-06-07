## PURPOSE: create heatmap of patient fidelity across miRNA sorted by RPM
## INPUT: isomiR summary data   tumor.isomir.tsv
##        manifest data         all-tumor-manifest.csv
## OUTPUT: heatmap of patient fidelity in a tumor

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

require(gdata)

manifest_file<-'/Users/chens22/Documents/miRNA/all-tumor-manifest.csv'
base_path<-'/Users/chens22/Documents/miRNA/'
summary_ext<-'.isomir.tsv'
patient_fidelity_ext<-'.patient.fidelity.by.tumor.tsv'
patient_rpm_ext<-'.patient.rpm.by.tumor.tsv'



getPatientIDFromSampleType<-function(tumor, sample_type, manifest_data){

  temp_manifest_data<-manifest_data[which(manifest_data$DISEASE.ABBV==tumor & manifest_data$SAMPLE.TYPE==sample_type),]
  patients<-temp_manifest_data[,c('NAME', 'ID')]
  return(patients)
  
}

removeNameExtension<-function(summary_data){
  unique_patients<-unique(factor(summary_data$SAMPLE))

  patients_with_name<-data.frame(SAMPLE=unique_patients)

  for(i in 1:nrow(patients_with_name)){
    curr_name<-toString(patients_with_name[i, 'SAMPLE'])
    new_name<-strsplit(curr_name, '.', fixed=TRUE)[[1]][1]
    patients_with_name[i,'NAME']<-new_name
  }

  new_summary_data<-merge(summary_data, patients_with_name, by='SAMPLE', all=TRUE, na.rm=TRUE)

  return(new_summary_data)
  
}


getTumorSummaryData<-function(tumor, patient_names_ids){
  
  summary_file<-paste(base_path, tumor, '/summary_files/', tumor, summary_ext, sep='')
  summary_data<-read.table(file=summary_file, header=TRUE, sep='\t')
  
  summary_data<-removeNameExtension(summary_data)
  #write.table(summary_data, file=summary_file, row.names=FALSE, sep='\t')
  
  new_summary_data<-summary_data[which(summary_data$NAME %in% patient_names_ids$NAME),]
  new_summary_data<-merge(new_summary_data, patient_names_ids, by='NAME', all=TRUE)
  
  new_summary_data<-getRPM(new_summary_data)
  
  return(new_summary_data)
  
}


getRPM<-function(summary_data){
  
  summary_data[,'RPM']<-(summary_data$TOTAL.READS/summary_data$TOTAL.READS.IN.SAMPLE)*1000000
  
  return(summary_data)
  
}

getTopExpressedMirna<-function(summary_data, amount){
  
  unique_mirna<-unique(summary_data$MIRNA)
  mirna_rpm<-data.frame(MIRNA=unique_mirna)

  for(i in 1:nrow(mirna_rpm)){
    
    curr_mirna<-toString(mirna_rpm[i,'MIRNA'])
    temp_summary_data<-summary_data[which(summary_data$MIRNA==curr_mirna),]
    rpm<-sum(temp_summary_data$RPM)
    mirna_rpm[i, 'RPM']<-rpm
    
  }
  
  mirna_rpm<-mirna_rpm[complete.cases(mirna_rpm),]
  mirna_rpm<-mirna_rpm[order(-mirna_rpm[,'RPM']),]
  if(!is.na(amount)){
    return(factor(mirna_rpm[1:amount,'MIRNA']))
  }else{
    return(factor(mirna_rpm[,'MIRNA']))
  }
  
  
}

getPatientFidelityByMirna<-function(summ_data, number_of_mirna){
  
  top_expressed_mirna<-getTopExpressedMirna(tumor_summary_data, number_of_mirna)
  
  patient_fidelity<-data.frame(MIRNA=top_expressed_mirna)
  
  patient_id<-unique(summ_data$ID)
  
  for(p in patient_id){
    
    temp_fid_data<-summ_data[which(summ_data$ID==p), c('MIRNA', 'FIDELITY.5P')]
    names(temp_fid_data)[names(temp_fid_data)=='FIDELITY.5P']<-p
    
    patient_fidelity<-merge(patient_fidelity, temp_fid_data, by='MIRNA', all=TRUE, na.rm=TRUE)
  }
  
  patient_fidelity<-patient_fidelity[patient_fidelity$MIRNA%in%top_expressed_mirna,]
  return(patient_fidelity)
  
}

getPatientRPMByMirna<-function(summ_data, number_of_mirna){
  
  top_expressed_mirna<-getTopExpressedMirna(tumor_summary_data, number_of_mirna)
  
  patient_fidelity<-data.frame(MIRNA=top_expressed_mirna)
  
  patient_id<-unique(summ_data$ID)
  
  for(p in patient_id){
    
    temp_fid_data<-summ_data[which(summ_data$ID==p),c('MIRNA', 'RPM')]
    names(temp_fid_data)[names(temp_fid_data)=='RPM']<-p
    
    patient_fidelity<-merge(patient_fidelity, temp_fid_data, by='MIRNA', all=TRUE, na.rm=TRUE)
  }
  
  patient_fidelity<-patient_fidelity[patient_fidelity$MIRNA%in%top_expressed_mirna,]
  return(patient_fidelity)
  
}

plotClusterHeatmap<-function(data, tumor, data_type, sample_type){

  filename<-gsub(' ', '.', paste(base_path, tumor, '/plots/', tumor, '.', sample_type, '.', data_type, '.heatmap.pdf', sep=''))
  
  mirna<-data[,'MIRNA'] #names should be a factor!!
  
  data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
  rownames(data) <- mirna                  # assign row names
  
  palette <- colorRampPalette(c("white", "blue"))(n = 299)
  pdf(filename)
  hcluster <- heatmap.2(data,
                        #cellnote = mat_data,  # same data set for cell labels
                        main = paste(tumor, data_type), # heat map title
                        notecol="black",      # change font color of cell labels to black
                        density.info="none",  # turns off density plot inside color legend
                        trace="none",         # turns off trace lines inside the heat map
                        margins =c(12,9),     # widens margins around plot
                        col=palette,       # use on color palette defined earlier
                        #breaks=col_breaks,    # enable color transition at specified limits
                        dendrogram="both",     # only draw a row dendrogram
                        key.title="NA",
                        key.xlab=data_type,
                        keysize=1,              # size intensity key
                        cexCol=0.2,             # size x labels (samples)
                        Colv="Rowv")            # turn off column clustering
  dev.off()
  
}


writeFile<-function(data, tumor, sample_type, data_type){
  
  filename<-paste(base_path, tumor, '/', tumor, '.', sample_type, '.patient.', data_type, '.by.tumor.tsv', sep='')
  write.table(data, file=filename, row.names=FALSE, sep='\t')
  
}

manifest_data<-read.csv(file=manifest_file, header=TRUE)

tumors<-c('ACC', 'BLCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
          'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
          'SARC', 'SKCM', 'STAD','TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM') 
#BRCA, 
#Redo: all, LAML; solid, STAD KIRC

curr_sample_type<-'SOLID.TISSUE.NORMAL'
num_of_mirna<-NA

for(curr_tumor in tumors){
  filtered_patients<-getPatientIDFromSampleType(curr_tumor, curr_sample_type, manifest_data)
  tumor_summary_data<-getTumorSummaryData(curr_tumor, filtered_patients)
  
  if(is.data.frame(filtered_patients) & nrow(filtered_patients)>0){
    patient_fid_by_tumor<-getPatientFidelityByMirna(tumor_summary_data, num_of_mirna)
    #patient_rpm_by_tumor<-getPatientRPMByMirna(tumor_summary_data, num_of_mirna)
    patient_fid_by_tumor <- patient_fid_by_tumor[,colSums(is.na(patient_fid_by_tumor))<nrow(patient_fid_by_tumor)]
    
    #plotClusterHeatmap(patient_fid_by_tumor, curr_tumor, '5\' FIDELITY', curr_sample_type)
    #plotClusterHeatmap(patient_rpm_by_tumor, curr_tumor, 'RPM')
    patient_fid_by_tumor[,'FIDELITY.AVERAGE']<-rowMeans(patient_fid_by_tumor[,2:ncol(patient_fid_by_tumor)], na.rm=TRUE)
    fidelity_sd<-transform(patient_fid_by_tumor, FIDELITY.SD=apply(patient_fid_by_tumor[,2:(ncol(patient_fid_by_tumor)-1)], 1, sd, na.rm=TRUE)) #subtract one for new average column
    patient_fid_by_tumor<-merge(patient_fid_by_tumor, fidelity_sd[,c('MIRNA', 'FIDELITY.SD')], by='MIRNA', na.rm=TRUE, all=TRUE)
    writeFile(patient_fid_by_tumor, curr_tumor, curr_sample_type, 'fidelity')
    #writeFile(patient_rpm_by_tumor, curr_tumor, curr_sample_type, 'rpm')
  }
} 
