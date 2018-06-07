
## PURPOSE: get longest stem in each sliding window (sub-sequence) of a sequence
## INPUT: sequence  example.sequence.txt
## OUTPUT: table containing longest stem of each sub-sequence   longest-shifting-miRNA-stem.tsv

require('stringr')
require('plyr')

#####SET DIRECTORY FOLDER HERE#####
setwd('/Volumes/2TB (MAC)/Susanna/mRNA/')

#####INPUT FILE WITH SEQUENCE HERE#####
seq_file <- 'example.sequence.txt' 


getSecondaryStructure<-function(sequence){
  input_data<-c(paste('>', sep=''), sequence)
  input_file<-'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-'temp-miRNA-secondary-struct.txt'
  secondary_structure<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  secondary_structure<-strsplit(secondary_structure, ' ', fixed=T)[[1]][1]
  return(toString(secondary_structure))
}

getEnergy<-function(sequence){
  input_data<-c(paste('>', sep=''), sequence)
  input_file<- 'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-'temp-miRNA-secondary-struct.txt'
  output_data<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  energy<-as.numeric(strsplit(output_data, ' ', fixed=T)[[1]][4])
  energy<-gsub(')', '', toString(energy), fixed=T)
  energy<-gsub('(', '', toString(energy), fixed=T)
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][2]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][3]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  return(toString(energy))
}

getMaxStem <- function(sequence){
  max_stem <- ''
  max_start <- 0
  max_end <- 0
  struct <- toString(getSecondaryStructure(sequence))
  
  #start with right bracket, keep on pairing with left as long as it's continuous
  right_brack_pos <- gregexpr(')', struct, fixed=T)[[1]]
  right_brack_start <- as.integer(right_brack_pos[1])
  right_brack_end <- as.integer(right_brack_pos[length(right_brack_pos)])
  
  temp_struct <- substr(toString(struct), 1, right_brack_start)
  left_brack_pos <- gregexpr('(', temp_struct, fixed=T)[[1]]
  left_brack_end <- as.integer(left_brack_pos[length(left_brack_pos)])
  left_brack_start <- as.integer(left_brack_pos[1])
  stem_start <- left_brack_end
  stem_end <- right_brack_start
  while(stem_end <= nchar(sequence)){
    
    if(substr(struct, stem_end, stem_end) == '(' || substr(struct, stem_start, stem_start) == ')' || (stem_start <= left_brack_start & stem_end >= right_brack_end)){
      
      if(substr(struct, stem_end, stem_end) == '('){
        stem_end <- stem_end - 1
      }
      if(substr(struct, stem_start, stem_start) == ')'){
        stem_start <- stem_start + 1
      }
      
      stem <- substr(struct, stem_start, stem_end)
      stem <- toString(stem)
      
      left_paired <- gregexpr('(', stem, fixed=T)[[1]]
      start <- left_paired[1]
      right_paired <- gregexpr(')', stem, fixed=T)[[1]]
      end <- right_paired[length(right_paired)]
      
      stem <- substr(struct, stem_start, stem_end)
      
      left_paired_bases <- str_count(stem, '\\(')
      right_paired_bases <- str_count(stem, '\\)')
      if(left_paired_bases>right_paired_bases){
        diff <- left_paired_bases - right_paired_bases
        start <- gregexpr('(', stem, fixed=T)[[1]]
        start <- start[diff + 1]
        end <- gregexpr( ')', stem, fixed=T)[[1]]
        end <- end[length(end)]
      }
      
      if(right_paired_bases>left_paired_bases){
        diff <- right_paired_bases - left_paired_bases
        start <- gregexpr('(', stem, fixed=T)[[1]]
        start <- start[1]
        end <- gregexpr( ')', stem, fixed=T)[[1]]
        end <- end[length(end) - diff]
      }
      
      
      if(start == -1 || length(start) == 0){ 
        break 
      }
      if( end == -1 || length(end) == 0){ 
        break 
      }
      
      stem <- substr(stem, as.numeric(start), as.numeric(end))
      
      if(nchar(stem) > nchar(max_stem)){
        max_stem <- stem
        max_start <- stem_start
        max_end <- stem_end
      }
      
      
      temp_struct <- substr(struct, stem_end + 1, nchar(struct))
      right_brack_pos <- gregexpr(')', temp_struct, fixed=T)[[1]]
      right_brack_start <- as.integer(right_brack_pos[1])
      temp_struct <- substr(temp_struct, 1, right_brack_start)
      left_brack_pos <- gregexpr('(', temp_struct, fixed=T)[[1]]
      left_brack_end <- as.numeric(left_brack_pos[length(left_brack_pos)])
      stem_start <- left_brack_end + stem_end
      stem_end <- right_brack_start + stem_end
    }
    
    
    if(stem_start > 0){
      stem_start <- stem_start - 1
    }
    if(stem_end < nchar(sequence)){
      stem_end <- stem_end + 1
    }
    
    
  }
  
  return(c(max_stem, max_start, max_end))
}

getGCCount <- function(sequence){
  return(as.numeric(str_count(sequence, 'G')) + as.numeric(str_count(sequence, 'C')))
}


seq_data <- readLines(seq_file)

data <- data.frame(START = integer(),
                   END = integer(),
                   GC.COUNT = integer(),
                   ENERGY = double(),
                   SECONDARY.STRUCTURE = character(),
                   MAX.STEM = character(),
                   SEQUENCE = character(),
                   SEQ.START = integer,
                   stringsAsFactors = F)

window_len <- 100
seq <- seq_data[1]
for(i in 1:(nchar(seq) - window_len)){

  curr_seq <- substr(seq, i, i + window_len - 1)
  max_stem_data <- getMaxStem(curr_seq)
  max_stem <- max_stem_data[1]
  max_stem_start <- as.numeric(max_stem_data[2])
  max_stem_end <- as.numeric(max_stem_data[3])
  gc_count <- getGCCount(curr_seq)
  energy <- getEnergy(curr_seq)


  sub_seq <- substr(curr_seq, max_stem_start, max_stem_end)
  
  # find if UG motif exists
  
  brack_pos_5p <- gregexpr('(', max_stem, fixed=T)[[1]]
  last_brack_5p <- brack_pos_5p[length(brack_pos_5p)]
  temp_seq <- substr(curr_seq, max_stem_start - 1 + last_brack_5p - 43, max_stem_start - 1 + last_brack_5p - 33)
  temp_seq <- toupper(temp_seq)
  ug_found <- gregexpr('TG', temp_seq, fixed=T)[[1]][1]!=-1
  
  # find if CNNC motif exists
  first_brack_3p <- gregexpr(')', max_stem, fixed=T)[[1]][1]
  temp_seq <- substr(curr_seq, max_stem_start - 1 + first_brack_3p + 33, max_stem_start - 1 + first_brack_3p + 43)
  temp_seq <- toupper(temp_seq)
  c_pos <- gregexpr('C', temp_seq, fixed=T)[[1]]
  cnnc_found <- F
  if(length(c_pos)>1){
    past_pos <- c_pos[1]
    for(pos in c_pos[2:length(c_pos)]){
      if(pos-past_pos == 2){
        cnnc_found <- T
      }
    }
  }
  
  #find if UGU motif exists
  temp_seq <- substr(curr_seq, max_stem_start - 1 + last_brack_5p, max_stem_start - 1 + first_brack_3p)
  temp_seq <- toupper(temp_seq)
  ugu_found <- gregexpr('TGT', temp_seq, fixed=T)[[1]][1]!=-1
  



  data <- rbind(data, data.frame(STEM.START = max_stem_start, #within entire sequence
                                 STEM.END = max_stem_end,
                                 GC.COUNT = gc_count,
                                 ENERGY = energy,
                                 SECONDARY.STRUCTURE = toString(getSecondaryStructure(seq)),
                                 MAX.STEM = max_stem,
                                 SEQUENCE = sub_seq,
                                 SEQ.START = i,
                                 UG.FOUND = ug_found,
                                 UGU.FOUND = ugu_found,
                                 CNNC.FOUND = cnnc_found,
                                 stringsAsFactors = F))
  write.table(data, 'longest-shifting-mRNA-stem.tsv', sep='\t', row.names=F)

}
