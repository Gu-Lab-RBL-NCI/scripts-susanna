## PURPOSE: return list of positions where motif is found in vector of sequences
## INPUT: list of sequences and motif
## OUTPUT: position

locateIUPACMotif <- function(sequences, motif){
  
  iupac_key <- vector(mode='list', length=15)
  names(iupac_key) <- c('A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 
                        'M', 'B', 'D', 'H', 'V', 'N')
  iupac_key['A'] <- 'A'
  iupac_key['C'] <- 'C'
  iupac_key['G'] <- 'G'
  iupac_key['T'] <- 'T'
  iupac_key['R'] <- 'A,G'
  iupac_key['Y'] <- 'C,T'
  iupac_key['S'] <- 'C,G'
  iupac_key['W'] <- 'A,T'
  iupac_key['K'] <- 'T,G'
  iupac_key['M'] <- 'A,C'
  iupac_key['B'] <- 'T,C,G'
  iupac_key['D'] <- 'T,A,G'
  iupac_key['H'] <- 'T,C,A'
  iupac_key['V'] <- 'A,C,G'
  iupac_key['N'] <- 'T,C,G,A'
  

  motif <- toString(motif)
  motif_combs <- c(motif)
  for(i in 1:nchar(motif)){
    curr_base <- substr(motif, i, i)
    bases <- iupac_key[curr_base]
    
    if(toString(bases[1])=='NULL'){
      bases <- c(curr_base)
    }
    bases <- strsplit(toString(bases), ',')[[1]]
    
    new_motif_combs <- c()
    for(m in motif_combs){
      
      for(b in bases){
        new_motif <- paste0(substr(m, 1, i - 1), b, substr(m, i + 1, nchar(m)))
        new_motif_combs <- c(new_motif_combs, new_motif)
      }
      
    }
    motif_combs <- new_motif_combs
  }
  
  
  seq_found <- c()
  
  for(m in motif_combs){
    
    seq_found <- c(seq_found, grep(m, sequences))
    
  }
  
  seq_found <- unique(seq_found)
  
  return(seq_found)
  
}
