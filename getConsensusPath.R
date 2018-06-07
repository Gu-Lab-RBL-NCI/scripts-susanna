
## PURPOSE: get path of sequences to the consensus sequence from isomiR data
## INPUT: miRBase data                   miRBase21-master.tsv
##        isomiR expression data  tumor.isomir.expression.tsv
## OUTPUT: table containing list of sequences to consensus
  
getPathToConsensusData<-function(data, mirna_data){
  require('stringdist')
  require('plyr')
  
  getChild<-function(step, root, children, consensus){
    children<-as.character(children)
    closest_dist<-Inf
    closest_child<-c()
    root_consensus<-stringdist(root, consensus)
    
    for(c in children){
      dist<-stringdist(c, root)
      dist_consensus<-stringdist(c, consensus)
      if(dist_consensus<root_consensus && dist<=closest_dist){
        closest_dist<-dist
        updated_closest_child<-c()
        for(x in closest_child){
          curr_dist<-stringdist(x, root)
          if(curr_dist<=closest_dist){
            updated_closest_child<-c(x, updated_closest_child)
          }
        }
        closest_child<-updated_closest_child
        closest_child<-c(c, closest_child)
      }
    }
    
    if(length(closest_child)>1){
      closest_child<-paste(closest_child, collapse=',')
      return(c(toString(closest_child), closest_dist))
    }else if(length(closest_child)==1){
      return(c(toString(closest_child), closest_dist))
    }else{
      return(c(NA, NA))
    }
    
  }
  
  getConsensus<-function(mirna){
    return(toString(unique(mirna_data[which(mirna_data$MIRNA==mirna), 'SEQUENCE'])))
  }
  
  getPathToConsensus<-function(graph, root_seq, consensus){
    
    path<-data.frame(SEQUENCE=NA, CHILD=NA, CHILD.DISTANCE=NA, CONSENSUS.DISTANCE=NA)
    
    graph<-unique(graph)
    
    queue<-c()
    checked<-c()
    distance<-stringdist(root_seq, consensus)
    jumps<-0
    
    queue<-c(root_seq, queue)
    checked<-c(root_seq, checked)
    last_child<-root_seq
    
    while(length(queue)>0){
      curr_seq<-toString(queue[1])
      queue<-queue[-1]
      curr_dist<-stringdist(curr_seq, consensus)
      
      if(curr_seq==consensus){
        return(path)
      }
      
      if(curr_dist<=distance){
        distance<-curr_dist
        children<-as.character(graph[which(graph$SEQUENCE==curr_seq), 'CHILD'])
        if(length(children)>0){
          children<-strsplit(children, ',', fixed=TRUE)[[1]]
          for(child in children){
            child_dist<-stringdist(child, last_child)
            if(!(child%in%checked)){
              path<-rbind(path, data.frame(SEQUENCE=curr_seq, 
                                           CHILD=child, 
                                           CHILD.DISTANCE=child_dist, 
                                           CONSENSUS.DISTANCE=distance))
              checked<-c(toString(child), checked)
              queue<-c(toString(child), queue)
              jumps<-jumps+1
              last_child<-child
            }
          }
        }else{
          #print(as.character(curr_seq))
        }
      }
    }
    return(path)
  }
  
  #get child and child distance for each isomir from table
  unique_mirna<-unique(data$MIRNA)
  unique_mirna<-unique_mirna[!is.na(unique_mirna)]
  count<-0
  for(m in unique_mirna){
    count<-count+1
    temp_data<-data[which(data$MIRNA==m),]
    temp_data<-temp_data[!is.na(temp_data$SEQUENCE),]
    sequences<-temp_data$SEQUENCE
    sequences<-sequences[!is.na(sequences)]
    consensus<-getConsensus(m)
    
    if(nrow(temp_data)>1){
      for(i in 1:nrow(temp_data)){
        curr_seq<-temp_data[i, 'SEQUENCE']
        closest_seq<-getChild(1, curr_seq, sequences[-which(sequences%in%curr_seq)], consensus)
        temp_data[i, 'CHILD']<-closest_seq[1]
        temp_data[i, 'CHILD.DISTANCE']<-closest_seq[2]
      }
    }
    data<-data[-which(data$MIRNA==m),]
    data<-rbind.fill(data, temp_data)
  }
  return(data)
  #get number of steps to consensus and maximum jump for each 
  
  for(i in 1:nrow(data)){
    dist_to_consensus<-data[i, 'DISTANCE']
    if(dist_to_consensus>0){
      seq<-data[i, 'SEQUENCE']
      consensus<-getConsensus(toString(data[i, 'MIRNA']))
      path<-getPathToConsensus(data[,c('SEQUENCE', 'CHILD', 'CHILD.DISTANCE')], seq, consensus)
      steps<-nrow(path)-1
      if(length(path$CHILD.DISTANCE)>=1){
        max_jump<-max(path$CHILD.DISTANCE, na.rm=TRUE)
      }
      data[i, 'STEPS.TO.CONSENSUS']<-steps
      data[i, 'MAX.JUMP']<-max_jump
      last_dist<-as.numeric(path[nrow(path), 'CONSENSUS.DISTANCE'])
      if(!is.na(last_dist)){  
        if(last_dist==1){
          data[i, 'PATH.FOUND']<-TRUE
        }else{
          data[i, 'PATH.FOUND']<-FALSE
        }
      }
    }
  }
  return(data)
}
