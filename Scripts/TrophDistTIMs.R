TrophDistTIMs <- function(Model){
  
  Model$web %>%
    graph_from_adjacency_matrix() %>%
    distances(mode = 'all') -> TrophicDistancesMatrix
  
  TrophDistList<-c()
  
  for( Sp in 1:Model$StartingNum){
    postExtEnd <- Model$postExtEnd_List[[Sp]]
    ExtinctSp <- which(postExtEnd ==0)
    
    ExtinctSp <-ExtinctSp[ExtinctSp!=Sp] # don't look for self-links
    
    if(postExtEnd[1] != 'DNF' & length(ExtinctSp)>0){
      TrophDistList<-c(TrophDistList,TrophicDistancesMatrix[Sp,ExtinctSp ])
    }
  }
  ExtinctionTable<-  data.frame( TrophDistList) 
  return(ExtinctionTable)
}