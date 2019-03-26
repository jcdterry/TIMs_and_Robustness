NumExtinctCalc<-function(Model){
  
  ExtinctionTable<-  data.frame('ModelID'=rep(Model$ID,Model$StartingNum),
                                'Finished' = NA,'SpeciesHunted'= NA,
                                'Extinctions'=NA, 
                                'NumFuncExtinction'=NA, 
                                'NumPopExplosion'=NA,
                                'HuntedSpExtinct'=NA)
  
  for( Sp in 1:Model$StartingNum){
    ExtinctionTable$SpeciesHunted[Sp] <- Sp
    postExtEnd <- Model$postExtEnd_List[[Sp]]
    ExtinctSp <- which(postExtEnd ==0)
    if(postExtEnd[1] == 'DNF'){
      ExtinctionTable$Finished[Sp] <- FALSE
    }else{
      ExtinctionTable$Finished[Sp] <- TRUE
      ExtinctionTable$Extinctions[Sp]    <- Model$StartingNum - sum(postExtEnd>0)
      ExtinctionTable$NumFuncExtinction[Sp] <- sum((postExtEnd/Model$StartPoints)<0.1) # Population <1/10th of start
      ExtinctionTable$NumPopExplosion[Sp] <- sum((postExtEnd/Model$StartPoints)>10) # Population >10x Start 
      ExtinctionTable$HuntedSpExtinct[Sp]<-    Sp %in%which(postExtEnd==0)
    }
  }
  return(ExtinctionTable)
}
