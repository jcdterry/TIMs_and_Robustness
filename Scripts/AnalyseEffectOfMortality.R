AnalyseEffectOfMortality<-function(Model){
  
  
  ExtinctionTable<-  data.frame('ModelID'=rep(Model$ID,Model$StartingNum),
                                'Finished' = NA,'SpeciesHunted'= NA,
                                'HuntedSp_TL'=NA ,  'HuntedSp_OldTL'=NA ,
                                'Extinctions'=NA, 
                                'NumFuncExtinction'=NA, 'NumPopExplosion'=NA, 'HuntedSpExtinct'=NA, 'meanchange'=NA,
                                'HuntedNumPrey' =NA,
                                'HuntedNumPred'  =NA,
                                'HuntedBiConn' =NA, 
                                'MeanExtinctNumPrey' = NA,
                                'MeanExtinctNumPred' = NA,
                                'MeanExtinctBiConn' = NA, 
                                'MeanExtinctTL'= NA,
                                'MeanNumTIMsfromHtoEx' =NA,
                                'MeanBenefTIMsfromHtoEx' =NA,
                                'MeanDetriTIMsfromHtoEx'=NA,
                                'MeanNumTIMsfromAnytoEx'=NA,
                                'MeanBenefTIMsfromAnytoEx' =NA,
                                'MeanDetriTIMsfromAnytoEx'=NA,
                                'MeanTrophDistHtoEx' = NA,
                                'NumTIMsFromHunted' = NA,
                                'NumPosTIMsFromHunted' = NA,
                                'NumNegTIMsFromHunted' = NA)
  
  Model$web %>%
    graph_from_adjacency_matrix() %>%
    distances(mode = 'all') -> TrophicDistancesMatrix
  
  NewTL<-GetTL(Model$web)
  
  NumPred <- rowSums(Model$web) 
  NumPrey <- colSums(Model$web) 
  BiConn  <- NumPrey+NumPred
  
  
  Model$TIMlist -> TIMlist
  
  TIMlist %>%
    mutate(FinalSign =  slope_ijk*Flip ) -> TIMlist
  
  
  # BiDirectional Trophic Connectance
  TrophConnMat <- (Model$web  + t(Model$web))>0
  ExtinctionTable$TrophConn <- sum( TrophConnMat  )/ (Model$S*(Model$S-1 ))
  
  # Non-Trophic (TIM) Connectance
  TIMconnMat <- matrix(0, Model$S, Model$S)
  TIMconnMat[matrix(c(TIMlist$i, TIMlist$k), ncol=2)] <-1
  TIMconnMat[matrix(c(TIMlist$j, TIMlist$k), ncol=2)] <-1
  ExtinctionTable$TIMConn <- sum(TIMconnMat)/ (Model$S*(Model$S-1 ))
  
  # Full Dynamic Connectance
  ExtinctionTable$DynConn <-   sum((TrophConnMat+TIMconnMat)>0) / (Model$S*(Model$S-1 ))
  
  
  
  for( Sp in 1:Model$StartingNum){
    ExtinctionTable$SpeciesHunted[Sp] <- Sp
    
    postExtEnd <- Model$postExtEnd_List[[Sp]]
    
    ExtinctSp <- which(postExtEnd ==0)
    
    if(postExtEnd[1] == 'DNF'){
      ExtinctionTable$Finished[Sp] <- FALSE
    }else{
      ExtinctionTable$Finished[Sp] <- TRUE
      ExtinctionTable$Extinctions[Sp]    <- Model$StartingNum - sum(postExtEnd>0)
      ExtinctionTable$NumFuncExtinction[Sp] <- sum((postExtEnd/Model$StartPoints)<0.1) # Population less than 1/100th of start
      ExtinctionTable$NumPopExplosion[Sp] <- sum((postExtEnd/Model$StartPoints)>10) # Population more than 100x Start 
      ExtinctionTable$HuntedSpExtinct[Sp]<-    Sp %in%which(postExtEnd==0)
      ExtinctionTable$HuntedSp_OldTL[Sp]<-    Model$TL[Sp]
      ExtinctionTable$HuntedSp_TL[Sp]<-    NewTL[Sp]
      ExtinctionTable$meanchange[Sp] <- mean(abs(log10(postExtEnd[postExtEnd>0]/Model$StartPoints[postExtEnd>0])))
      
      
      ExtinctionTable$HuntedNumPrey[Sp] <-NumPrey[Sp]
      ExtinctionTable$HuntedNumPred[Sp] <- NumPred[Sp]
      ExtinctionTable$HuntedBiConn [Sp] <- BiConn[Sp]
      
      ## Number of TIMs from Hunted
      TIMlist %>% filter(k == Sp) %>% nrow  -> NumTIMsFromHunted
      ## Number of Strengthening TIMs from Hunted
      TIMlist %>% filter(k == Sp) %>% filter(sign(FinalSign)==1) %>% nrow  -> NumPosTIMsFromHunted
      ## Number of Weakening TIMs from Hunted
      TIMlist %>% filter(k == Sp) %>% filter(sign(FinalSign)==-1)%>% nrow  -> NumNegTIMsFromHunted
      
      ExtinctionTable$NumTIMsFromHunted[Sp] <- NumTIMsFromHunted
      ExtinctionTable$NumPosTIMsFromHunted[Sp] <- NumPosTIMsFromHunted
      ExtinctionTable$NumNegTIMsFromHunted[Sp] <- NumNegTIMsFromHunted
      
      # Excluding the Hunted Species as an extinct species for links between
      ExtinctSp <-ExtinctSp[ExtinctSp!=Sp]
      
      if(length(ExtinctSp)>0){
        TIMsfromHuntedToExtinct_list<-c()
        BenefTIMsfromHtoEx_list<-c()
        DetriTIMsfromHtoEx_list<-c()
        TIMsfromAnyToExtinct_list<-c()
        BenefTIMsfromAnytoEx_list<-c()
        DetriTIMsfromAnytoEx_list<-c()
        
        
        for( ExSp in ExtinctSp){
          
          ## Mean TIMs from hunted onto each extinct
          TIMlist %>%filter(k == Sp) %>%
            filter(i == ExSp| j == ExSp) %>% nrow-> TIMsfromHuntedToExtinct   
          
          ## Beneficial for ExtinctSp (Ex = resource and negative TIM, or Ex = Consumer and positive TIM)
          TIMlist %>%filter(k == Sp) %>% 
            filter((i == ExSp&   sign(FinalSign)==-1 )  | (j == ExSp &   sign(FinalSign)==1)  ) %>% nrow-> BenefTIMsfromHtoEx
          
          ## Detrimental for ExtinctSp (Ex = resource and positive TIM, or Ex = Consumer and negative TIM)
          TIMlist %>%filter(k == Sp) %>% 
            filter((i == ExSp&   sign(FinalSign)==1 )  | (j == ExSp &   sign(FinalSign)==-1)  ) %>% nrow-> DetriTIMsfromHtoEx
          
          ## TIMs from any onto extinct
          TIMlist %>%
            filter(i == ExSp| j == ExSp ) %>% nrow-> TIMsfromAnyToExtinct   
          
          ## Beneficial for ExtinctSp (Ex = resource and negative TIM, or Ex = Consumer and positive TIM)
          TIMlist %>%
            filter((i == ExSp&   sign(FinalSign)==-1 )  | (j == ExSp &   sign(FinalSign)==1)  ) %>% nrow-> BenefTIMsfromAnytoEx
          
          ## Detrimental for ExtinctSp (Ex = resource and positive TIM, or Ex = Consumer and negative TIM)
          TIMlist %>%
            filter((i == ExSp&   sign(FinalSign)==1 )  | (j == ExSp &   sign(FinalSign)==-1)  )%>% nrow -> DetriTIMsfromAnytoEx
          
          TIMsfromHuntedToExtinct_list<-c(TIMsfromHuntedToExtinct_list,TIMsfromHuntedToExtinct)
          BenefTIMsfromHtoEx_list<-c(BenefTIMsfromHtoEx_list,BenefTIMsfromHtoEx)
          DetriTIMsfromHtoEx_list<-c(DetriTIMsfromHtoEx_list,DetriTIMsfromHtoEx)
          TIMsfromAnyToExtinct_list<-c(TIMsfromAnyToExtinct_list,TIMsfromAnyToExtinct)
          BenefTIMsfromAnytoEx_list<-c(BenefTIMsfromAnytoEx_list,BenefTIMsfromAnytoEx)
          DetriTIMsfromAnytoEx_list<-c(DetriTIMsfromAnytoEx_list,DetriTIMsfromAnytoEx)
          }
        ExtinctionTable$MeanNumTIMsfromHtoEx[Sp] <-  mean(TIMsfromHuntedToExtinct_list)
        ExtinctionTable$MeanBenefTIMsfromHtoEx[Sp] <- mean(BenefTIMsfromHtoEx_list)
        ExtinctionTable$MeanDetriTIMsfromHtoEx[Sp]<-  mean(DetriTIMsfromHtoEx_list)
        ExtinctionTable$MeanNumTIMsfromAnytoEx[Sp] <-  mean(TIMsfromAnyToExtinct_list)
        ExtinctionTable$MeanBenefTIMsfromAnytoEx[Sp] <- mean(BenefTIMsfromAnytoEx_list)
        ExtinctionTable$MeanDetriTIMsfromAnytoEx[Sp]<-  mean(DetriTIMsfromAnytoEx_list)
        
        ## Trophic Distance between hunted and extinct
        ExtinctionTable$MeanTrophDistHtoEx[Sp] <-   mean(TrophicDistancesMatrix[ExtinctSp,Sp])
        
        ## Connectance of Extinct Species
        ExtinctionTable$MeanExtinctNumPrey[Sp] <-   mean(NumPrey[ExtinctSp])
        ExtinctionTable$MeanExtinctNumPred[Sp] <-   mean(NumPred[ExtinctSp])
        ExtinctionTable$MeanExtinctBiConn[Sp] <-   mean(BiConn[ExtinctSp])
        
        ExtinctionTable$MeanExtinctTL[Sp] <- mean(NewTL[ExtinctSp])
        
        
      }
    }
  }
  
  return(ExtinctionTable)
}





