TestConstMortality<- function(Model, AddedMortality=1){
  
  source(file = paste0( Model$FolderString,Model$ID,'.R')) # Adds LVModelTIM() and ModelIIPars object to environment 
  ModelIIPars$TIMs=Model$TIMlist
  
    ModelIIPars$Mort<- rep(AddedMortality, Model$StartingNum)  
  
  postExtEnd<-'DNF'
  try({
    postExtEnd<- TestExtinctions( Starts = Model$StartPoints,
                                  S = Model$StartingNum,
                                  web =  Model$web,
                                  TL = Model$TL,
                                  Pars = ModelIIPars,
                                  FUNC = LVModel)
  }, silent = TRUE)
  
  
  ExtinctionTable<-  data.frame('ModelID'=Model$ID, 
                                'AddedMortality' = AddedMortality,
                                'Finished' = NA,
                                'Extinctions'=NA,
                                'NumFuncExtinction'  =NA,   
                                'NumPopExplosion' =NA,   
                                'meanchange'  =NA,   
                                'MeanNumTIMsfromAnytoEx' =NA,   
                                'MeanBenefTIMsfromAnytoEx'= NA,
                                'MeanDetriTIMsfromAnytoEx' = NA) 
  
  
  
  Model$TIMlist  %>%
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
  
  ExtinctSp <- which(postExtEnd ==0)
  
  if(postExtEnd[1] == 'DNF'){
    ExtinctionTable$Finished <- FALSE
  }else{
    ExtinctionTable$Finished <- TRUE
    ExtinctionTable$Extinctions    <- Model$StartingNum - sum(postExtEnd>0)
    ExtinctionTable$NumFuncExtinction <- sum((postExtEnd/Model$StartPoints)<0.01) # Population less than 1/100th of start
    ExtinctionTable$NumPopExplosion <- sum((postExtEnd/Model$StartPoints)>100) # Population more than 100x Start 
    ExtinctionTable$meanchange <- mean(abs(log10(postExtEnd[postExtEnd>0]/Model$StartPoints[postExtEnd>0])))

    if(length(ExtinctSp)>0){
      TIMsfromAnyToExtinct_list<-c()
      BenefTIMsfromAnytoEx_list<-c()
      DetriTIMsfromAnytoEx_list<-c()
      for( ExSp in ExtinctSp){
        
        ## TIMs from any onto extinct
        TIMlist %>%
          filter(i == ExSp| j == ExSp ) %>% nrow-> TIMsfromAnyToExtinct   
        
        ## Beneficial for ExtinctSp (Ex = resource and negative TIM, or Ex = Consumer and positive TIM)
        TIMlist %>%
          filter((i == ExSp&   sign(FinalSign)==-1 )  | (j == ExSp &   sign(FinalSign)==1)  ) %>% nrow-> BenefTIMsfromAnytoEx
        
        ## Detrimental for ExtinctSp (Ex = resource and positive TIM, or Ex = Consumer and negative TIM)
        TIMlist %>%
          filter((i == ExSp&   sign(FinalSign)==1 )  | (j == ExSp &   sign(FinalSign)==-1)  )%>% nrow -> DetriTIMsfromAnytoEx
        
        TIMsfromAnyToExtinct_list<-c(TIMsfromAnyToExtinct_list,TIMsfromAnyToExtinct)
        BenefTIMsfromAnytoEx_list<-c(BenefTIMsfromAnytoEx_list,BenefTIMsfromAnytoEx)
        DetriTIMsfromAnytoEx_list<-c(DetriTIMsfromAnytoEx_list,DetriTIMsfromAnytoEx)
      }
      ExtinctionTable$MeanNumTIMsfromAnytoEx <-  mean(TIMsfromAnyToExtinct_list)
      ExtinctionTable$MeanBenefTIMsfromAnytoEx <- mean(BenefTIMsfromAnytoEx_list)
      ExtinctionTable$MeanDetriTIMsfromAnytoEx<-  mean(DetriTIMsfromAnytoEx_list)
    }
  }
  
  return(ExtinctionTable)
  
}