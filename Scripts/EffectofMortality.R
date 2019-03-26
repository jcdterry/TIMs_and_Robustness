## NB this functions role has now been split into 'Run Mortaility' and 'Analyse Effect of Mortality'

EffectofMortality<-function(Model, AddedMortality=0){
  
  
  ExtinctionTable<-  data.frame('ModelID'=rep(Model$ID,Model$StartingNum),
                                'Finished' = NA,'SpeciesHunted'= NA,'HuntedSp_TL'=NA ,'Extinctions'=NA, 
                                'NumFuncExtinction'=NA, 'NumPopExplosion'=NA, 'HuntedSpExtinct'=NA, 'meanchange'=NA)
  
  source(file = paste0( Model$FolderString,Model$ID,'.R')) # Adds LVModelTIM() and ModelIIPars object to environment 
  ModelIIPars$TIMs=Model$TIMlist
  
  for( Sp in 1:Model$StartingNum){
    
    ExtinctionTable$SpeciesHunted[Sp] <- Sp
    ModelIIPars$Mort<- rep(0, Model$StartingNum)  
    ModelIIPars$Mort[Sp] <- AddedMortality
    
    postExtEnd<-'DNF'
    try({
      postExtEnd<- TestExtinctions( Starts = Model$StartPoints,
                                    S = Model$StartingNum,
                                    web =  Model$web,
                                    TL = Model$TL,
                                    Pars = ModelIIPars,
                                    FUNC = LVModel)
    }, silent = TRUE)
    
    if(postExtEnd[1] == 'DNF'){
      ExtinctionTable$Finished[Sp] <- FALSE
    }else{
      ExtinctionTable$Finished[Sp] <- TRUE
      ExtinctionTable$Extinctions[Sp]    <- Model$StartingNum - sum(postExtEnd>0)
      ExtinctionTable$NumFuncExtinction[Sp] <- sum((postExtEnd/Model$StartPoints)<0.1) # Population less than 1/100th of start
      ExtinctionTable$NumPopExplosion[Sp] <- sum((postExtEnd/Model$StartPoints)>10) # Population more than 100x Start 
      ExtinctionTable$HuntedSpExtinct[Sp]<-    Sp %in%which(postExtEnd==0)
      ExtinctionTable$HuntedSp_TL[Sp]<-    Model$TL[Sp]
      ExtinctionTable$meanchange[Sp] <- mean(abs(log10(postExtEnd[postExtEnd>0]/Model$StartPoints[postExtEnd>0])))
    }
    cat('.')
  }
  cat('\n')
  return(ExtinctionTable)
}
