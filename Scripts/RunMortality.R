RunMortality<- function(Model, AddedMortality=1){
  
  source(file = paste0( Model$FolderString,Model$ID,'.R')) # Adds LVModelTIM() and ModelIIPars object to environment 
  ModelIIPars$TIMs <- Model$TIMlist
  
  postExtEnd_List <- list()
  
  for( Sp in 1:Model$StartingNum){
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
    postExtEnd_List[[Sp]]<- postExtEnd
    cat('.')
  }
  cat('\n')
  
  
  Model$postExtEnd_List <- postExtEnd_List
  Model$AddedMortality <- AddedMortality
  
  return(Model)
  
}