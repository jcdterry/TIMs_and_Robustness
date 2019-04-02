
TestInitialResponse<- function(Model, AddedMortality=1){
  
  source(file = paste0( Model$FolderString,Model$ID,'.R')) # Adds LVModel() and ModelIIPars object to environment 
  ModelIIPars$TIMs <- Model$TIMlist
  
  DirectionChange <- list()
  
  for( Sp in 1:Model$StartingNum){
    
    ModelIIPars$Mort<- rep(0, Model$StartingNum)  
    ModelIIPars$Mort[Sp] <- 1
    FirstNudge <- Model$StartPoints
    FirstNudge[Sp] <- FirstNudge[Sp] *0.99
    responses<-unlist(LVModel(t=1, Bs = FirstNudge,pars = ModelIIPars))
    SignResponse <-  sign(round(responses, digits = 6))
    SignResponse <-SignResponse[-Sp]
    
    DirectionChange[[Sp]]<- SignResponse
  }
#  cat('.')
  return(unlist(DirectionChange))
}