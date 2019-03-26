#### Generating TIMs
# TIM density = the fraction of possible TIMs that exist:

GeneratingTIMs <- function(Model,  TIM_Density, TIM_Slopes, TIM_Inflexion,TIM_Range, LimitTIMsPerInteraction){
  
  
  Model$TIM_Density_Setting <- TIM_Density
  Model$TIM_Slope_Setting <- TIM_Slopes
  Model$TIM_Inflexion_Setting <- TIM_Inflexion
  Model$TIM_Range_Setting <-TIM_Range
  Model$TIMperInt_Setting <- LimitTIMsPerInteraction
  
  ID <- Model$ID
  web<- Model$web
  S<- Model$S
  StartPoints <- Model$StartPoints
  
  set.seed(ID)
  
  if(LimitTIMsPerInteraction){
    TIMlist<-as.data.frame(which(web==1, arr.ind = TRUE))
    colnames(TIMlist)  <-c('i','j')
    TIMlist$k<-NA
    NumTIMsWanted <-  floor(TIM_Density * sum(web) *(S-2))
    
    if( NumTIMsWanted >  nrow(TIMlist) ){
      warning('Requested TIM Density too high') # num mods must be smaller than
      NumTIMsWanted <-  nrow(TIMlist)
    }
    
    TIMlist <- TIMlist[sample(1:nrow(TIMlist), NumTIMsWanted, replace=FALSE),]
    for(xx in 1:nrow(TIMlist)){
      TIMlist$k[xx]<-sample(   (1:S)[-c(TIMlist$i[xx], TIMlist$j[xx])] , 1   )
    }
  }else{
    Present<-array(rep(web,S) * (runif(S^3)<TIM_Density), c(S,S,S))
    TIMlist<-as.data.frame(which(Present==1, arr.ind = TRUE))
    colnames(TIMlist)  <-c('i','j','k')
  }
  set.seed(ID)
  
  TIMlist$slope_ijk   <- runif(nrow(TIMlist), min = -TIM_Slopes, max=TIM_Slopes)
  TIMlist$Inflect_ijk <- runif(nrow(TIMlist), min = 0-TIM_Inflexion, max=TIM_Inflexion)
  TIMlist$Range_ijk   <- runif(nrow(TIMlist), 0.1,TIM_Range)
  TIMlist$Flip        <- sample(c(-1,1), nrow(TIMlist), replace=TRUE) # this terms evens out positive and negative TIMs by flipping half
  
  
  ## Inferred Params
  TIMlist$gomp_c<-  TIMlist$slope_ijk*   exp(1) /  TIMlist$Range_ijk
  TIMlist$gomp_b <-  exp( TIMlist$Inflect_ijk *  TIMlist$gomp_c )
  TIMlist$gomp_f0 <-  TIMlist$Range_ijk * exp(-  TIMlist$gomp_b)
  TIMlist$Start <- StartPoints[TIMlist$k]
  
  return(TIMlist)
}
