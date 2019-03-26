
AddTIMsMort <- function(Model, TIM_Density, TIM_Slopes, TIM_Inflexion,
                        TIM_Range, LimitTIMsPerInteraction, FolderName = 'TIMmods'){
  
  Model <- CutModelDown(Model)
  Model$TIMlist<-  GeneratingTIMs(Model,  TIM_Density, TIM_Slopes, TIM_Inflexion,TIM_Range, LimitTIMsPerInteraction)
  Model$FolderString <- paste0(FolderName,'/Model_')
  
  S<- Model$S

  cat(paste0('

  ModelIIPars <- list(Aij= matrix(c(',paste(as.vector(signif(Model$Aij,4)), collapse=','), '),nrow=',S,',ncol=',S,'),
                    Es= matrix(c(',paste(as.vector(signif(Model$Es,4)), collapse=','),'),nrow=',S,',ncol=',S,'),
                    r = c(', paste(signif(Model$r,4), collapse=','), '),
                    K = c(', paste(signif(Model$K,4), collapse=','), '))
  
  LVModel <- function(t, Bs, pars=NULL){
   TIMlist<- pars$TIMs
   ## ID = ',Model$ID,' 
   
   if(any(is.nan(Bs))){Bs[is.nan(Bs)]<- 0}
   if(any(Bs<', Model$min,')){Bs[Bs<',Model$min,']<- 0}

   LogRatioVect   <-  log10( Bs[TIMlist$k] / TIMlist$Start) 
   LogRatioVect[LogRatioVect< -10]<- -10  # Put a floor on change to stop infinity issues

   LogmuVect = TIMlist$Flip * (TIMlist$Range_ijk *  exp(-(TIMlist$gomp_b) * exp(- (TIMlist$gomp_c)* LogRatioVect  )) - TIMlist$gomp_f0)
   LogmuVect[is.na(LogmuVect)]<-0
   Logmu_ijk<- array(0, c(',S,',',S,',',S,'))
   Logmu_ijk[ as.matrix(TIMlist[,1:3])  ] <- LogmuVect
   mu_ij<- 10^(rowSums(Logmu_ijk, dims = 2))

   Amodij <- mu_ij*pars$Aij

   Yij <-  t(Amodij)*Bs  # losses rate per loser
   Zij <-  pars$Es*Amodij*Bs # gain rate per gainer
             
   dBs<- Bs* ( pars$r - pars$Mort  -(Bs/pars$K)   -colSums(Yij)+   colSums(Zij))
             
   dBs[is.nan(Bs)]<-0
   dBs[Bs<',Model$min,']<-0
   
   return(list(as.vector(dBs)))
  }
             '), file= paste0(Model$FolderString,Model$ID,'.R'))
  
  return(Model)
}
