
Type1OrigModelMaker <- function(ID, S=20, C=0.15, min= 10^-10, k =1, CarryCapMult=1){
  
  Model <- MakeModelCore(ID, S, C, min, k,CarryCapMult)
  
  cat(paste0('

BasicPars <- list(Aij= matrix(c(',paste(as.vector(signif(Model$Aij,4)), collapse=','), '),nrow=',S,',ncol=',S,'),
             Es= matrix(c(',paste(as.vector(signif(Model$Es,4)), collapse=','),'),nrow=',S,',ncol=',S,'),
             r = c(', paste(signif(Model$r,4), collapse=','), '),
             K = c(', paste(signif(Model$K,4), collapse=','), '))

LVModel <- function(t, Bs, pars=NULL){
            ## ID = ',ID,' 

             if(any(is.nan(Bs))){Bs[is.nan(Bs)]<- 0}
             if(any(Bs<', min,')){Bs[Bs<',min,']<- 0}

            
             Yij <-  t(pars$Aij)*Bs  # losses rate per loser
             Zij <-  pars$Es*pars$Aij*Bs # gain rate per gainer

             dBs<- Bs* ( pars$r  -(Bs/pars$K)   -colSums(Yij)+   colSums(Zij))
             
             dBs[is.nan(Bs)]<-0
             dBs[Bs<',min,']<-0
             
             return(list(as.vector(dBs)))
}
'), file= paste0('BasicOrigMods/Model',ID,'.R'))
  return(Model)
}