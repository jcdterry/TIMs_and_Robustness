CreateModel <- function(ID=1, S=25, C=0.2, min=10^-10, smallestsize=10, k=1,CarryCapMult=2, MakePlots=FALSE){
  
  model<-Type1OrigModelMaker(ID,S,C, min, k, CarryCapMult=CarryCapMult)
  source(file = paste0('BasicOrigMods/Model',ID,'.R')) # This adds LVModel and BasicPars to Environment
  
  Out <- Safely_find_eqm( S=S, min = min,smallestsize=smallestsize,
                          web= model$web, TL=model$TL, MakePlots=MakePlots,
                          parms = BasicPars)$result
  if(is.null(Out)){
    print('Fail')
    Success<-FALSE
  }else{
    eqm<- unname(unlist(Out))
    if(all(is.na(eqm))){
      Success<-FALSE
    }else{
      Success=TRUE
      
      ## Perfecting Eqm with root finding
      RootFunc<-function(ShortBs, parms){
        LongBs <- rep(0, S)
        LongBs[eqm!=0]<-ShortBs
        dBs_long <- unlist(LVModel(1, LongBs,parms))
        dBs_short<- dBs_long[eqm!=0]
        return(dBs_short)
      }
      eqm_perf<- eqm
      eqm_perf[eqm!=0]<- multiroot(RootFunc,unname(eqm[eqm!=0]),
                                   parms =BasicPars  ,
                                   positive = TRUE)$root
      if(any(  eqm_perf[eqm!=0] < min)){
        warning('RootFinder is killing off species!! ID:,',ID,',. Counting as a Fail')
        Success <- FALSE
      }
    }
  }
  
  model$smallestsize= smallestsize
  model$modeltype = 'HollingType1'

  
  if(Success){
    model$Success= Success
    model$OrigEqm = eqm
    model$Start =  eqm_perf
    model$StartingNum = sum(eqm>0)
  }else{
    model$Success= FALSE
    model$OrigEqm = NA
    model$Start =  rep(NA, S)
    model$StartingNum = NA
  }
  return(model)
}

