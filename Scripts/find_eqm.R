
find_eqm <- function(S,min=10^-10, smallestsize=10, web, TL, MakePlots=FALSE, parms){
  result<- NULL
  NEED_TO_RERUN <- TRUE
  FAIL <- FALSE
  
  output<-ode(rep(10, S), times = seq(0, 10, l =200), LVModel,parms = parms)
  mean<-rep(1,S)
  Bs <-output[200,2:(S+1)]
  
  while(NEED_TO_RERUN) { 
    if(FAIL){  return(rep(NA,S))}  
    
    for(i in 1:10){
      Bs[is.nan(Bs)]<- 0
      Bs[Bs<min]<-0
      output<-ode(Bs, times = seq(0, 100, l =200), LVModel,parms = parms)
      meandiff<-abs((colMeans(output[,-1]) - mean)/mean)
      if(all(meandiff[is.finite(meandiff)] <0.01)){
        
        result<- output[200,2:(S+1)]
        break
      }else{
        mean<- colMeans(output[,-1])
        Bs <-output[200,2:(S+1)]
      }
    }
    result<- output[200,2:(S+1)]
    result[is.nan(Bs)]<- 0
    
    if(all(is.na(result)) |sum(result>0) < smallestsize ){
      FAIL<-TRUE
      NEED_TO_RERUN <- FALSE
    }else{
      
      NEED_TO_RERUN <- FALSE
      # Checks changes are small enough (relative)
      if(any(meandiff[is.finite(meandiff)]>0.01)){
        NEED_TO_RERUN <- TRUE}

      
      persist <- result>0 # update persisters 
      trimweb <- web
      trimweb[,!persist]<-0
      trimweb[!persist,]<-0  # remove non-persisters and disconnected
      HaveFood<-colSums(trimweb)>0 
      
      if(any( TL>1&!HaveFood  & persist )){ ## Kill   any consumers that don't have food and keep running
        result[TL>1&!HaveFood] <-0
        NEED_TO_RERUN <- TRUE
        Bs<-result
      }
      ## Check connected:
      persist <- result>0
      
      trimweb <- web
      trimweb[,!persist]<-0
      trimweb[!persist,]<-0 # Start with original web, remove non-persisters
      
      cc<- components(graph_from_adjacency_matrix(trimweb)) # Get Components
      
      if(MakePlots){
        # plot(graph_from_adjacency_matrix(trimweb), vertex.color = persist+1)
        WebPlotter(web, persist)
        WebPlotter(trimweb, persist, remove=TRUE)
      }
      
      
      InMainComponent <- cc$membership == which.max(cc$csize)
      
      if(any( !InMainComponent  & persist )){
        result[!InMainComponent] <-0
        NEED_TO_RERUN <- TRUE
        Bs<-result
        if( sum(result>0, na.rm = TRUE) < smallestsize){
          FAIL<-TRUE
          NEED_TO_RERUN <- FALSE
        }  
      }
    }
  }
  
  if(FAIL){
    return(as.data.frame(rep(NA,S)))
  }else{
    return(as.data.frame(result))} 
}

Safely_find_eqm<- safely(find_eqm)