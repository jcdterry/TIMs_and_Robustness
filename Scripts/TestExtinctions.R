
TestExtinctions <- function(Starts, web, TL, Pars, S, Threshold= 10^-4, FUNC){
  
  result<- NULL
  NEED_TO_RERUN <- TRUE
  FAIL <- FALSE
  counter=0
  
  
  output<-ode(Starts, times = seq(0, 10, l =200), FUNC,parms = Pars)
  
  Bs <-output[200,2:(S+1)]
  mean<-colMeans(output[,-1])
  
  while(NEED_TO_RERUN) { 
    if(FAIL){  return(rep(NA,S))}  
    
    ## Main Running Loop
    for(i in 1:10){
      Bs[is.nan(Bs)]<- 0
      Bs[Bs<Threshold]<-0
      output<-ode(Bs, times = seq(0, 100, l =200), FUNC,parms = Pars)
      meandiff<-abs((colMeans(output[,-1]) - mean)/mean)
      if(all(meandiff[is.finite(meandiff)] <0.01)){
        result<- output[200,2:(S+1)]
        break
      }else{
        Bs<- output[200,2:(S+1)]
        Bs[is.nan(Bs)]<- 0
        mean<- colMeans(output[,-1])
      }
    }
    if(all(is.na(Bs))){
      NEED_TO_RERUN <- FALSE# Give Up if all Na
    }else{
      
      NEED_TO_RERUN <- FALSE
      if(any(meandiff[is.finite(meandiff)]>0.01)){
        NEED_TO_RERUN <- TRUE
        }
      ## Check all consumers have food
      persist <- Bs>0 # update persisters 
      trimweb <- web
      trimweb[,!persist]<-0
      trimweb[!persist,]<-0  
      HaveFood<-colSums(trimweb)>0 
      
      if(any( TL>1&!HaveFood & persist )){ ## Kill any consumers that don't have food and keep running
        Bs[TL>1&!HaveFood] <-0
        NEED_TO_RERUN <- TRUE
      }
    } 
    counter= counter+1
    if(counter == 50){ 
     # print('counter limit reached')
     # print(tail(output))
      return('DNF')}
  }
 # print(FUNC(1, as.vector(unname(Bs)), pars =  list('TIMs'=TIMs)))
return(as.vector(unname(Bs)))
}
