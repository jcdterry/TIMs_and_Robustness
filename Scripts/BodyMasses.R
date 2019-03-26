BodyMasses<- function(TL){
  
  # following Iles and Novak (2016) Am.Nat. to assign body masses 
  #  based on trophic level using ratio distributions
  
  body.masses<-rep(NA,length(TL))
  body.masses[TL==1]<-0 # NB All body masses are on log10 scale, until later!
  body.masses[TL==2]<- rnorm(sum(TL==2),  0.65, 0.5) 

  # Remove the possiblity that consumers are smaller than herbivores
  while(any(  body.masses[TL==2]<0)){
    body.masses[TL==2][body.masses[TL==2]<0] <- rnorm(sum(body.masses[TL==2]<0),  0.65, 0.5) 
  }
  
  for(k in 1:length(TL)){
    if(TL[k] >2) body.masses[k]<- 0.65 + sum(rnorm(TL[k]-2,2.73,0.5/sqrt(TL[k]-2))) 
  }
  body.masses<- 10^body.masses # Back to linear scale
  
  return(body.masses)  
}

