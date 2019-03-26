MakeModelCore <- function(    ID, S=20, C=0.15, min= 10^-10, k =10,CarryCapMult=2){
  
  
  
  set.seed(ID)
  makeweb=TRUE
  while(makeweb){
    web<-niche.model(S = S, C = C)
    TL<-  GetTL(web)
    if(all(!is.na(TL))){
      makeweb <- FALSE
      if(sum(TL==1)>10){makeweb <- TRUE}} # Make sure less than 10 producers
  }
  intTL<-  round(TL+0.0001) # Adding 0.0001 deals with issue of R rounding .5 towards even  digit...
  body.masses<-BodyMasses(intTL)
  
  Ws <-signif((1/colSums(web)),3) # Calculating attack rate reduction due to prey diversity
  Ws[is.infinite(Ws)]<-0 # set non-consumer attack rates to zero
  
  #Mi <-  matrix(body.masses,S,S, byrow = FALSE) 
  Mj <-  matrix(body.masses,S,S, byrow =TRUE) 
  wj <-  matrix(Ws,S,S        , byrow =TRUE) 
  
  Aij <-   web * wj * k * (Mj^(-0.25))
  Es <- matrix(round(runif(S^2, min = 0.05, max=0.15), 3), nrow = S, ncol=S)
  
  r<- K<-rep(0, S)
  r[TL==1]<-1
  r[TL>1]<-  -0.1*(body.masses^(-0.25))[TL>1]
  K[TL==1]<-  CarryCapMult*runif(S, 1,10)[TL==1] # Define K as for producers
  K[TL> 1]<-  10^(runif(S, 2, 3)[TL>1])  # Density Dependence of consumers
  
  
  return(list('ID'=ID, 'web'=web, 'Ws'=Ws, 'Mj'=Mj, 'TL'=intTL, 'Aij'=Aij, 'Es'=Es, 'r'=r, 'K'=K, 'S'=S, 'C'=C, 'min'=min, 'BMs' = body.masses))
  
}