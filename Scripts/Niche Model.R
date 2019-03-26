niche.model<-function(S,C){ 
  # S = the number of species, C = target connectance
  # from Klecka Plos one
  
  C.web=2 # C.web = actual connectance, set absurd to initiate  while loop
  TL=TRUE # initiatior
  ProdFrac=0 # Producer fraction, set to iniate while loop
  components=FALSE # False Initiator of number of components
  
  # Keep running so that:
  # C.web should be +- 3% of the target value of C
  # TL should be determinable
  # The graph should only have one component (ie all linked together)
  while((C.web>1.03*C)|(C.web<0.97*C)|any(is.na(TL))|components!=1 ){
    con<-rep(0,times=S) # vector containing 0 for disconnected species
    
    # while some elements of con are zero, try again until  no disconnected species
    while(length(con[con==0])>0|ProdFrac<0.2){
      beta<-(1-2*C)/(2*C)
      sp.niche<-sort(runif(n=S,min=0,max=1)) # generate niche values
      sp.ri<-c(1-(1-runif(S))^(1/beta))*sp.niche # niche widths
      sp.ri[1]<-0 # smallest species always primary producer - zero niche width
      for (n in 1:S){
        while (sp.ri[n]/2>sp.niche[n]){# s r_i must satisfy  condition
          sp.ri[n]<-(1-(1-runif(1))^(1/beta))}
      } 
      #  find centre of feeding niches:
      sp.cen<-sapply(1:S,function (x) runif(1,min=sp.ri[x]/2,max=sp.niche[x])) 
      
      web<-matrix(0,nrow=S,ncol=S) # initiate the food web matrix
      
      # Add the interactions:
      for (j in 1:S){
        for (i in 1:S){
          if ((sp.niche[i]<(sp.cen[j]+sp.ri[j]/2))&(sp.niche[i]>(sp.cen[j]-sp.ri[j]/2))){
            web[i,j]<-1} 
        }
      } 
      #             Purely cannibalistic species are counted as disconnected:
      web.zero.diagonal<-web*(matrix(1,nrow=S,ncol=S)-diag(S))
      #             Disconnected species have a value of con=0:
      con<-colSums(web.zero.diagonal)+rowSums(web.zero.diagonal) 
      #             Determine fraction of species that are producers L
      ProdFrac<- length(which(colSums(web.zero.diagonal)==0))/S 
    }
    C.web=sum(web)/S^2 

    
    # Make sure only one component to the web (no fragments)
    components<-(web %>%graph_from_adjacency_matrix(mode='undirected') %>%count_components) 
    
  }
  
  diag(web)<-0
  
  return(web)
}
