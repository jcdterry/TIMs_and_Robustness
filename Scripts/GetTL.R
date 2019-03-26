GetTL <- function(web){
  # Baesd on  Petchey Github. Follows Levine 1980
  # Finds Trophic Levels of a species in a web 
  # takes predation matrix with consumers in columns
  
  # Can get a bit over-excited about loops and rapidly lead to v.high trophic levels
  # Perhaps wpould be good to cap the TL?
  
  diag(web)<-0 # Remove cannibal links
  
  tweb <- t(web)
  sp<- length(tweb[,1])
  rs <- rowSums(tweb)
  for(i in 1:sp) tweb[i,tweb[i,]==1] = 1/rs[i]    ## make the rows add to one
  nb.TL <- try(solve(diag(sp) - tweb), T)
  if(class(nb.TL)=="try-error")return(rep(NA, sp))
  if(class(nb.TL)!="try-error")return(rowSums(nb.TL))
}

