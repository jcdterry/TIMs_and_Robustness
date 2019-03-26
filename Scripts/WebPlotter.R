WebPlotter<-function(web, persist=NULL, TLs=NULL, remove=FALSE){
  
  web %>% graph_from_adjacency_matrix()->gr
  
  lay<-matrix(nrow=dim(web)[1],ncol=2) # create a matrix with one column as runif, the other as trophic level
  lay[,1]<-sample((1:(dim(web)[1]))*4000, dim(web)[1] )
  lay[,2]<-GetTL(web)

  if(is.null(persist)){
  plot.igraph(gr,layout=lay)
  }else{
    
    if(remove){
    lay[!persist,1]<-0
    lay[!persist,2]<-0
    }
    plot.igraph(gr,layout=lay, vertex.color = persist+1)
  }
}
