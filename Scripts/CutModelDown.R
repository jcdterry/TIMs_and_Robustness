### Cut Down To Size

CutModelDown <- function( Model){

Starts <-which(Model$Start>0)
Model$S<- sum(Model$Start>0)

Model$Aij <- Model$Aij[Starts,];
Model$Aij <- Model$Aij[,Starts]
Model$Es <- Model$Es[Starts,]; Model$Es <- Model$Es[,Starts]
Model$web<-Model$web[Starts,];  Model$web<-Model$web[,Starts]
Model$r<- Model$r[Starts]
Model$K<- Model$K[Starts]   
Model$TL<- Model$TL[Starts]
Model$StartPoints<- Model$Start[Starts]

return(Model)
}