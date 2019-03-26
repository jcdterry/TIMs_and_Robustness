CalcExpTIMs <- function(Model, TIM_Density){
  
  NumPred <- rowSums(Model$web) 
  NumPrey <- colSums(Model$web) 
  BiConn  <- NumPrey+NumPred
  
  Model$TIMlist -> TIMlist
  
  TIMlist %>%
    mutate(FinalSign =  slope_ijk*Flip ) -> TIMlist
  
  
  AllExpNum <- c()
  HuntedExpNum <-c()
  mean_BenefTIMsfromHtoEx_list<-c()
  mean_DetriTIMsfromHtoEx_list<-c()
  mean_BenefTIMsfromAnytoEx_list<-c()
  mean_DetriTIMsfromAnytoEx_list<-c()
  
  for( Sp in 1:Model$StartingNum){
    postExtEnd <- Model$postExtEnd_List[[Sp]]
    ExtinctSp <- which(postExtEnd ==0)
    
    ExtinctSp <-ExtinctSp[ExtinctSp!=Sp] # don't look for self-links
    
    
    BenefTIMsfromHtoEx_list<-c()
    DetriTIMsfromHtoEx_list<-c()
    BenefTIMsfromAnytoEx_list<-c()
    DetriTIMsfromAnytoEx_list<-c()
    
    if(postExtEnd[1] != 'DNF' & length(ExtinctSp)>0){
      TIMsfromHuntedToExtinct_list<-c()
      
     
      for( ExSp in ExtinctSp){
        
        ## Beneficial for ExtinctSp (Ex = resource and negative TIM, or Ex = Consumer and positive TIM)
        TIMlist %>%filter(k == Sp) %>% 
          filter((i == ExSp&   sign(FinalSign)==-1 )  | (j == ExSp &   sign(FinalSign)==1)  ) %>% nrow-> BenefTIMsfromHtoEx
        
        ## Detrimental for ExtinctSp (Ex = resource and positive TIM, or Ex = Consumer and negative TIM)
        TIMlist %>%filter(k == Sp) %>% 
          filter((i == ExSp&   sign(FinalSign)==1 )  | (j == ExSp &   sign(FinalSign)==-1)  ) %>% nrow-> DetriTIMsfromHtoEx
        ## Beneficial for ExtinctSp (Ex = resource and negative TIM, or Ex = Consumer and positive TIM)
        TIMlist %>%
          filter((i == ExSp&   sign(FinalSign)==-1 )  | (j == ExSp &   sign(FinalSign)==1)  ) %>% nrow-> BenefTIMsfromAnytoEx
        
        ## Detrimental for ExtinctSp (Ex = resource and positive TIM, or Ex = Consumer and negative TIM)
        TIMlist %>%
          filter((i == ExSp&   sign(FinalSign)==1 )  | (j == ExSp &   sign(FinalSign)==-1)  )%>% nrow -> DetriTIMsfromAnytoEx
        
        BenefTIMsfromHtoEx_list<-c(BenefTIMsfromHtoEx_list,BenefTIMsfromHtoEx)
        DetriTIMsfromHtoEx_list<-c(DetriTIMsfromHtoEx_list,DetriTIMsfromHtoEx)
        BenefTIMsfromAnytoEx_list<-c(BenefTIMsfromAnytoEx_list,BenefTIMsfromAnytoEx)
        DetriTIMsfromAnytoEx_list<-c(DetriTIMsfromAnytoEx_list,DetriTIMsfromAnytoEx)
        
      }
       mean_BenefTIMsfromHtoEx_list  <-c( mean_BenefTIMsfromHtoEx_list  ,mean(BenefTIMsfromHtoEx_list) )
    mean_DetriTIMsfromHtoEx_list  <-c( mean_DetriTIMsfromHtoEx_list  ,mean(DetriTIMsfromHtoEx_list) )
    mean_BenefTIMsfromAnytoEx_list<-c( mean_BenefTIMsfromAnytoEx_list,mean(BenefTIMsfromAnytoEx_list) )
    mean_DetriTIMsfromAnytoEx_list<-c( mean_DetriTIMsfromAnytoEx_list,mean(DetriTIMsfromAnytoEx_list) )
    
    AllExpNum <-c(AllExpNum,  mean(BiConn[ExtinctSp] * TIM_Density * ((Model$S)-2) *0.5)   ) #  *0.5  is for each interaction giving either benefical or negative
    HuntedExpNum <-c(HuntedExpNum, mean(BiConn[ExtinctSp] * TIM_Density *1  *0.5  ))
      
    }
    
   
    
    
    
  }

  return( data.frame( AllExpNum,
                                 HuntedExpNum,
                                 mean_BenefTIMsfromHtoEx_list,
                                 mean_DetriTIMsfromHtoEx_list,
                                 mean_BenefTIMsfromAnytoEx_list,
                                 mean_DetriTIMsfromAnytoEx_list, 
                      'ID'= Model$ID) )
}