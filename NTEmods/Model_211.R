
             
             ModelIIPars <- list(Aij= matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.3476,0.3476,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.09486,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07382,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.03146,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02659,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.003958,0.003958,0.003958,0.003958,0,0.003958,0.003958,0,0,0,0,0,0,0,0,0,0.01779,0.01779,0.01779,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.03963,0.03963,0.03963,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00164,0.00164,0.00164,0.00164,0.00164,0.00164,0.00164,0.00164,0.00164,0.00164,0,0.00164,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.006034,0,0.006034,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00664,0.00664,0.00664,0,0),nrow=20,ncol=20),
             Es= matrix(c(0.068,0.147,0.091,0.13,0.149,0.064,0.076,0.145,0.073,0.055,0.082,0.087,0.144,0.076,0.082,0.094,0.062,0.132,0.079,0.06,0.122,0.085,0.127,0.088,0.083,0.057,0.135,0.141,0.082,0.086,0.071,0.051,0.088,0.068,0.095,0.109,0.054,0.06,0.127,0.111,0.098,0.095,0.068,0.119,0.123,0.103,0.13,0.109,0.082,0.121,0.088,0.076,0.114,0.082,0.113,0.123,0.129,0.085,0.148,0.118,0.134,0.134,0.111,0.058,0.074,0.109,0.088,0.088,0.077,0.102,0.063,0.131,0.128,0.132,0.15,0.075,0.112,0.127,0.138,0.095,0.115,0.128,0.146,0.064,0.063,0.072,0.081,0.147,0.083,0.099,0.15,0.107,0.12,0.121,0.079,0.111,0.144,0.054,0.068,0.07,0.108,0.079,0.089,0.131,0.108,0.126,0.13,0.113,0.124,0.09,0.145,0.103,0.069,0.119,0.098,0.109,0.117,0.08,0.055,0.115,0.093,0.122,0.139,0.097,0.079,0.118,0.147,0.082,0.147,0.066,0.144,0.127,0.134,0.096,0.15,0.088,0.093,0.098,0.051,0.054,0.086,0.103,0.124,0.054,0.07,0.133,0.116,0.117,0.091,0.056,0.087,0.078,0.127,0.089,0.145,0.079,0.061,0.109,0.147,0.112,0.135,0.066,0.052,0.108,0.137,0.065,0.088,0.091,0.098,0.138,0.095,0.08,0.058,0.134,0.143,0.129,0.113,0.15,0.141,0.145,0.119,0.107,0.079,0.114,0.075,0.105,0.08,0.144,0.129,0.097,0.063,0.054,0.126,0.133,0.116,0.127,0.138,0.051,0.101,0.067,0.066,0.085,0.1,0.115,0.067,0.086,0.083,0.056,0.093,0.076,0.135,0.127,0.1,0.061,0.09,0.083,0.078,0.124,0.106,0.089,0.084,0.122,0.055,0.068,0.092,0.122,0.134,0.101,0.135,0.15,0.089,0.102,0.11,0.109,0.059,0.082,0.089,0.078,0.063,0.064,0.103,0.141,0.086,0.058,0.079,0.109,0.095,0.125,0.051,0.087,0.064,0.05,0.088,0.066,0.056,0.051,0.123,0.088,0.059,0.074,0.102,0.134,0.1,0.058,0.082,0.115,0.056,0.074,0.1,0.082,0.117,0.117,0.092,0.084,0.134,0.059,0.138,0.053,0.099,0.079,0.057,0.082,0.075,0.124,0.138,0.083,0.137,0.12,0.053,0.108,0.055,0.109,0.14,0.079,0.087,0.108,0.096,0.083,0.149,0.089,0.066,0.096,0.075,0.129,0.08,0.074,0.104,0.059,0.053,0.08,0.115,0.059,0.055,0.09,0.05,0.075,0.141,0.109,0.063,0.1,0.096,0.104,0.148,0.097,0.147,0.082,0.148,0.11,0.064,0.145,0.08,0.087,0.138,0.069,0.14,0.073,0.139,0.112,0.134,0.1,0.108,0.108,0.055,0.084,0.116,0.065,0.079,0.111,0.127,0.08,0.134,0.06,0.124,0.072,0.06,0.08,0.127,0.119,0.056,0.126,0.079,0.108,0.145,0.065,0.134,0.093,0.074,0.104,0.07,0.09,0.098,0.074,0.109,0.105,0.138,0.097,0.108,0.104,0.107,0.087,0.057,0.054,0.064,0.141,0.133,0.099,0.078,0.116,0.082,0.106,0.093,0.062,0.12,0.14,0.066,0.081,0.106,0.121,0.113,0.137),nrow=20,ncol=20),
             r = c(1,1,-0.06952,1,1,1,-0.01897,-0.01476,-0.022,-0.002659,-0.003166,-0.01423,1,1,1,-0.02373,-0.002459,1,-0.004219,-0.003976),
             K = c(17.09,18.24,730.9,6.081,10.18,2.604,246.9,463.2,233.5,294.7,590.8,102.8,8.833,7.196,5.52,107.3,660.9,18.75,171.2,191),
             Starts =  c(3.001,3.203,2.372,5.258,8.195,2.199,0.4926,1.473,1.897,0.05817,0.8439,0.4332,8.801,7.194,5.518,3.219,0.1895,18.48,1.343,0.8852))
             
ModelIIPars$NTEonResConst <- t(t(ModelIIPars$Aij)*ModelIIPars$Starts)
ModelIIPars$NTEonConConst <- ModelIIPars$Es*ModelIIPars$Aij*ModelIIPars$Starts


             LVModel <- function(t, Bs, pars=NULL){
   TIMlist<- pars$TIMs
             ## ID = 211 
             
             if(any(is.nan(Bs))){Bs[is.nan(Bs)]<- 0}
             if(any(Bs<1e-10)){Bs[Bs<1e-10]<- 0}
             
             LogRatioVect   <-  log10( Bs[TIMlist$k] / TIMlist$Start) 
             LogRatioVect[LogRatioVect< -10]<- -10  # Put a floor on change to stop infinity issues
             
             LogmuVect = TIMlist$Flip * (TIMlist$Range_ijk *  exp(-(TIMlist$gomp_b) * exp(- (TIMlist$gomp_c)* LogRatioVect  )) - TIMlist$gomp_f0)
             LogmuVect[is.na(LogmuVect)]<-0
             Logmu_ijk<- array(0, c(20,20,20))
             Logmu_ijk[ as.matrix(TIMlist[,1:3])  ] <- LogmuVect
             musub1<-  rowSums((10^(Logmu_ijk))-1, dims=2)
             
             NTEViaBeingRes <- rowSums(ModelIIPars$NTEonResConst * musub1)
             NTEViaBeingCon <- colSums(ModelIIPars$NTEonConConst * musub1) 


             Yij <-  t(pars$Aij)*Bs  # losses rate per loser
             Zij <-  pars$Es*pars$Aij*Bs # gain rate per gainer
             
             dBs<- Bs* ( pars$r - pars$Mort  -(Bs/pars$K)   -colSums(Yij)+   colSums(Zij)  - NTEViaBeingRes + NTEViaBeingCon )
             
             dBs[is.nan(Bs)]<-0
             dBs[Bs<1e-10]<-0
             
             return(list(as.vector(dBs)))
             }
             