
jobs=seq(5,8,1)
setting<-function(job){
  if(job ==1) {
  signal.var = 1e-2
  corr_signal=1e-2
  pi=0.2
} else if (job == 2) {
  signal.var = 1e-2
  corr_signal=2e-1
  pi=0.2
} else if (job == 3) {
  signal.var = 1e-2
  corr_signal=3e-1
  pi=0.2
} else if (job == 4) {
  signal.var = 1e-2
  corr_signal=4e-1
  pi=0.2
} else if (job == 5) {
  signal.var = 1e-2
  corr_signal=5e-2
  pi=0.01
} else if (job == 6) {
  signal.var = 1e-2
  corr_signal=1e-2
  pi=0.01
} else if (job ==7) {
  signal.var = 1e-2
  corr_signal=1e-2
  pi=0.2
} else if (job ==8) {
  signal.var = 1e-2
  corr_signal=2e-1
  pi=0.01
}
  return(c(signal.var,corr_signal,pi))
  }
require(TwoSampleMR)
require(data.table)
require(dplyr)
require(mr.divw)
require(IDPmisc)
#sample size
reprecessing<-function(result,iter)
{
  library(IDPmisc)
  IVW_TRUE=NULL
  IVW_OVB=NULL
  DIVW_TRUE=NULL
  DIVW_OVB=NULL
  IVW_TRUE_sd=NULL
  IVW_OVB_sd=NULL
  DIVW_TRUE_sd=NULL
  DIVW_OVB_sd=NULL
  IVW_TRUE_iv=NULL
  IVW_OVB_iv=NULL
  DIVW_TRUE_iv=NULL
  DIVW_OVB_iv=NULL
  IVW_TRUE_f=NULL
  IVW_OVB_f=NULL
  DIVW_TRUE_f=NULL
  DIVW_OVB_f=NULL
  for(i in(1:(iter-1)))
  {
    IVW_TRUE[i]=result[[i]]$IVW_true[[1]]$beta
    DIVW_TRUE[i]=result[[i]]$DIVW_true$beta.hat
    IVW_OVB[i]=result[[i]]$IVW_ovb[[1]]$beta
    DIVW_OVB[i]=result[[i]]$DIVW_ovb$beta.hat
    IVW_TRUE_sd[i]=result[[i]]$IVW_true[[1]]$sd
    IVW_OVB_sd[i]=result[[i]]$IVW_ovb[[1]]$sd
    DIVW_TRUE_sd[i]=result[[i]]$DIVW_true$beta.se
    DIVW_OVB_sd[i]=result[[i]]$DIVW_ovb$beta.se
    IVW_TRUE_iv[i]=result[[i]]$IVW_true[[1]]$numIV_sel
    IVW_OVB_iv[i]=result[[i]]$IVW_ovb[[1]]$numIV_sel
    DIVW_TRUE_iv[i]=result[[i]]$DIVW_true$n.IV
    DIVW_OVB_iv[i]=result[[i]]$DIVW_ovb$n.IV
    IVW_TRUE_f[i]=result[[i]]$IVW_true[[1]]$f
    IVW_OVB_f[i]=result[[i]]$IVW_ovb[[1]]$f
    DIVW_TRUE_f[i]=result[[i]]$DIVW_true$f
    DIVW_OVB_f[i]=result[[i]]$DIVW_ovb$f
  }
  IVW_true=cbind( mean(na.omit(IVW_TRUE)), sd(na.omit(IVW_TRUE)),mean(na.omit(IVW_TRUE_sd)),  mean(na.omit(IVW_TRUE_iv)),
                  mean(na.omit( IVW_TRUE_f))  )
  IVW_OVB=cbind( mean( na.omit(IVW_OVB)),sd( na.omit(IVW_OVB)),mean( na.omit(IVW_OVB_sd)), mean(na.omit(IVW_OVB_iv)),
                 mean(na.omit(IVW_OVB_f)))
DIVW_true=cbind( mean(na.omit( DIVW_TRUE)),sd( na.omit(DIVW_TRUE)),mean( na.omit(DIVW_TRUE_sd)),mean(na.omit(DIVW_TRUE_iv)),mean(na.omit( DIVW_TRUE_f)))

DIVW_ovb=cbind( mean(na.omit(DIVW_OVB)),sd(na.omit(DIVW_OVB)), mean(na.omit(DIVW_OVB_sd)), mean( na.omit(DIVW_OVB_iv)), mean(na.omit(DIVW_OVB_f)))
result=rbind(IVW_true,IVW_OVB,DIVW_true,DIVW_ovb)
result=data.frame(result)
names(result)=c("beta","MCSD","SE","nIV","f")
return(result)
}
IVW_correct<-function(gamma.exp,gamma.out,se.exp,se.out,pthr=5e-8)
{
  C_sel = qnorm(pthr/2,lower.tail = FALSE)
  ind_filter = which(abs(gamma.exp / se.exp) >= C_sel)
  numIV_sel = length(ind_filter)
  gamma.exp_sel = gamma.exp[ind_filter]
  gamma.out_sel = gamma.out[ind_filter]
  se.exp_sel = se.exp[ind_filter]
  se.out_sel = se.out[ind_filter]
  
  gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
  
  beta = sum(gamma.exp_sel*gamma.out_sel *(1/se.out_sel^2)) / sum(gamma.exp_sel^2/se.out_sel^2)
  sd = TwoSampleMR::mr_ivw(gamma.exp_sel,gamma.out_sel,se.exp_sel,se.out_sel)$se
  p = pnorm(abs(beta/sd), lower.tail = F) * 2
  f=sum(gamma.exp_sel^2 / se.exp_sel^2)/ length(se.exp_sel) - 1
  
  prop_psmall = round(sum(gamma.exp.pval_sel < 5e-8 & gamma.exp.pval_sel > 5e-10)/numIV_sel * 100, 2)
  
  output = data.frame(beta, sd, p,f,numIV_sel,  prop_psmall)
  
  return(list(output))
} 
for (job in jobs){
signal_val=  setting(job)[1]
corr_signal=  setting(job)[2]
pi= setting(job)[3]
n=2e4
#gammaX
nsnp=50000
correlation=rep(0,nsnp)
gamma_x=rep(0,nsnp)
gamma_y=rep(0,nsnp)

gamma_x_ovb=rep(0,nsnp)
gamma_y_ovb=rep(0,nsnp)


for ( i in (1:nsnp))
{
  rand_number=runif(1,0,1)
  if(rand_number<0.2)
  {
    correlation[i]=rnorm(1,0,sd=corr_signal)
    if(abs(correlation[i])>0.5)
      correlation[i]=runif(1,0,0.5)
  }
  else
    correlation[i]=0
}
#correlation_marginal=matrix(rep(0,nsnp*nsnp),nsnp,nsnp)
#diag(correlation_marginal)=1/n
#for(i in (1:nsnp))
#{
 # for(j in (1:nsnp))
  #{
   # if(i!=j)
    #{correlation_marginal[i,j]=-1/n*correlation[i]*correlation[j]
    #correlation_marginal[j,i]=-1/n*correlation[i]*correlation[j]
    #}
  #}
#}

#nsnp2=nsnp+1
#correlation_joint=matrix(rep(0,nsnp2*nsnp2),nsnp2,nsnp2)
#diag(correlation_joint)=1
#for(i in (2:nsnp2))
#{
 # correlation_joint[1,i]=correlation[i-1]
#  correlation_joint[i,1]=correlation[i-1]
#}
#correlation_joint=1/n*solve(correlation_joint)[2:nsnp2,2:nsnp2]
for ( i in (1:nsnp))
{
  rand_number=runif(1,0,1)
  if(rand_number<pi)
    gamma_x[i]=rnorm(1,0,sd=signal_val)
  else
   gamma_x[i]=0
}
beta=0.2
gamma_y=beta*gamma_x
#ovb

for ( i in (1:nsnp))
{
  print(i)
  gamma_temp=gamma_x[-i]
  correlation_temp=correlation[-i]
  gamma_x_ovb[i]=gamma_x[i]-correlation[i]*t(gamma_temp)%*%correlation_temp
}
gamma_y_ovb=beta*gamma_x_ovb
iter=1
repeat_times=500
sd_x=rep(1/sqrt(n),nsnp)
sd_y=rep(1/sqrt(n),nsnp)
Result=list()
while(iter<repeat_times)
{

  gamma_joint_x=rep(0,nsnp)
  gamma_joint_y=rep(0,nsnp)
  gamma_marginal_x=rep(0,nsnp)
  gamma_marginal_y=rep(0,nsnp)
  for(i in(1:nsnp))
  {  gamma_joint_x[i]=gamma_x[i]+rnorm(1,0,sd_x[i])
  gamma_joint_y[i]=gamma_y[i]+rnorm(1,0,sd_y[i])
  gamma_marginal_x[i]=gamma_x_ovb[i]+rnorm(1,0,sd_x[i])
  gamma_marginal_y[i]=gamma_y_ovb[i]+rnorm(1,0,sd_y[i])
 
  }# gamma_joint_x=mvrnorm(1,gamma_x,correlation_joint)
  #gamma_joint_y=mvrnorm(1,gamma_y,correlation_joint)
  #gamma_marginal_x=mvrnorm(1,gamma_x_ovb,correlation_marginal)
  #gamma_marginal_y=mvrnorm(1,gamma_y_ovb,correlation_marginal)
 IVW_true= IVW_correct(gamma_joint_x,gamma_joint_y,sd_x,sd_y)
 IVW_ovb= IVW_correct(gamma_marginal_x,gamma_marginal_y,sd_x,sd_y)
 divw_true=  mr.divw(gamma_joint_x,gamma_joint_y,sd_x,sd_y)
  divw_ovb=mr.divw(gamma_marginal_x,gamma_marginal_y,sd_x,sd_y)
  divw_f_joint=sum(gamma_joint_x^2/sd_x^2)/length(gamma_joint_x)-1
  divw_f_ovb=sum(gamma_marginal_x^2/sd_x^2)/length(gamma_marginal_x)-1
  divw_ovb$f=divw_f_ovb
  divw_true$f=divw_f_joint
 Result[[iter]]=list(IVW_true=IVW_true,IVW_ovb=IVW_ovb,DIVW_true=divw_true,
                  DIVW_ovb=divw_ovb)
 iter=iter+1
 print(iter)
  
}
result=reprecessing(Result,iter)
file_path=paste0("C:\\Users\\25110\\Desktop\\50000",pi,signal_val,corr_signal,".RData")
save.image(file_path)
file_path2=paste0("C:\\Users\\25110\\Desktop\\50000",pi,signal_val,corr_signal,".csv")
write.csv(result,file_path2)

}