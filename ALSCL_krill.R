# stock assessment of krill
mat_func<-function(L50,L95,length)
{
  b1=log(0.95/0.05)/(L95-L50)
  bo = -L50*b1
  logit_pt =  bo+b1*length
  matp= exp(logit_pt)/(1+exp(logit_pt))
  return(matp)
}
# input data
setwd("E:\\ALSCL new")
cl<-read.table("new alscl.csv",header=T,sep=",")
wgt<-read.table("wgt-2mm-all.csv",header=T,sep=",")
mat<-read.table("mat-2mm-all.csv",header=T,sep=",")
logN_at_len<-as.matrix(log(cl[,3:22]+1e-5))
na_matrix<-matrix(1,nrow=24,ncol=20)
na_matrix[which(cl[,3:22]==0)]=0
len_mid=seq(13,59,2)
len_lower=c(-Inf,seq(14,58,2))
len_upper=seq(14,60,2)
len_border=seq(14,58,2)
log_q<-log(mat_func(30,38,len_upper))
age=c(1:7)
Y=20
A=7
L=24
weight=wgt[,3:22]
mat=mat[,3:22]
M=0
growth_step=1

# prepare data
tmb.data<-list(
  logN_at_len=logN_at_len,
  na_matrix=na_matrix,
  log_q=log_q,
  len_mid=len_mid,
  len_lower=len_lower,
  len_upper=len_upper,
  len_border=len_border,
  age=age,
  Y=Y,
  A=A,
  L=L,
  weight=as.matrix(weight),
  mat=as.matrix(mat),
  M=M,
  growth_step=growth_step
)

parameters = list(
  log_init_Z =log(0.5),
  log_sigma_log_N0 = log(1),
  
  mean_log_R = 5,
  log_sigma_log_R = log(1),
  logit_log_R = log(0.01/0.99),
  
  mean_log_F = log(0.3),
  log_sigma_log_F = log(1),
  logit_log_F_y = log(0.75/0.25),
  logit_log_F_l = log(0.75/0.25),
  
  log_vbk = log(0.4),
  log_Linf = log(60),
  log_t0 = log(1/60),
  log_cv_len = log(0.3),
  log_cv_grow = log(0.3),
  
  log_sigma_index = log(0.1),
  
  # random effects
  dev_log_R = rep(0,tmb.data$Y),
  dev_log_F = array(0,c(tmb.data$L,tmb.data$Y)),
  dev_log_N0 = rep(0,(tmb.data$A-1))
)

parameters.L = list(
  # fixed effects
  log_init_Z = log(0.01),
  log_sigma_log_N0 =-Inf,
  
  mean_log_R = log(10),
  log_sigma_log_R = log(0.01),
  logit_log_R = -30,
  
  mean_log_F = log(0.01),
  #log_sigma_log_F = log(0.01),
  logit_log_F_y = -20,
  logit_log_F_l = -10,
  
  log_vbk = log(0.01),
  log_Linf = log(40),
  #log_t0 = -20,
  log_cv_len = log(0.01),
  log_cv_grow = log(0.01),
  
  log_sigma_index = log(0.01)
)

parameters.U = list(
  # fixed effects
  log_init_Z = log(10),
  log_sigma_log_N0 = log(10),
  
  mean_log_R = 20,
  log_sigma_log_R = log(10),
  logit_log_R = 20,
  
  mean_log_F = log(2),
  #log_sigma_log_F = log(10),
  logit_log_F_y = 20,
  logit_log_F_l = 10,
  
  log_vbk = log(1),
  log_Linf = log(70),
  #log_t0 = 0,
  log_cv_len = log(1),
  log_cv_grow = log(1),
  
  log_sigma_index = log(1.5)
)

lower=unlist(parameters.L)
upper=unlist(parameters.U)

map = list(
  #log_sigma_log_N0 = factor(NA),
  #log_sigma_log_R = factor(NA),
  log_sigma_log_F = factor(NA),
  #log_Linf=factor(NA),
  log_t0=factor(NA)
  #log_sigma_index=factor(NA)
)

rnames=c("dev_log_R","dev_log_F","dev_log_N0")

library("TMB")
compile("ALSCL_krill.cpp")
dyn.load("ALSCL_krill")
obj<-MakeADFun(tmb.data,parameters,random=rnames,map=map,DLL="ALSCL_krill",inner.control=list(trace=F, maxit=500))
opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))

opt
obj$gr(opt$par)
cbind(opt$par,lower,upper)
exp(opt$par) 

report<-obj$report()
bound_check<-c((as.vector(opt$par)-as.vector(lower)),(as.vector(upper)-as.vector(opt$par)))
bound_hit<-min(bound_check)==0
final_outer_mgc<-max(abs(obj$gr(opt$par)))
l50=30
l95=38

sdresult=sdreport(obj)
ss=summary(sdresult)

result<-list(obj=obj,opt=opt,report=report,est_std=ss,bound_hit=bound_hit,
             bound_check=bound_check,converge=opt$message,
             final_outer_mgc=final_outer_mgc,l50=l50,l95=l95,year=c(1992:2011),
             len_border=len_border,len_mid=len_mid)
