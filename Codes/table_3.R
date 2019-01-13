#The following codes produces the raw RData used in Table 3. "table_1_2_3_merge.R" uses the RData produced by "table_3.R" to produce Table 3 in our manuscript.
#Cluster server is used to run 1 simulation per 1 job with 500 simulations at the same time.
#Note that the code source for our implementation of MR has requested we not share the codes.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

#Takes in the job number assigned from the server
server_num=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

set.seed(server_num)

library(BayesTree)
library(nlme)
library(mgcv)

expit=function(x){
	return(exp(x)/(1+exp(x)))
}

n=500
D=200
p=4
nknots=20

Z=matrix(rnorm(n*p),ncol=p)
obsprob=expit(-Z[,1]+0.5*Z[,2]-0.25*Z[,3]-0.1*Z[,4])
M=rbinom(n,1,obsprob)
nobs=sum(M==1)	
Xobs=cbind(exp(Z[,1]/2),Z[,2]/(1+exp(Z[,1])),((Z[,1]*Z[,3])/25+0.6)^3,(Z[,2]+Z[,4]+20)^2)
E=rnorm(n)
G=210+27.4*Z[,1]+13.7*Z[,2]+13.7*Z[,3]+13.7*Z[,4]
y=G+E
yobs=y
yobs[M==0]=NA
	
est_mean_bd=rep(mean(y),4)
est_mean_cc=rep(mean(yobs,na.rm=T),4)
est_mean_mi=matrix(NA,nrow=D,ncol=4)
est_mean_np=matrix(NA,nrow=D,ncol=2)
est_mean_pspp=matrix(NA,nrow=D,ncol=4)
est_mean_psbpp=matrix(NA,nrow=D,ncol=4)
est_mean_rrz=rep(NA,4)
est_mean_rrz_bart=rep(NA,4)
est_mean_bart=matrix(NA,nrow=D,ncol=4)
est_mean_bartps=matrix(NA,nrow=D,ncol=4)

bs_est_mean_bd=rep(NA,D)
bs_est_mean_cc=rep(NA,D)
est_var_mi=matrix(NA,nrow=D,ncol=4)
est_var_np=matrix(NA,nrow=D,ncol=2)
est_var_pspp=matrix(NA,nrow=D,ncol=4)
est_var_psbpp=matrix(NA,nrow=D,ncol=4)
bs_est_mean_rrz_s1=rep(NA,D)
bs_est_mean_rrz_s2=rep(NA,D)
bs_est_mean_rrz_s3=rep(NA,D)
bs_est_mean_rrz_s4=rep(NA,D)
bs_est_mean_rrz_bart_s1=rep(NA,D)
bs_est_mean_rrz_bart_s2=rep(NA,D)
bs_est_mean_rrz_bart_s3=rep(NA,D)
bs_est_mean_rrz_bart_s4=rep(NA,D)
est_var_bart=matrix(NA,nrow=D,ncol=4)
est_var_bartps=matrix(NA,nrow=D,ncol=4)

bootdat=data.frame(y,yobs,M,Z,Xobs)

for(boot in 1:D){
	bsdat=bootdat[sample(n,n,replace=T),]
	y=bsdat$y
	yobs=bsdat$yobs
	M=bsdat$M
	Z=cbind(bsdat$X1,bsdat$X2,bsdat$X3,bsdat$X4)
	Xobs=cbind(bsdat$X1.1,bsdat$X2.1,bsdat$X3.1,bsdat$X4.1)
	nobs=sum(M==1)

	#CC and BD
	bs_est_mean_bd[boot]=mean(y)
	bs_est_mean_cc[boot]=mean(yobs[M==1])

	#Models for PSPP, AIPWT, PSBPP, and BARTps
	logit_model_s1=glm(M~Z[,1]+Z[,2]+Z[,3]+Z[,4],family="binomial")
	logit_ps_s1=logit_model_s1$fitted.values
	logit_ps1_s1=logit_ps_s1[M==1]
	logit_ps0_s1=logit_ps_s1[M==0]

	bart_model_s1=bart(Z,M,verbose=F)
	tmp_bart_ps_s1=apply(pnorm(bart_model_s1$yhat.train),2,mean)
	bart_ps_s1=tmp_bart_ps_s1
	bart_ps1_s1=bart_ps_s1[M==1]
	bart_ps0_s1=bart_ps_s1[M==0]

	logit_model_s3=glm(M~Xobs[,1]+Xobs[,2]+Xobs[,3]+Xobs[,4],family="binomial")
	logit_ps_s3=logit_model_s3$fitted.values
	logit_ps1_s3=logit_ps_s3[M==1]
	logit_ps0_s3=logit_ps_s3[M==0]	

	bart_model_s3=bart(Xobs,M,verbose=F)
	tmp_bart_ps_s3=apply(pnorm(bart_model_s3$yhat.train),2,mean)
	bart_ps_s3=tmp_bart_ps_s3
	bart_ps1_s3=bart_ps_s3[M==1]
	bart_ps0_s3=bart_ps_s3[M==0]	

	#scenario 1
	#PSPP
	knots=quantile(c(min(logit_ps1_s1),max(logit_ps1_s1)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1_s1-knots[i+1])*(logit_ps1_s1>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0_s1-knots[i+1])*(logit_ps0_s1>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psppdat=data.frame(yobs[M==1],logit_ps1_s1,Z[M==1,],linear_basis_s1)
	colnames(psppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s1,Z[M==0,],linear_basis0_s1)
	colnames(predpsppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,1]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,1]=var(c(pspppred,yobs[M==1]))

	#PSBPP
	knots=quantile(c(min(bart_ps1_s1),max(bart_ps1_s1)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_bart_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps1_s1-knots[i+1])*(bart_ps1_s1>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_bart_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps0_s1-knots[i+1])*(bart_ps0_s1>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psbppdat=data.frame(yobs[M==1],bart_ps1_s1,Z[M==1,],linear_basis_bart_s1)
	colnames(psbppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s1,Z[M==0,],linear_basis0_bart_s1)
	colnames(predpsbppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,1]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,1]=var(c(psbpppred,yobs[M==1]))
	
	#AIPWT
	rrz=lm(yobs~X1+X2+X3+X4,data=psppdat)
	bs_est_mean_rrz_s1[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s1)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))
	
	#MLR
	est_mean_mi[boot,c(1,3)]=mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	est_var_mi[boot,c(1,3)]=var(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))

	#NP
	gamdat=data.frame(yobs[M==1],Z[M==1,1],Z[M==1,2],Z[M==1,3],Z[M==1,4])
	colnames(gamdat)=c("yobs","X1","X2","X3","X4")
	predgamdat=data.frame(Z[M==0,1],Z[M==0,2],Z[M==0,3],Z[M==0,4])
	colnames(predgamdat)=c("X1","X2","X3","X4")
	gammod=gam(yobs~s(X1)+s(X2)+s(X3)+s(X4),data=gamdat)
	predy=c(as.numeric(predict(gammod,newdata=predgamdat)+rnorm((n-nobs),0,sqrt(gammod$sig2))),yobs[M==1])
	est_mean_np[boot,1]=mean(predy)
	est_var_np[boot,1]=var(predy)

	#BART
	mismod=bart(x.train=Z[M==1,],y.train=yobs[M==1],x.test=Z[M==0,],verbose=F)
	bartpred=mismod$yhat.test.mean
	est_mean_bart[boot,c(1,3)]=mean(c(bartpred,yobs[M==1]))
	est_var_bart[boot,c(1,3)]=var(c(bartpred,yobs[M==1]))

	#AIPWT_BART
	bs_est_mean_rrz_bart_s1[boot]=sum(as.numeric(yobs[M==1]-mismod$yhat.train.mean)/bart_ps1_s1)/n+mean(c(as.numeric(bartpred),as.numeric(mismod$yhat.train.mean)))
	bs_est_mean_rrz_bart_s3[boot]=sum(as.numeric(yobs[M==1]-mismod$yhat.train.mean)/bart_ps1_s3)/n+mean(c(as.numeric(bartpred),as.numeric(mismod$yhat.train.mean)))

	#BARTps
	mismod=bart(x.train=cbind(bart_ps1_s1,Z[M==1,]),y.train=yobs[M==1],x.test=cbind(bart_ps0_s1,Z[M==0,]),verbose=F)
	bartpspred=mismod$yhat.test.mean
	est_mean_bartps[boot,1]=mean(c(bartpspred+rnorm((n-nobs),0,mean(mismod$sigma)),yobs[M==1]))
	est_var_bartps[boot,1]=var(c(bartpspred+rnorm((n-nobs),0,mean(mismod$sigma)),yobs[M==1]))

	#scenario 3
	#PSPP
	knots=quantile(c(min(logit_ps1_s3),max(logit_ps1_s3)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1_s3-knots[i+1])*(logit_ps1_s3>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0_s3-knots[i+1])*(logit_ps0_s3>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psppdat=data.frame(yobs[M==1],logit_ps1_s3,Z[M==1,],linear_basis_s3)
	colnames(psppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s3,Z[M==0,],linear_basis0_s3)
	colnames(predpsppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,3]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,3]=var(c(pspppred,yobs[M==1]))

	#AIPWT
	bs_est_mean_rrz_s3[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s3)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))

	#PSBPP
	knots=quantile(c(min(bart_ps1_s3),max(bart_ps1_s3)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_bart_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps1_s3-knots[i+1])*(bart_ps1_s3>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_bart_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps0_s3-knots[i+1])*(bart_ps0_s3>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psbppdat=data.frame(yobs[M==1],bart_ps1_s3,Z[M==1,],linear_basis_bart_s3)
	colnames(psbppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s3,Z[M==0,],linear_basis0_bart_s3)
	colnames(predpsbppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,3]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,3]=var(c(psbpppred,yobs[M==1]))

	#BARTps
	mismod=bart(x.train=cbind(bart_ps1_s3,Z[M==1,]),y.train=yobs[M==1],x.test=cbind(bart_ps0_s3,Z[M==0,]),verbose=F)
	bartpspred=mismod$yhat.test.mean
	est_mean_bartps[boot,3]=mean(c(bartpspred,yobs[M==1]))
	est_var_bartps[boot,3]=var(c(bartpspred,yobs[M==1]))

	#scenario 2	
	#PSPP
	psppdat=data.frame(yobs[M==1],logit_ps1_s1,Xobs[M==1,],linear_basis_s1)
	colnames(psppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s1,Xobs[M==0,],linear_basis0_s1)
	colnames(predpsppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,2]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,2]=var(c(pspppred,yobs[M==1]))

	#AIPWT
	rrz=lm(yobs~+X1+X2+X3+X4,data=psppdat)
	bs_est_mean_rrz_s2[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s1)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))
	
	#MLR
	est_mean_mi[boot,c(2,4)]=mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	est_var_mi[boot,c(2,4)]=var(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	
	#NP
	gamdat=data.frame(yobs[M==1],Xobs[M==1,1],Xobs[M==1,2],Xobs[M==1,3],Xobs[M==1,4])
	colnames(gamdat)=c("yobs","X1","X2","X3","X4")
	predgamdat=data.frame(Xobs[M==0,1],Xobs[M==0,2],Xobs[M==0,3],Xobs[M==0,4])
	colnames(predgamdat)=c("X1","X2","X3","X4")
	gammod=gam(yobs~s(X1)+s(X2)+s(X3)+s(X4),data=gamdat)
	predy=c(as.numeric(predict(gammod,newdata=predgamdat)+rnorm((n-nobs),0,sqrt(gammod$sig2))),yobs[M==1])
	est_mean_np[boot,2]=mean(predy)
	est_var_np[boot,2]=var(predy)

	#PSBPP
	psbppdat=data.frame(yobs[M==1],bart_ps1_s1,Xobs[M==1,],linear_basis_bart_s1)
	colnames(psbppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s1,Xobs[M==0,],linear_basis0_bart_s1)
	colnames(predpsbppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,c(2,4)]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,c(2,4)]=var(c(psbpppred,yobs[M==1]))

	#BART
	mismod=bart(x.train=Xobs[M==1,],y.train=yobs[M==1],x.test=Xobs[M==0,],verbose=F)
	bartpred=mismod$yhat.test.mean
	est_mean_bart[boot,c(2,4)]=mean(c(bartpred,yobs[M==1]))
	est_var_bart[boot,c(2,4)]=var(c(bartpred,yobs[M==1]))
	
	#AIPWT_BART
	bs_est_mean_rrz_bart_s2[boot]=sum(as.numeric(yobs[M==1]-mismod$yhat.train.mean)/bart_ps1_s1)/n+mean(c(as.numeric(bartpred),as.numeric(mismod$yhat.train.mean)))
	bs_est_mean_rrz_bart_s4[boot]=sum(as.numeric(yobs[M==1]-mismod$yhat.train.mean)/bart_ps1_s3)/n+mean(c(as.numeric(bartpred),as.numeric(mismod$yhat.train.mean)))
	
	#BARTps
	mismod=bart(x.train=cbind(bart_ps1_s1,Xobs[M==1,]),y.train=yobs[M==1],x.test=cbind(bart_ps0_s1,Xobs[M==0,]),verbose=F)
	bartpspred=mismod$yhat.test.mean
	est_mean_bartps[boot,2]=mean(c(bartpspred,yobs[M==1]))
	est_var_bartps[boot,2]=var(c(bartpspred,yobs[M==1]))
	
	#scenario 4
	#PSPP
	psppdat=data.frame(yobs[M==1],logit_ps1_s3,Xobs[M==1,],linear_basis_s3)
	colnames(psppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s3,Xobs[M==0,],linear_basis0_s3)
	colnames(predpsppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,4]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,4]=mean(c(pspppred,yobs[M==1]))

	#AIPWT
	bs_est_mean_rrz_s4[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s3)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))

	#PSBPP
	psbppdat=data.frame(yobs[M==1],bart_ps1_s3,Xobs[M==1,],linear_basis_bart_s3)
	colnames(psbppdat)=c("yobs","ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2+X3+X4,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s3,Xobs[M==0,],linear_basis0_bart_s3)
	colnames(predpsbppdat)=c("ps","X1","X2","X3","X4",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,4]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,4]=var(c(psbpppred,yobs[M==1]))

	#BARTps
	mismod=bart(x.train=cbind(bart_ps1_s3,Xobs[M==1,]),y.train=yobs[M==1],x.test=cbind(bart_ps0_s3,Xobs[M==0,]),verbose=F)
	bartpspred=mismod$yhat.test.mean
	est_mean_bartps[boot,4]=mean(c(bartpspred,yobs[M==1]))
	est_var_bartps[boot,4]=var(c(bartpspred,yobs[M==1]))	
	
}

lower_rrz=c(as.numeric(quantile(bs_est_mean_rrz_s1,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s2,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s3,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s4,probs=0.025)))
lower_rrz_bart=c(as.numeric(quantile(bs_est_mean_rrz_bart_s1,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_bart_s2,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_bart_s3,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_bart_s4,probs=0.025)))
upper_rrz=c(as.numeric(quantile(bs_est_mean_rrz_s1,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s2,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s3,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s4,probs=0.975)))
upper_rrz_bart=c(as.numeric(quantile(bs_est_mean_rrz_bart_s1,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_bart_s2,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_bart_s3,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_bart_s4,probs=0.975)))

mean_np=colMeans(est_mean_np)
upper_np=rep(NA,2)
lower_np=rep(NA,2)

for(i in 1:2){
	td=mean(est_var_np[,i]/n)+(B+1)/B*var(est_mean_np[,i])
	qnum=qt(0.975,(B-1)*(1+1/(B+1)*(mean(est_var_np[,i]/n))/(var(est_mean_np[,i])))^2)
	upper_np[i]=mean_np[i]+qnum*sqrt(td)
	lower_np[i]=mean_np[i]-qnum*sqrt(td)
}

mean_mi=numeric(4)
mean_pspp=numeric(4)
mean_psbpp=numeric(4)
mean_bart=numeric(4)
mean_bartps=numeric(4)

lower_mi=numeric(4)
lower_pspp=numeric(4)
lower_psbpp=numeric(4)
lower_bart=numeric(4)
lower_bartps=numeric(4)

upper_mi=numeric(4)
upper_pspp=numeric(4)
upper_psbpp=numeric(4)
upper_bart=numeric(4)
upper_bartps=numeric(4)

for(i in 1:4){
	mean_mi[i]=mean(est_mean_mi[,i])
	td=mean(est_var_mi[,i]/n)+(D+1)/D*var(est_mean_mi[,i])
	qnum=qt(0.975,(D-1)*(1+1/(D+1)*(mean(est_var_mi[,i]/n))/(var(est_mean_mi[,i])))^2)
	upper_mi[i]=mean_mi[i]+qnum*sqrt(td)
	lower_mi[i]=mean_mi[i]-qnum*sqrt(td)

	mean_pspp[i]=mean(est_mean_pspp[,i])
	td=mean(est_var_pspp[,i]/n)+(D+1)/D*var(est_mean_pspp[,i])
	qnum=qt(0.975,(D-1)*(1+1/(D+1)*(mean(est_var_pspp[,i]/n))/(var(est_mean_pspp[,i])))^2)
	upper_pspp[i]=mean_pspp[i]+qnum*sqrt(td)
	lower_pspp[i]=mean_pspp[i]-qnum*sqrt(td)

	mean_psbpp[i]=mean(est_mean_psbpp[,i])
	td=mean(est_var_psbpp[,i]/n)+(D+1)/D*var(est_mean_psbpp[,i])
	qnum=qt(0.975,(D-1)*(1+1/(D+1)*(mean(est_var_psbpp[,i]/n))/(var(est_mean_psbpp[,i])))^2)
	upper_psbpp[i]=mean_psbpp[i]+qnum*sqrt(td)
	lower_psbpp[i]=mean_psbpp[i]-qnum*sqrt(td)

	mean_bart[i]=mean(est_mean_bart[,i])
	td=mean(est_var_bart[,i]/n)+(D+1)/D*var(est_mean_bart[,i])
	qnum=qt(0.975,(D-1)*(1+1/(D+1)*(mean(est_var_bart[,i]/n))/(var(est_mean_bart[,i])))^2)
	upper_bart[i]=mean_bart[i]+qnum*sqrt(td)
	lower_bart[i]=mean_bart[i]-qnum*sqrt(td)

	mean_bartps[i]=mean(est_mean_bartps[,i])
	td=mean(est_var_bartps[,i]/n)+(D+1)/D*var(est_mean_bartps[,i])
	qnum=qt(0.975,(D-1)*(1+1/(D+1)*(mean(est_var_bartps[,i]/n))/(var(est_mean_bartps[,i])))^2)
	upper_bartps[i]=mean_bartps[i]+qnum*sqrt(td)
	lower_bartps[i]=mean_bartps[i]-qnum*sqrt(td)
}

keep <- function(..., x = c())
{
  
  if (length(x) > 0)
  {
    if (!is.character(x)) 
      stop("x must contain character vector")
    
    L <- ls(name = parent.frame())
    rm(list = L[!L %in% x], pos = parent.frame())
    
    return(invisible(ls(name = parent.frame())))
  }
  
  dots <- match.call(expand.dots = FALSE)$...
  
  if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                  is.character(x)))) 
    stop("... must contain names or character strings")
  
  names <- sapply(dots, as.character)
  L <- ls(name = parent.frame())
  rm(list = L[!L %in% names], pos = parent.frame())
  
  return(invisible(ls(name = parent.frame())))
  
}

keep(server_num,est_mean_bd,est_mean_cc,mean_mi,mean_pspp,mean_psbpp,est_mean_rrz,est_mean_rrz_bart,mean_bart,mean_bartps,upper_bd,upper_cc,upper_mi,upper_pspp,upper_psbpp,upper_rrz,upper_rrz_bart,upper_bart,upper_bartps,lower_bd,lower_cc,lower_mi,lower_pspp,lower_psbpp,lower_rrz,lower_rrz_bart,lower_bart,lower_bartps,mean_np,upper_np,lower_np)

save.image(paste0("table_3_",server_num,".RData"))
