#The following codes produces the raw RData used in Table 2. "table_1_2_3_merge.R" uses the RData produced by "table_2.R" to produce Table 2 in our manuscript.
#Cluster server is used to run 1 simulation per 1 job with 500 simulations at the same time.

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

#Change n to obtain the results for the Tables 8 and 9 in the Appendix.
n=1000
D=200
b0=0.15
b1=0.75
b2=0.75
b3=-2
s1=0.5
s2=0.5

nknots=20

x1=rnorm(n,0,sqrt(s1))
x2=x1+rnorm(n,.25,sqrt(s2))
Xobs=cbind(x1,x2)
bx=(b0+b1*x1+b2*x2+b3*x1*x2)/3
obsprob=expit(bx)
M=rbinom(n,1,obsprob)
nobs=sum(M==1)
	
E=rnorm(n,mean=0,sd=2)
G=11.875+b1*x1+b2*x2+b3*(x1*x2)^2
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
est_mean_mr=matrix(NA,nrow=D,ncol=4)

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
bs_est_mean_rrz_bart=rep(NA,D)
est_var_bart=matrix(NA,nrow=D,ncol=4)
est_var_bartps=matrix(NA,nrow=D,ncol=4)

obj.1111=function(rho){ 
	ee=0
        for (i in 1:n) ee=ee+M[i]*g.hat[i,]/(1+sum(rho*g.hat[i,]))
        return(sum(ee^2))
}

bootdat=data.frame(y,yobs,M,Xobs)

for(boot in 1:D){
	bsdat=bootdat[sample(n,n,replace=T),]
	y=bsdat$y
	yobs=bsdat$yobs
	M=bsdat$M
	Xobs=cbind(bsdat$x1,bsdat$x2)
	nobs=sum(M==1)

	#For BD and CC
	bs_est_mean_bd[boot]=mean(y)
	bs_est_mean_cc[boot]=mean(yobs[M==1])

	#models for PSPP, AIPWT, PSBPP, and BARTps
	logit_model_s1=glm(M~Xobs[,1]+Xobs[,2]+Xobs[,1]:Xobs[,2],family="binomial")
	logit_ps_s1=logit_model_s1$fitted.values
	logit_ps1_s1=logit_ps_s1[M==1]
	logit_ps0_s1=logit_ps_s1[M==0]
	
	bart_model_s1=bart(Xobs,M,verbose=F)
	tmp_bart_ps_s1=apply(pnorm(bart_model_s1$yhat.train),2,mean)
	bart_ps_s1=tmp_bart_ps_s1
	bart_ps1_s1=bart_ps_s1[M==1]
	bart_ps0_s1=bart_ps_s1[M==0]

	logit_model_s3=glm(M~Xobs[,1]+Xobs[,2],family="binomial")
	logit_ps_s3=logit_model_s3$fitted.values
	logit_ps1_s3=logit_ps_s3[M==1]
	logit_ps0_s3=logit_ps_s3[M==0]		

	#Models for MR
	MRobs=cbind(bsdat$x1,bsdat$x2,bsdat$x1^2,bsdat$x2^2,(bsdat$x1*bsdat$x2)^2)
	propensity.1=glm(M~MRobs[,1]+MRobs[,2]+MRobs[,1]:MRobs[,2],family="binomial") #estimates propensity model 1  
	pai.1=propensity.1$fitted.values
	
    	propensity.2=glm(M~MRobs[,1]+MRobs[,2],family="binomial") #estimates propensity model 2
    	pai.2=propensity.2$fitted.values  
	
	propensity.3=glm(M~MRobs[,3]+MRobs[,4],family="binomial") #estimates propensity model 3
    	pai.3=propensity.3$fitted.values  

    	augmentation.1=glm(yobs~MRobs[,1]+MRobs[,2]+MRobs[,5],family=gaussian)
    	aug.1=(cbind(rep(1,n),MRobs[,1],MRobs[,2],MRobs[,5])%*%as.matrix(augmentation.1$coefficients,,1))[,1] #estimates mean model 1	

	augmentation.2=glm(yobs~MRobs[,1]+MRobs[,2],family=gaussian)
    	aug.2=(cbind(rep(1,n),MRobs[,1],MRobs[,2])%*%as.matrix(augmentation.2$coefficients,,1))[,1] #estimates mean model 1
	
	augmentation.3=glm(yobs~MRobs[,3]+MRobs[,4],family=gaussian)
    	aug.3=(cbind(rep(1,n),MRobs[,3],MRobs[,4])%*%as.matrix(augmentation.3$coefficients,,1))[,1] #estimates mean model 1

	#Ordering of methods and Scenarios based on computation ease.
	#scenario 1
	#PSPP
	knots=quantile(c(min(logit_ps1_s1),max(logit_ps1_s1)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1_s1-knots[i+1])*(logit_ps1_s1>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0_s1-knots[i+1])*(logit_ps0_s1>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psppdat=data.frame(yobs[M==1],logit_ps1_s1,Xobs[M==1,1],Xobs[M==1,2],(Xobs[M==1,1]*Xobs[M==1,2])^2,linear_basis_s1)
	colnames(psppdat)=c("yobs","ps","X1","X2","int",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+int,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s1,Xobs[M==0,1],Xobs[M==0,2],(Xobs[M==0,1]*Xobs[M==0,2])^2,linear_basis0_s1)
	colnames(predpsppdat)=c("ps","X1","X2","int",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,1]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,1]=var(c(pspppred,yobs[M==1]))

	#AIPWT
	rrz=lm(yobs~X1+X2+int,data=psppdat)
	bs_est_mean_rrz_s1[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s1)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))
	
	#MLR
	est_mean_mi[boot,c(1,3)]=mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	est_var_mi[boot,c(1,3)]=var(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))

	#NP
	gamdat=data.frame(yobs[M==1],Xobs[M==1,1],Xobs[M==1,2],(Xobs[M==1,1]*Xobs[M==1,2])^2)
	colnames(gamdat)=c("yobs","X1","X2","int")
	predgamdat=data.frame(Xobs[M==0,1],Xobs[M==0,2],(Xobs[M==0,1]*Xobs[M==0,2])^2)
	colnames(predgamdat)=c("X1","X2","int")
	gammod=gam(yobs~s(X1)+s(X2)+s(int),data=gamdat)
	predy=c(as.numeric(predict(gammod,newdata=predgamdat)+rnorm((n-nobs),0,sqrt(gammod$sig2))),yobs[M==1])
	est_mean_np[boot,1]=mean(predy)
	est_var_np[boot,1]=var(predy)

	#PSBPP
	knots=quantile(c(min(bart_ps1_s1),max(bart_ps1_s1)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_bart_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps1_s1-knots[i+1])*(bart_ps1_s1>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_bart_s1=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps0_s1-knots[i+1])*(bart_ps0_s1>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psbppdat=data.frame(yobs[M==1],bart_ps1_s1,Xobs[M==1,1],Xobs[M==1,2],(Xobs[M==1,1]*Xobs[M==1,2])^2,linear_basis_bart_s1)
	colnames(psbppdat)=c("yobs","ps","X1","X2","int",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2+int,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s1,Xobs[M==0,1],Xobs[M==0,2],(Xobs[M==0,1]*Xobs[M==0,2])^2,linear_basis0_bart_s1)
	colnames(predpsbppdat)=c("ps","X1","X2","int",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,c(1,3)]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,c(1,3)]=var(c(psbpppred,yobs[M==1]))

	#BART
	mismod=bart(x.train=Xobs[M==1,],y.train=yobs[M==1],x.test=Xobs[M==0,],verbose=F)
	bartpred=mismod$yhat.test.mean
	est_mean_bart[boot,]=rep(mean(c(bartpred,yobs[M==1])),4)
	est_var_bart[boot,]=rep(var(c(bartpred,yobs[M==1])),4)

	#AIPWT_BART
	bs_est_mean_rrz_bart[boot]=sum(as.numeric(yobs[M==1]-mismod$yhat.train.mean)/bart_ps1_s1)/n+mean(c(as.numeric(bartpred),as.numeric(mismod$yhat.train.mean)))

	#BARTps
	mismod=bart(x.train=cbind(bart_ps1_s1,Xobs[M==1,]),y.train=yobs[M==1],x.test=cbind(bart_ps0_s1,Xobs[M==0,]),verbose=F)
	bartpspred=mismod$yhat.test.mean
	est_mean_bartps[boot,]=rep(mean(c(bartpspred,yobs[M==1])),4)
	est_var_bartps[boot,]=rep(var(c(bartpspred,yobs[M==1])),4)

	#MR
	g.hat=cbind(pai.1-mean(pai.1),pai.2-mean(pai.2),aug.1-mean(aug.1),aug.2-mean(aug.2))
    	out1=optim(c(0,min(n/nobs,-1/min(g.hat[,2])),0,0),obj.1111,control=list(maxit=1000))
    	out2=optim(c(min(n/nobs,-1/min(g.hat[,1])),0,0,0),obj.1111,control=list(maxit=1000))
    	out3=optim(c(n/(2*nobs),n/(2*nobs),0,0),obj.1111,control=list(maxit=1000))
    	f.1111=min(out1[[2]],out2[[2]],out3[[2]])
    	if(f.1111==out1[[2]]){	
		out=out1
	} else if(f.1111==out2[[2]]){
		out=out2
	} else if(f.1111==out3[[2]]){
		out=out3
	}
    	rho.1111=out[[1]]
    	weight=(1/(nobs*(1+g.hat%*%matrix(rho.1111,,1))))[,1]
    	weight=weight/sum(M*weight)
    	est_mean_mr[boot,1]=sum(M*yobs*weight,na.rm=T)

	#scenario 3
	#PSPP
	knots=quantile(c(min(logit_ps1_s3),max(logit_ps1_s3)), probs=seq(0,1,1/(nknots+1)))
	linear_basis_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1_s3-knots[i+1])*(logit_ps1_s3>knots[i+1]))),nrow=nobs,ncol=nknots)
	linear_basis0_s3=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0_s3-knots[i+1])*(logit_ps0_s3>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
	psppdat=data.frame(yobs[M==1],logit_ps1_s3,Xobs[M==1,1],Xobs[M==1,2],(Xobs[M==1,1]*Xobs[M==1,2])^2,linear_basis_s3)
	colnames(psppdat)=c("yobs","ps","X1","X2","int",paste0("basis",1:nknots))
	psppdat$id=1
	pspp=lme(yobs~ps+X1+X2+int,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s3,Xobs[M==0,1],Xobs[M==0,2],(Xobs[M==0,1]*Xobs[M==0,2])^2,linear_basis0_s3)
	colnames(predpsppdat)=c("ps","X1","X2","int",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,3]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,3]=var(c(pspppred,yobs[M==1]))

	#AIPWT_BART
	bs_est_mean_rrz_s3[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s3)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))

	#MR
	g.hat=cbind(pai.3-mean(pai.3),pai.2-mean(pai.2),aug.1-mean(aug.1),aug.2-mean(aug.2))
	out1=optim(c(0,min(n/nobs,-1/min(g.hat[,2])),0,0),obj.1111,control=list(maxit=1000))
	out2=optim(c(min(n/nobs,-1/min(g.hat[,1])),0,0,0),obj.1111,control=list(maxit=1000))
	out3=optim(c(n/(2*nobs),n/(2*nobs),0,0),obj.1111,control=list(maxit=1000))
	f.1111=min(out1[[2]],out2[[2]],out3[[2]])
	if(f.1111==out1[[2]]){	
		out=out1
	} else if(f.1111==out2[[2]]){
		out=out2
	} else if(f.1111==out3[[2]]){
		out=out3
	}
	rho.1111=out[[1]]
	weight=(1/(nobs*(1+g.hat%*%matrix(rho.1111,,1))))[,1]
	weight=weight/sum(M*weight)
	est_mean_mr[boot,3]=sum(M*yobs*weight,na.rm=T)

	#scenario 2
	#PSPP	
	psppdat=data.frame(yobs[M==1],logit_ps1_s1,Xobs[M==1,1],Xobs[M==1,2],linear_basis_s1)
	colnames(psppdat)=c("yobs","ps","X1","X2",paste0("basis",1:nknots))
	psppdat$id=1	
	pspp=lme(yobs~ps+X1+X2,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s1,Xobs[M==0,1],Xobs[M==0,2],linear_basis0_s1)
	colnames(predpsppdat)=c("ps","X1","X2",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,2]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,2]=var(c(pspppred,yobs[M==1]))

	#AIPWT
	rrz=lm(yobs~X1+X2,data=psppdat)
	bs_est_mean_rrz_s2[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s1)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))
	
	#MLR
	est_mean_mi[boot,c(2,4)]=mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	est_var_mi[boot,c(2,4)]=var(c(as.numeric(predict(rrz,newdata=predpsppdat)),yobs[M==1]))	
	
	#NP
	gammod=gam(yobs~s(X1)+s(X2),data=gamdat)
	predy=c(as.numeric(predict(gammod,newdata=predgamdat)+rnorm((n-nobs),0,sqrt(gammod$sig2))),yobs[M==1])
	est_mean_np[boot,2]=mean(predy)
	est_var_np[boot,2]=var(predy)

	#PSBPP
	psbppdat=data.frame(yobs[M==1],bart_ps1_s1,Xobs[M==1,1],Xobs[M==1,2],linear_basis_bart_s1)
	colnames(psbppdat)=c("yobs","ps","X1","X2",paste0("basis",1:nknots))
	psbppdat$id=1
	psbpp=lme(yobs~ps+X1+X2,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
	predpsbppdat=data.frame(bart_ps0_s1,Xobs[M==0,1],Xobs[M==0,2],linear_basis0_bart_s1)
	colnames(predpsbppdat)=c("ps","X1","X2",paste0("basis",1:nknots))
	predpsbppdat$id=1
	psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))
	est_mean_psbpp[boot,c(2,4)]=mean(c(psbpppred,yobs[M==1]))
	est_var_psbpp[boot,c(2,4)]=var(c(psbpppred,yobs[M==1]))

	#MR
	g.hat=cbind(pai.1-mean(pai.1),pai.2-mean(pai.2),aug.3-mean(aug.3),aug.2-mean(aug.2))
    	out1=optim(c(0,min(n/nobs,-1/min(g.hat[,2])),0,0),obj.1111,control=list(maxit=1000))
    	out2=optim(c(min(n/nobs,-1/min(g.hat[,1])),0,0,0),obj.1111,control=list(maxit=1000))
    	out3=optim(c(n/(2*nobs),n/(2*nobs),0,0),obj.1111,control=list(maxit=1000))
    	f.1111=min(out1[[2]],out2[[2]],out3[[2]])
    	if(f.1111==out1[[2]]){	
		out=out1
	} else if(f.1111==out2[[2]]){
		out=out2
	} else if(f.1111==out3[[2]]){
		out=out3
	}
    	rho.1111=out[[1]]
    	weight=(1/(nobs*(1+g.hat%*%matrix(rho.1111,,1))))[,1]
    	weight=weight/sum(M*weight)
    	est_mean_mr[boot,2]=sum(M*yobs*weight,na.rm=T)

	#scenario 4
	#PSPP
	psppdat=data.frame(yobs[M==1],logit_ps1_s3,Xobs[M==1,1],Xobs[M==1,2],linear_basis_s3)
	colnames(psppdat)=c("yobs","ps","X1","X2",paste0("basis",1:nknots))
	psppdat$id=1	
	pspp=lme(yobs~ps+X1+X2,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
	predpsppdat=data.frame(logit_ps0_s3,Xobs[M==0,1],Xobs[M==0,2],linear_basis0_s3)
	colnames(predpsppdat)=c("ps","X1","X2",paste0("basis",1:nknots))
	predpsppdat$id=1
	pspppred=as.numeric(predict(pspp,newdata=predpsppdat))
	est_mean_pspp[boot,4]=mean(c(pspppred,yobs[M==1]))
	est_var_pspp[boot,4]=var(c(pspppred,yobs[M==1]))

	#AIPWT
	bs_est_mean_rrz_s4[boot]=sum(as.numeric(rrz$residuals)/logit_ps1_s3)/n+mean(c(as.numeric(predict(rrz,newdata=predpsppdat)),as.numeric(rrz$fitted.values)))

	#MR
	g.hat=cbind(pai.3-mean(pai.3),pai.2-mean(pai.2),aug.3-mean(aug.3),aug.2-mean(aug.2))
	out1=optim(c(0,min(n/nobs,-1/min(g.hat[,2])),0,0),obj.1111,control=list(maxit=1000))
	out2=optim(c(min(n/nobs,-1/min(g.hat[,1])),0,0,0),obj.1111,control=list(maxit=1000))
	out3=optim(c(n/(2*nobs),n/(2*nobs),0,0),obj.1111,control=list(maxit=1000))
	f.1111=min(out1[[2]],out2[[2]],out3[[2]])
	if(f.1111==out1[[2]]){	
		out=out1
	} else if(f.1111==out2[[2]]){
		out=out2
	} else if(f.1111==out3[[2]]){
		out=out3
	}
	rho.1111=out[[1]]
	weight=(1/(nobs*(1+g.hat%*%matrix(rho.1111,,1))))[,1]
	weight=weight/sum(M*weight)
	est_mean_mr[boot,4]=sum(M*yobs*weight,na.rm=T)
}

lower_rrz=c(as.numeric(quantile(bs_est_mean_rrz_s1,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s2,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s3,probs=0.025)),as.numeric(quantile(bs_est_mean_rrz_s4,probs=0.025)))
lower_rrz_bart=rep(as.numeric(quantile(bs_est_mean_rrz_bart,probs=0.025)),4)
upper_rrz=c(as.numeric(quantile(bs_est_mean_rrz_s1,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s2,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s3,probs=0.975)),as.numeric(quantile(bs_est_mean_rrz_s4,probs=0.975)))
upper_rrz_bart=rep(as.numeric(quantile(bs_est_mean_rrz_bart,probs=0.975)),4)
lower_mr=c(as.numeric(quantile(est_mean_mr[,1],probs=0.025)),as.numeric(quantile(est_mean_mr[,2],probs=0.025)),as.numeric(quantile(est_mean_mr[,3],probs=0.025)),as.numeric(quantile(est_mean_mr[,4],probs=0.025)))	
upper_mr=c(as.numeric(quantile(est_mean_mr[,1],probs=0.975)),as.numeric(quantile(est_mean_mr[,2],probs=0.975)),as.numeric(quantile(est_mean_mr[,3],probs=0.975)),as.numeric(quantile(est_mean_mr[,4],probs=0.975)))

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

keep(server_num,est_mean_bd,est_mean_cc,mean_mi,mean_pspp,mean_psbpp,est_mean_rrz,est_mean_rrz_bart,mean_bart,mean_bartps,upper_bd,upper_cc,upper_mi,upper_pspp,upper_psbpp,upper_rrz,upper_rrz_bart,upper_bart,upper_bartps,lower_bd,lower_cc,lower_mi,lower_pspp,lower_psbpp,lower_rrz,lower_rrz_bart,lower_bart,lower_bartps,est_mean_mr,lower_mr,upper_mr,mean_np,upper_np,lower_np)

save.image(paste0("table_2_",server_num,".RData"))
