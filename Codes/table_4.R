#The following codes produces the raw RData used in Table 4. "table_4_merge.R" uses the RData produced by "table_4.R" to produce Table 4 in our manuscript.
#Cluster server is used to run 1 simulation per 1 job with 250 simulations at the same time.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

#Takes in the job number assigned from the server
server_num=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

options(java.parameters = "-Xmx14g")

#load the dataset provided. Original 2014 NASS-CDS is available online. Loaded dataset has already been pre-processed. 
load("findat.RData")

library(survey)
library(bartMachine)
library(nlme)
library(mice)
library(polyapost)

set.seed(93467)
N=80000
Bt1=50
Bt2=5
Mt=5
nknots=20

##Step 1: Generate synthetic populations with missing data; 
#Stage 1: Create bootstrap samples from the parent sample;
dsgn=svydesign(ids=~PSU,strata=~CASEID,nest=TRUE,data=findat,weights=~RATWGT) 
dsgn.RW=as.svrepdesign(design=dsgn,type="bootstrap",replicates=Bt1) 
repwt=as.matrix(dsgn.RW$repweights) 

#parallelize Bt1
server_num=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
col_num=ceiling(server_num/5)

st.BB=cbind(findat,repwt[,col_num]) 
colnames(st.BB)=c(colnames(findat),"REPWT")
#delete those units with zero weights for each bootstrap sample; 
st.BB=st.BB[st.BB$REPWT!=0,]
#need to calculate the replicate weights; 
Samwt=st.BB$RATWGT*st.BB[,ncol(st.BB)] 
#normalize again the adjusted weights; 
Samwts=Samwt*N/sum(Samwt) 
np=nrow(st.BB) 
ids=seq(np) 
ns=N-np 

##Stage 2: Create unweighted synthetic populations within each bootstrap sample; 
set.seed(server_num)

l=vector() 
smp=wtpolyap(ids,Samwts,ns) 
#input the adjusted weights in the weighted Polya sampling algorithm; 
for (k in 1:np) l=c(l,length(smp[smp==k])) 

bootdat=matrix(0,nrow=N,ncol=ncol(findat))
bootdat=data.frame(bootdat)
for(loop in 1:(ncol(st.BB)-1)) bootdat[,loop]=rep(st.BB[,loop],l)
colnames(bootdat)=colnames(findat)

Xdat=bootdat[,8:ncol(bootdat)]
yout=bootdat$DVTOTAL
rout=bootdat$RESPOND

temp1=data.frame(yout,Xdat)
temp1_imp=mice(temp1,method="norm", m=Mt) 
ml=complete(temp1_imp, "long") 
	
ccXdat=Xdat[rout==1,]
ccyout=yout[rout==1]
missdat=Xdat[rout==0,]

n=nrow(bootdat)
nobs=length(ccyout)
	
logit_model=glm(rout~.,family="binomial",data=cbind.data.frame(rout,Xdat))
logit_ps=logit_model$fitted.values
logit_ps1=logit_ps[rout==1]
logit_ps0=logit_ps[rout==0]
	
knots=quantile(c(min(logit_ps1),max(logit_ps1)), probs=seq(0,1,1/(nknots+1)))
linear_basis=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1-knots[i+1])*(logit_ps1>knots[i+1]))),nrow=nobs,ncol=nknots)
linear_basis0=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0-knots[i+1])*(logit_ps0>knots[i+1]))),nrow=(n-nobs),ncol=nknots)

psppdat=data.frame(ccyout,ccXdat,logit_ps1,linear_basis)
colnames(psppdat)=c("DVTOTAL",colnames(ccXdat),"ps",paste0("basis",1:nknots))		
psppdat$id=1
pspp=lme(DVTOTAL~ps+ACCTYPE+ANGTHIS+CLIMATE+BODYTYPE+CURBWGT+DOCTRAJ+DRIVDIST+DRINKING+PREILOC+LANES+LGTCOND+MAKE+MANEUVER+MODELYR+OCUPANTS+PREEVENT+PREMOVE+PREISTAB+ALIGNMNT+SURCOND+SURTYPE+RACE+RELINTER+SPECOTH+TRAFCONT+TRAVELSP+TRAFFLOW+DRUGS+ROOF1+ANTILOCK+DAYRUNLT+otbdytyp+DIRDAMW+EXTENT1+EXTENT2+OBJCONT2+PDOF1+AINJSER+AGE+BAGAVRPT+HEIGHT+INJSEV+PARUSE+SEX+WEIGHT,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)
predpsppdat=data.frame(missdat,logit_ps0,linear_basis0)
colnames(predpsppdat)=c(colnames(missdat),"ps",paste0("basis",1:nknots))
predpsppdat$id=1
pspppred=as.numeric(predict(pspp,newdata=predpsppdat))

rrz=lm(DVTOTAL~ACCTYPE+ANGTHIS+CLIMATE+BODYTYPE+CURBWGT+DOCTRAJ+DRIVDIST+DRINKING+PREILOC+LANES+LGTCOND+MAKE+MANEUVER+MODELYR+OCUPANTS+PREEVENT+PREMOVE+PREISTAB+ALIGNMNT+SURCOND+SURTYPE+RACE+RELINTER+SPECOTH+TRAFCONT+TRAVELSP+TRAFFLOW+DRUGS+ROOF1+ANTILOCK+DAYRUNLT+otbdytyp+DIRDAMW+EXTENT1+EXTENT2+OBJCONT2+PDOF1+AINJSER+AGE+BAGAVRPT+HEIGHT+INJSEV+PARUSE+SEX+WEIGHT,data=psppdat)

#note that bartMachine treats 1 as reference, 0 as response. 
bart_model=bartMachine(X=Xdat,y=factor(rout),verbose=F,run_in_sample=F)
bartpsmodel=bart_machine_get_posterior(bart_model,Xdat)
bart_ps=1-bartpsmodel$y_hat	
bart_ps1=bart_ps[rout==1]
bart_ps0=bart_ps[rout==0]

knots=quantile(c(min(bart_ps1),max(bart_ps1)), probs=seq(0,1,1/(nknots+1)))
linear_basis_bart=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps1-knots[i+1])*(bart_ps1>knots[i+1]))),nrow=nobs,ncol=nknots)
linear_basis0_bart=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps0-knots[i+1])*(bart_ps0>knots[i+1]))),nrow=(n-nobs),ncol=nknots)
	
psbppdat=data.frame(ccyout,ccXdat,bart_ps1,linear_basis_bart)
colnames(psbppdat)=c("DVTOTAL",colnames(ccXdat),"ps",paste0("basis",1:nknots))
psbppdat$id=1
psbpp=lme(DVTOTAL~ps+ACCTYPE+ANGTHIS+CLIMATE+BODYTYPE+CURBWGT+DOCTRAJ+DRIVDIST+DRINKING+PREILOC+LANES+LGTCOND+MAKE+MANEUVER+MODELYR+OCUPANTS+PREEVENT+PREMOVE+PREISTAB+ALIGNMNT+SURCOND+SURTYPE+RACE+RELINTER+SPECOTH+TRAFCONT+TRAVELSP+TRAFFLOW+DRUGS+ROOF1+ANTILOCK+DAYRUNLT+otbdytyp+DIRDAMW+EXTENT1+EXTENT2+OBJCONT2+PDOF1+AINJSER+AGE+BAGAVRPT+HEIGHT+INJSEV+PARUSE+SEX+WEIGHT,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)
predpsbppdat=data.frame(missdat,bart_ps0,linear_basis0_bart)
colnames(predpsbppdat)=c(colnames(missdat),"ps",paste0("basis",1:nknots))
predpsbppdat$id=1
psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))

tmpcc=data.frame(bart_ps1,ccXdat)
tmpms=data.frame(bart_ps0,missdat)
colnames(tmpms)=colnames(tmpcc)
bartps=bartMachine(X=tmpcc,y=ccyout,run_in_sample=F,verbose=F)
bartpspred=bart_machine_get_posterior(bartps,tmpms)$y_hat

CC_mean=mean(bootdat$DVTOTAL,na.rm=T)
MI_mean=mean(ml$yout)
PSPP_mean=mean(c(ccyout,pspppred))
AIPWT_mean=mean(c(as.numeric(rrz$residuals)/logit_ps1+as.numeric(rrz$fitted.values),as.numeric(predict(rrz,newdata=predpsppdat))))
PSBPP_mean=mean(c(ccyout,psbpppred))
BARTps_mean=mean(c(ccyout,bartpspred))

CCdat=data.frame(bootdat$DVTOTAL,ifelse(bootdat$INJSEV!=0,1,0),ifelse((bootdat$INJSEV>=3)&(bootdat$INJSEV<=4),1,0))
MIdat=data.frame(ml$yout,ifelse(ml$INJSEV!=0,1,0),ifelse((ml$INJSEV>=3)&(ml$INJSEV<=4),1,0))
PSPPdat=data.frame(c(ccyout,pspppred),c(ifelse(ccXdat$INJSEV!=0,1,0),ifelse(missdat$INJSEV!=0,1,0)),c(ifelse((ccXdat$INJSEV>=3)&(ccXdat$INJSEV<=4),1,0),ifelse((missdat$INJSEV>=3)&(missdat$INJSEV<=4),1,0)))
AIPWTdat=data.frame(c(as.numeric(rrz$residuals)/logit_ps1+as.numeric(rrz$fitted.values),as.numeric(predict(rrz,newdata=predpsppdat))),c(ifelse(ccXdat$INJSEV!=0,1,0),ifelse(missdat$INJSEV!=0,1,0)),c(ifelse((ccXdat$INJSEV>=3)&(ccXdat$INJSEV<=4),1,0),ifelse((missdat$INJSEV>=3)&(missdat$INJSEV<=4),1,0)))
PSBPPdat=data.frame(c(ccyout,psbpppred),c(ifelse(ccXdat$INJSEV!=0,1,0),ifelse(missdat$INJSEV!=0,1,0)),c(ifelse((ccXdat$INJSEV>=3)&(ccXdat$INJSEV<=4),1,0),ifelse((missdat$INJSEV>=3)&(missdat$INJSEV<=4),1,0)))
BARTpsdat=data.frame(c(ccyout,bartpspred),c(ifelse(ccXdat$INJSEV!=0,1,0),ifelse(missdat$INJSEV!=0,1,0)),c(ifelse((ccXdat$INJSEV>=3)&(ccXdat$INJSEV<=4),1,0),ifelse((missdat$INJSEV>=3)&(missdat$INJSEV<=4),1,0)))

colnames(CCdat)=c("DV","IS_Y","IS_sev")
colnames(MIdat)=c("DV","IS_Y","IS_sev")
colnames(PSPPdat)=c("DV","IS_Y","IS_sev")
colnames(AIPWTdat)=c("DV","IS_Y","IS_sev")
colnames(PSBPPdat)=c("DV","IS_Y","IS_sev")
colnames(BARTpsdat)=c("DV","IS_Y","IS_sev")

tmpdvcc=numeric(length(CCdat$DV))
tmpdvcc[is.na(CCdat$DV)]=NA
tmpdvcc[(CCdat$DV>=15)&(CCdat$DV<=35)]=1
tmpdvcc[CCdat$DV>35]=2	

tmpdvmi=numeric(length(MIdat$DV))
tmpdvmi[is.na(MIdat$DV)]=NA
tmpdvmi[(MIdat$DV>=15)&(MIdat$DV<=35)]=1
tmpdvmi[MIdat$DV>35]=2	

tmpdvpspp=numeric(length(PSPPdat$DV))
tmpdvpspp[is.na(PSPPdat$DV)]=NA
tmpdvpspp[(PSPPdat$DV>=15)&(PSPPdat$DV<=35)]=1
tmpdvpspp[PSPPdat$DV>35]=2

tmpdvaipwt=numeric(length(AIPWTdat$DV))
tmpdvaipwt[is.na(AIPWTdat$DV)]=NA
tmpdvaipwt[(AIPWTdat$DV>=15)&(AIPWTdat$DV<=35)]=1
tmpdvaipwt[AIPWTdat$DV>35]=2

tmpdvpsbpp=numeric(length(PSBPPdat$DV))
tmpdvpsbpp[is.na(PSBPPdat$DV)]=NA
tmpdvpsbpp[(PSBPPdat$DV>=15)&(PSBPPdat$DV<=35)]=1
tmpdvpsbpp[PSBPPdat$DV>35]=2

tmpdvbartps=numeric(length(BARTpsdat$DV))
tmpdvbartps[is.na(BARTpsdat$DV)]=NA
tmpdvbartps[(BARTpsdat$DV>=15)&(BARTpsdat$DV<=35)]=1
tmpdvbartps[BARTpsdat$DV>35]=2

CC_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvcc),family="binomial",data=CCdat)$coefficients[2])
MI_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvmi),family="binomial",data=MIdat)$coefficients[2])
PSPP_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvpspp),family="binomial",data=PSPPdat)$coefficients[2])
AIPWT_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvaipwt),family="binomial",data=AIPWTdat)$coefficients[2])
PSBPP_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvpsbpp),family="binomial",data=PSBPPdat)$coefficients[2])
BARTps_OR_IS_any1=exp(glm(IS_Y~as.factor(tmpdvbartps),family="binomial",data=BARTpsdat)$coefficients[2])

CC_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvcc),family="binomial",data=CCdat)$coefficients[3])
MI_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvmi),family="binomial",data=MIdat)$coefficients[3])
PSPP_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvpspp),family="binomial",data=PSPPdat)$coefficients[3])
AIPWT_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvaipwt),family="binomial",data=AIPWTdat)$coefficients[3])
PSBPP_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvpsbpp),family="binomial",data=PSBPPdat)$coefficients[3])
BARTps_OR_IS_any2=exp(glm(IS_Y~as.factor(tmpdvbartps),family="binomial",data=BARTpsdat)$coefficients[3])

CC_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvcc),family="binomial",data=CCdat)$coefficients[2])
MI_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvmi),family="binomial",data=MIdat)$coefficients[2])
PSPP_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvpspp),family="binomial",data=PSPPdat)$coefficients[2])
AIPWT_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvaipwt),family="binomial",data=AIPWTdat)$coefficients[2])
PSBPP_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvpsbpp),family="binomial",data=PSBPPdat)$coefficients[2])
BARTps_OR_IS_sev1=exp(glm(IS_sev~as.factor(tmpdvbartps),family="binomial",data=BARTpsdat)$coefficients[2])

CC_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvcc),family="binomial",data=CCdat)$coefficients[3])
MI_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvmi),family="binomial",data=MIdat)$coefficients[3])
PSPP_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvpspp),family="binomial",data=PSPPdat)$coefficients[3])
AIPWT_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvaipwt),family="binomial",data=AIPWTdat)$coefficients[3])
PSBPP_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvpsbpp),family="binomial",data=PSBPPdat)$coefficients[3])
BARTps_OR_IS_sev2=exp(glm(IS_sev~as.factor(tmpdvbartps),family="binomial",data=BARTpsdat)$coefficients[3])

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

keep(direct,server_num,CC_mean,MI_mean,PSPP_mean,AIPWT_mean,PSBPP_mean,BARTps_mean,CC_OR_IS_any1,MI_OR_IS_any1,PSPP_OR_IS_any1,AIPWT_OR_IS_any1,PSBPP_OR_IS_any1,BARTps_OR_IS_any1,CC_OR_IS_sev1,MI_OR_IS_sev1,PSPP_OR_IS_sev1,AIPWT_OR_IS_sev1,PSBPP_OR_IS_sev1,BARTps_OR_IS_sev1,CC_OR_IS_any2,MI_OR_IS_any2,PSPP_OR_IS_any2,AIPWT_OR_IS_any2,PSBPP_OR_IS_any2,BARTps_OR_IS_any2,CC_OR_IS_sev2,MI_OR_IS_sev2,PSPP_OR_IS_sev2,AIPWT_OR_IS_sev2,PSBPP_OR_IS_sev2,BARTps_OR_IS_sev2)

save.image(paste0("table_4_",server_num,".RData"))
