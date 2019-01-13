#The following codes produces the raw RData used in Table 5. "table_5_merge.R" uses the RData produced by "table_5.R" to produce Table 5 in our manuscript.
#Cluster server is used to run 1 simulation per 1 job with 200 simulations at the same time.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

#Takes in the job number assigned from the server
server_num=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#load the dataset provided. Original 2015 FARS is available online. Loaded dataset has already been pre-processed. 
load("findat_FARS.RData")

options(java.parameters = "-Xmx8g")

library(bartMachine)
library(MASS)
library(nlme)

n=nrow(findat)
tmpdat=findat
findat=tmpdat[sample(n,n,replace=T),]
ccbac=findat$BAC[findat$RESPOND==1]
mibac=rowMeans(cbind(findat$P1,findat$P2,findat$P3,findat$P4,findat$P5,findat$P6,findat$P7,findat$P8,findat$P9,findat$P10))

ccdat=findat[findat$RESPOND==1,]
BACgeq0=as.numeric(ccdat$BAC>0)

Xobs=findat[,16:ncol(findat)]
modelgeq0=bartMachine(X=Xobs[findat$RESPOND==1,],y=factor(BACgeq0),verbose=F,run_in_sample=F)
bartpsmodel=bart_machine_get_posterior(modelgeq0,Xobs[findat$RESPOND==0,])
probgeq0=1-bartpsmodel$y_hat	

geq0=numeric(length(probgeq0))
for(i in 1:length(probgeq0)) geq0[i]=rbinom(1,1,probgeq0[i])
missdat=findat[findat$RESPOND==0,]
ccdatgeq0=ccdat[ccdat$BAC>0,]
missdatgeq0=missdat[geq0>0,]
tmpdat=rbind(ccdatgeq0,missdatgeq0)
bcmodel=boxcox(BAC~HOUR+DAY_WEEK+TYP_INT+WRK_ZONE+REL_ROAD+WEATHER+FATALS+VE_FORMS+ATST_TYP+MONTH+FUNC_SYS+MAN_COLL+MAKE+MOD_YEAR+FIRE_EXP+AGE+SEX+INJ_SEV+REST_USE+AIR_BAG+EXTRICAT+DRINKING+ALC_DET+DRUGS+DRUG_DET+NUMOCCS+HIT_RUN+OWNER+TRAV_SP+UNDERIDE+ROLINLOC+TOWED+M_HARM+VEH_SC1+L_STATUS+L_RESTRI+DR_HGT+DR_WGT+PREV_ACC+SPEEDREL+VTRAFWAY+VNUM_LAN+VSPD_LIM+VALIGN+VPROFILE+VPAVETYP+VSURCOND+VTRAFCON+P_CRASH1+P_CRASH2+P_CRASH3+PCRASH4+PCRASH5+ACC_TYPE+DR_DRINK+DRUGRES+CF+DR_SF,data=tmpdat,plotit=F)
bclambda=bcmodel$x[bcmodel$y==max(bcmodel$y)]+1
if(bclambda!=0){
	tmpdat$BCBAC=(tmpdat$BAC^bclambda-1)/bclambda
	ccdatgeq0$BCBAC=(ccdatgeq0$BAC^bclambda-1)/bclambda
} else {
	tmpdat$BCBAC=log(tmpdat$BAC)
	ccdatgeq0$BCBAC=log(ccdatgeq0$BAC)
}

#PSPP
nknots=20
n=nrow(tmpdat)
nobs=sum(tmpdat$RESPOND)
logit_model=glm(RESPOND~HOUR+DAY_WEEK+TYP_INT+WRK_ZONE+REL_ROAD+WEATHER+FATALS+VE_FORMS+ATST_TYP+MONTH+FUNC_SYS+MAN_COLL+MAKE+MOD_YEAR+FIRE_EXP+AGE+SEX+INJ_SEV+REST_USE+AIR_BAG+EXTRICAT+DRINKING+ALC_DET+DRUGS+DRUG_DET+NUMOCCS+HIT_RUN+OWNER+TRAV_SP+UNDERIDE+ROLINLOC+TOWED+M_HARM+VEH_SC1+L_STATUS+L_RESTRI+DR_HGT+DR_WGT+PREV_ACC+SPEEDREL+VTRAFWAY+VNUM_LAN+VSPD_LIM+VALIGN+VPROFILE+VPAVETYP+VSURCOND+VTRAFCON+P_CRASH1+P_CRASH2+P_CRASH3+PCRASH4+PCRASH5+ACC_TYPE+DR_DRINK+DRUGRES+CF+DR_SF,family="binomial",data=tmpdat)
logit_ps=logit_model$fitted.values
logit_ps1=logit_ps[tmpdat$RESPOND==1]
logit_ps0=logit_ps[tmpdat$RESPOND==0]

knots=quantile(c(min(logit_ps1),max(logit_ps1)), probs=seq(0,1,1/(nknots+1)))
linear_basis=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps1-knots[i+1])*(logit_ps1>knots[i+1]))),nrow=nobs,ncol=nknots)
linear_basis0=matrix(do.call(cbind, lapply(1:nknots, function(i) (logit_ps0-knots[i+1])*(logit_ps0>knots[i+1]))),nrow=(n-nobs),ncol=nknots)		
psppdat=data.frame(ccdatgeq0,logit_ps1,linear_basis)
colnames(psppdat)=c(colnames(ccdatgeq0),"ps",paste0("basis",1:nknots))		
psppdat$id=1
pspp=lme(BCBAC~ps+HOUR+DAY_WEEK+TYP_INT+WRK_ZONE+REL_ROAD+WEATHER+FATALS+VE_FORMS+MONTH+ATST_TYP+FUNC_SYS+MAN_COLL+MAKE+MOD_YEAR+FIRE_EXP+AGE+SEX+INJ_SEV+REST_USE+AIR_BAG+EXTRICAT+DRINKING+ALC_DET+DRUGS+DRUG_DET+NUMOCCS+HIT_RUN+OWNER+TRAV_SP+UNDERIDE+ROLINLOC+TOWED+M_HARM+VEH_SC1+L_STATUS+L_RESTRI+DR_HGT+DR_WGT+PREV_ACC+SPEEDREL+VTRAFWAY+VNUM_LAN+VSPD_LIM+VALIGN+VPROFILE+VPAVETYP+VSURCOND+VTRAFCON+P_CRASH1+P_CRASH2+P_CRASH3+PCRASH4+PCRASH5+ACC_TYPE+DR_DRINK+DRUGRES+CF+DR_SF,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psppdat)		
predpsppdat=data.frame(missdatgeq0,logit_ps0,linear_basis0)
colnames(predpsppdat)=c(colnames(missdatgeq0),"ps",paste0("basis",1:nknots))
predpsppdat$id=1
pspppred=as.numeric(predict(pspp,newdata=predpsppdat))+rnorm(nrow(predpsppdat),mean=0,sd=pspp$sigma)
pspppredgeq0=pspppred[pspppred>0]
if(bclambda!=0){
	psppbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(pspppred<=0)),ccdatgeq0$BAC,(pspppredgeq0*bclambda+1)^(1/bclambda))
} else {
	psppbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(pspppred<=0)),ccdatgeq0$BAC,exp(pspppredgeq0))
}	

#AIPWT
rrz=lm(BCBAC~HOUR+DAY_WEEK+TYP_INT+WRK_ZONE+REL_ROAD+WEATHER+FATALS+VE_FORMS+MONTH+FUNC_SYS+MAN_COLL+ATST_TYP+MAKE+MOD_YEAR+FIRE_EXP+AGE+SEX+INJ_SEV+REST_USE+AIR_BAG+EXTRICAT+DRINKING+ALC_DET+DRUGS+DRUG_DET+NUMOCCS+HIT_RUN+OWNER+TRAV_SP+UNDERIDE+ROLINLOC+TOWED+M_HARM+VEH_SC1+L_STATUS+L_RESTRI+DR_HGT+DR_WGT+PREV_ACC+SPEEDREL+VTRAFWAY+VNUM_LAN+VSPD_LIM+VALIGN+VPROFILE+VPAVETYP+VSURCOND+VTRAFCON+P_CRASH1+P_CRASH2+P_CRASH3+PCRASH4+PCRASH5+ACC_TYPE+DR_DRINK+DRUGRES+CF+DR_SF,data=psppdat)
rrzpred=c(as.numeric(rrz$residuals)/logit_ps1+as.numeric(rrz$fitted.values),as.numeric(predict(rrz,newdata=predpsppdat)))+rnorm(nrow(predpsppdat)+nrow(psppdat),mean=0,sd=summary(rrz)$sigma)
rrzpredgeq0=rrzpred[rrzpred>0]
if(bclambda!=0){
	aipwtbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(rrzpred<=0)),(rrzpredgeq0*bclambda+1)^(1/bclambda))
} else {
	aipwtbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(rrzpred<=0)),exp(rrzpredgeq0))
}	

#PSBPP
Xobs=tmpdat[,16:(ncol(tmpdat)-1)]
bart_model=bartMachine(X=Xobs,y=factor(tmpdat$RESPOND),verbose=F,run_in_sample=F)
bartpos=bart_machine_get_posterior(bart_model,Xobs)
bart_ps=1-bartpos$y_hat
bart_ps1=bart_ps[tmpdat$RESPOND==1]
bart_ps0=bart_ps[tmpdat$RESPOND==0]

knots=quantile(c(min(bart_ps1),max(bart_ps1)), probs=seq(0,1,1/(nknots+1)))
linear_basis_bart=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps1-knots[i+1])*(bart_ps1>knots[i+1]))),nrow=nobs,ncol=nknots)
linear_basis0_bart=matrix(do.call(cbind, lapply(1:nknots, function(i) (bart_ps0-knots[i+1])*(bart_ps0>knots[i+1]))),nrow=(n-nobs),ncol=nknots)	

psbppdat=data.frame(ccdatgeq0,bart_ps1,linear_basis_bart)
colnames(psbppdat)=c(colnames(ccdatgeq0),"ps",paste0("basis",1:nknots))
psbppdat$id=1
psbpp=lme(BCBAC~ps+HOUR+DAY_WEEK+TYP_INT+WRK_ZONE+REL_ROAD+WEATHER+FATALS+VE_FORMS+MONTH+ATST_TYP+FUNC_SYS+MAN_COLL+MAKE+MOD_YEAR+FIRE_EXP+AGE+SEX+INJ_SEV+REST_USE+AIR_BAG+EXTRICAT+DRINKING+ALC_DET+DRUGS+DRUG_DET+NUMOCCS+HIT_RUN+OWNER+TRAV_SP+UNDERIDE+ROLINLOC+TOWED+M_HARM+VEH_SC1+L_STATUS+L_RESTRI+DR_HGT+DR_WGT+PREV_ACC+SPEEDREL+VTRAFWAY+VNUM_LAN+VSPD_LIM+VALIGN+VPROFILE+VPAVETYP+VSURCOND+VTRAFCON+P_CRASH1+P_CRASH2+P_CRASH3+PCRASH4+PCRASH5+ACC_TYPE+DR_DRINK+DRUGRES+CF+DR_SF,random=list(id=pdIdent(~0+basis1+basis2+basis3+basis4+basis5+basis6+basis7+basis8+basis9+basis10+basis11+basis12+basis13+basis14+basis15+basis16+basis17+basis18+basis19+basis20)),data=psbppdat)			
predpsbppdat=data.frame(missdatgeq0,bart_ps0,linear_basis0_bart)
colnames(predpsbppdat)=c(colnames(missdatgeq0),"ps",paste0("basis",1:nknots))
predpsbppdat$id=1
psbpppred=as.numeric(predict(psbpp,newdata=predpsbppdat))+rnorm(nrow(predpsbppdat),mean=0,sd=psbpp$sigma)
psbpppredgeq0=psbpppred[psbpppred>0]
if(bclambda!=0){
	psbppbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(psbpppred<=0)),ccdatgeq0$BAC,(psbpppredgeq0*bclambda+1)^(1/bclambda))
} else {
	psbppbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(psbpppred<=0)),ccdatgeq0$BAC,exp(psbpppredgeq0))
}	
				
#BARTps
tmpcc=data.frame(bart_ps1,Xobs[tmpdat$RESPOND==1,])
tmpms=data.frame(bart_ps0,Xobs[tmpdat$RESPOND==0,])
colnames(tmpms)=colnames(tmpcc)
bartps=bartMachine(X=tmpcc,y=as.numeric(ccdatgeq0$BCBAC),run_in_sample=F,verbose=F)
bartpspred=bart_machine_get_posterior(bartps,tmpms)$y_hat
bartpspredgeq0=bartpspred[bartpspred>0]
if(bclambda!=0){
	bartpsbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(bartpspred<=0)),ccdatgeq0$BAC,(bartpspredgeq0*bclambda+1)^(1/bclambda))
} else {
	bartpsbac=c(numeric(sum(ccdat$BAC==0)+sum(geq0==0)+sum(bartpspred<=0)),ccdatgeq0$BAC,exp(bartpspredgeq0))
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

keep(direct,server_num,ccbac,mibac,psppbac,aipwtbac,psbppbac,bartpsbac)

save.image(paste0("table_5_",server_num,".RData"))
