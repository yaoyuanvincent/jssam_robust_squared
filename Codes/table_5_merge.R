#Takes in results from table_5.R to produce Table 5.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

D=200

tmppsppmean=numeric(D)
tmppsppvar=numeric(D)
tmpaipwtmean=numeric(D)
tmpaipwtvar=numeric(D)
tmppsbppmean=numeric(D)
tmppsbppvar=numeric(D)
tmpbartpsmean=numeric(D)
tmpbartpsvar=numeric(D)

tmppsppmean1=numeric(D)
tmppsppvar1=numeric(D)
tmpaipwtmean1=numeric(D)
tmpaipwtvar1=numeric(D)
tmppsbppmean1=numeric(D)
tmppsbppvar1=numeric(D)
tmpbartpsmean1=numeric(D)
tmpbartpsvar1=numeric(D)

tmppsppmean10=numeric(D)
tmppsppvar10=numeric(D)
tmpaipwtmean10=numeric(D)
tmpaipwtvar10=numeric(D)
tmppsbppmean10=numeric(D)
tmppsbppvar10=numeric(D)
tmpbartpsmean10=numeric(D)
tmpbartpsvar10=numeric(D)

for(iter in 1:D){
	load(paste0("table_5_",iter,".RData"))
	
	#mean
	tmppsppmean[iter]=mean(psppbac)
	tmppsppvar[iter]=var(psppbac)
	tmpaipwtmean[iter]=mean(aipwtbac)
	tmpaipwtvar[iter]=var(aipwtbac)
	tmppsbppmean[iter]=mean(psbppbac)
	tmppsbppvar[iter]=var(psbppbac)
	tmpbartpsmean[iter]=mean(bartpsbac)
	tmpbartpsvar[iter]=var(bartpsbac)

	#geq1
	tmppsppmean1[iter]=mean(psppbac>1)
	tmppsppvar1[iter]=tmppsppmean1[iter]*(1-tmppsppmean1[iter])
	tmpaipwtmean1[iter]=mean(aipwtbac>1)
	tmpaipwtvar1[iter]=tmpaipwtmean1[iter]*(1-tmpaipwtmean1[iter])
	tmppsbppmean1[iter]=mean(psbppbac>1)
	tmppsbppvar1[iter]=tmppsbppmean1[iter]*(1-tmppsbppmean1[iter])
	tmpbartpsmean1[iter]=mean(bartpsbac>1)
	tmpbartpsvar1[iter]=tmpbartpsmean1[iter]*(1-tmpbartpsmean1[iter])

	#geq10
	tmppsppmean10[iter]=mean(psppbac>10)
	tmppsppvar10[iter]=tmppsppmean10[iter]*(1-tmppsppmean10[iter])
	tmpaipwtmean10[iter]=mean(aipwtbac>10)
	tmpaipwtvar10[iter]=tmpaipwtmean10[iter]*(1-tmpaipwtmean10[iter])
	tmppsbppmean10[iter]=mean(psbppbac>10)
	tmppsbppvar10[iter]=tmppsbppmean10[iter]*(1-tmppsbppmean10[iter])
	tmpbartpsmean10[iter]=mean(bartpsbac>10)
	tmpbartpsvar10[iter]=tmpbartpsmean10[iter]*(1-tmpbartpsmean10[iter])

}

load("findat_FARS.RData")

ccbac=findat$BAC[findat$RESPOND==1]
num=nrow(findat)

#mean
ccmean=mean(ccbac)
ccvar=var(ccbac)
ccub1=ccmean+qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))
cclb=ccmean-qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))

mimean=mean(c(mean(findat$P1),mean(findat$P2),mean(findat$P3),mean(findat$P4),mean(findat$P5),mean(findat$P6),mean(findat$P7),mean(findat$P8),mean(findat$P9),mean(findat$P10)))
miwd=1/num*mean(c(var(findat$P1),var(findat$P2),var(findat$P3),var(findat$P4),var(findat$P5),var(findat$P6),var(findat$P7),var(findat$P8),var(findat$P9),var(findat$P10)))
mithetahat=var(c(mean(findat$P1),mean(findat$P2),mean(findat$P3),mean(findat$P4),mean(findat$P5),mean(findat$P6),mean(findat$P7),mean(findat$P8),mean(findat$P9),mean(findat$P10)))
mivar=miwd+(D+1)/D*mithetahat
miqnum=qt(0.975,(D-1)*(1+1/(D+1)*miwd/mithetahat)^2)
miub=mimean+miqnum*sqrt(mivar)
milb=mimean-miqnum*sqrt(mivar)

#geq1
ccmean1=mean(ccbac>1)
ccvar=ccmean1*(1-ccmean1)
ccub1=ccmean1+qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))
cclb1=ccmean1-qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))

tmpmimean1=c(mean(findat$P1>1),mean(findat$P2>1),mean(findat$P3>1),mean(findat$P4>1),mean(findat$P5>1),mean(findat$P6>1),mean(findat$P7>1),mean(findat$P8>1),mean(findat$P9>1),mean(findat$P10>1))
tmpmivar=tmpmimean1*(1-tmpmimean1)
mimean1=mean(tmpmimean1)
miwd=1/num*mean(tmpmivar)
mithetahat=var(tmpmimean1)
mivar=miwd+(D+1)/D*mithetahat
miqnum=qt(0.975,(D-1)*(1+1/(D+1)*miwd/mithetahat)^2)
miub1=mimean1+miqnum*sqrt(mivar)
milb1=mimean1-miqnum*sqrt(mivar)

#geq10
ccmean10=mean(ccbac>10)
ccvar=ccmean10*(1-ccmean10)
ccub10=ccmean10+qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))
cclb10=ccmean10-qt(0.975,length(ccbac)-1)*sqrt(ccvar/length(ccbac))

tmpmimean10=c(mean(findat$P1>10),mean(findat$P2>10),mean(findat$P3>10),mean(findat$P4>10),mean(findat$P5>10),mean(findat$P6>10),mean(findat$P7>10),mean(findat$P8>10),mean(findat$P9>10),mean(findat$P10>10))
tmpmivar=tmpmimean10*(1-tmpmimean10)
mimean10=mean(tmpmimean10)
miwd=1/num*mean(tmpmivar)
mithetahat=var(tmpmimean10)
mivar=miwd+(D+1)/D*mithetahat
miqnum=qt(0.975,(D-1)*(1+1/(D+1)*miwd/mithetahat)^2)
miub10=mimean10+miqnum*sqrt(mivar)
milb10=mimean10-miqnum*sqrt(mivar)

psppmean=mean(tmppsppmean)
psppwd=1/num*mean(tmppsppvar)
psppthetahat=var(tmppsppmean)
psppvar=psppwd+(D+1)/D*psppthetahat
psppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psppwd/psppthetahat)^2)
psppub=psppmean+psppqnum*sqrt(psppvar)
pspplb=psppmean-psppqnum*sqrt(psppvar)

aipwtmean=mean(tmpaipwtmean)
aipwtwd=1/num*mean(tmpaipwtvar)
aipwtthetahat=var(tmpaipwtmean)
aipwtvar=aipwtwd+(D+1)/D*aipwtthetahat
aipwtqnum=qt(0.975,(D-1)*(1+1/(D+1)*aipwtwd/aipwtthetahat)^2)
aipwtub=aipwtmean+aipwtqnum*sqrt(aipwtvar)
aipwtlb=aipwtmean-aipwtqnum*sqrt(aipwtvar)

psbppmean=mean(tmppsbppmean)
psbppwd=1/num*mean(tmppsbppvar)
psbppthetahat=var(tmppsbppmean)
psbppvar=psbppwd+(D+1)/D*psbppthetahat
psbppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psbppwd/psbppthetahat)^2)
psbppub=psbppmean+psbppqnum*sqrt(psbppvar)
psbpplb=psbppmean-psbppqnum*sqrt(psbppvar)

bartpsmean=mean(tmpbartpsmean)
bartpswd=1/num*mean(tmpbartpsvar)
bartpsthetahat=var(tmpbartpsmean)
bartpsvar=bartpswd+(D+1)/D*bartpsthetahat
bartpsqnum=qt(0.975,(D-1)*(1+1/(D+1)*bartpswd/bartpsthetahat)^2)
bartpsub=bartpsmean+bartpsqnum*sqrt(bartpsvar)
bartpslb=bartpsmean-bartpsqnum*sqrt(bartpsvar)

psppmean1=mean(tmppsppmean1)
psppwd=1/num*mean(tmppsppvar)
psppthetahat=var(tmppsppmean1)
psppvar=psppwd+(D+1)/D*psppthetahat
psppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psppwd/psppthetahat)^2)
psppub1=psppmean1+psppqnum*sqrt(psppvar)
pspplb1=psppmean1-psppqnum*sqrt(psppvar)

aipwtmean1=mean(tmpaipwtmean1)
aipwtwd=1/num*mean(tmpaipwtvar)
aipwtthetahat=var(tmpaipwtmean1)
aipwtvar=aipwtwd+(D+1)/D*aipwtthetahat
aipwtqnum=qt(0.975,(D-1)*(1+1/(D+1)*aipwtwd/aipwtthetahat)^2)
aipwtub1=aipwtmean1+aipwtqnum*sqrt(aipwtvar)
aipwtlb1=aipwtmean1-aipwtqnum*sqrt(aipwtvar)

psbppmean1=mean(tmppsbppmean1)
psbppwd=1/num*mean(tmppsbppvar)
psbppthetahat=var(tmppsbppmean1)
psbppvar=psbppwd+(D+1)/D*psbppthetahat
psbppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psbppwd/psbppthetahat)^2)
psbppub1=psbppmean1+psbppqnum*sqrt(psbppvar)
psbpplb1=psbppmean1-psbppqnum*sqrt(psbppvar)

bartpsmean1=mean(tmpbartpsmean1)
bartpswd=1/num*mean(tmpbartpsvar)
bartpsthetahat=var(tmpbartpsmean1)
bartpsvar=bartpswd+(D+1)/D*bartpsthetahat
bartpsqnum=qt(0.975,(D-1)*(1+1/(D+1)*bartpswd/bartpsthetahat)^2)
bartpsub1=bartpsmean1+bartpsqnum*sqrt(bartpsvar)
bartpslb1=bartpsmean1-bartpsqnum*sqrt(bartpsvar)

psppmean10=mean(tmppsppmean10)
psppwd=1/num*mean(tmppsppvar)
psppthetahat=var(tmppsppmean10)
psppvar=psppwd+(D+1)/D*psppthetahat
psppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psppwd/psppthetahat)^2)
psppub10=psppmean10+psppqnum*sqrt(psppvar)
pspplb10=psppmean10-psppqnum*sqrt(psppvar)

aipwtmean10=mean(tmpaipwtmean10)
aipwtwd=1/num*mean(tmpaipwtvar)
aipwtthetahat=var(tmpaipwtmean10)
aipwtvar=aipwtwd+(D+1)/D*aipwtthetahat
aipwtqnum=qt(0.975,(D-1)*(1+1/(D+1)*aipwtwd/aipwtthetahat)^2)
aipwtub10=aipwtmean10+aipwtqnum*sqrt(aipwtvar)
aipwtlb10=aipwtmean10-aipwtqnum*sqrt(aipwtvar)

psbppmean10=mean(tmppsbppmean10)
psbppwd=1/num*mean(tmppsbppvar)
psbppthetahat=var(tmppsbppmean10)
psbppvar=psbppwd+(D+1)/D*psbppthetahat
psbppqnum=qt(0.975,(D-1)*(1+1/(D+1)*psbppwd/psbppthetahat)^2)
psbppub10=psbppmean10+psbppqnum*sqrt(psbppvar)
psbpplb10=psbppmean10-psbppqnum*sqrt(psbppvar)

bartpsmean10=mean(tmpbartpsmean10)
bartpswd=1/num*mean(tmpbartpsvar)
bartpsthetahat=var(tmpbartpsmean10)
bartpsvar=bartpswd+(D+1)/D*bartpsthetahat
bartpsqnum=qt(0.975,(D-1)*(1+1/(D+1)*bartpswd/bartpsthetahat)^2)
bartpsub10=bartpsmean10+bartpsqnum*sqrt(bartpsvar)
bartpslb10=bartpsmean10-bartpsqnum*sqrt(bartpsvar)

print(rbind(c("CC",round(ccmean,2),paste0("(",round(cclb,2),", ",round(ccub,2),")"),round(ccmean1,2),paste0("(",round(cclb1,2),", ",round(ccub1,2),")"),round(ccmean10,2),paste0("(",round(cclb10,2),", ",round(ccub10,2),")")),c("MI",round(mimean,2),paste0("(",round(milb,2),", ",round(miub,2),")"),round(mimean1,2),paste0("(",round(milb1,2),", ",round(miub1,2),")"),round(mimean10,2),paste0("(",round(milb10,2),", ",round(miub10,2),")")),c("PSPP",round(psppmean,2),paste0("(",round(pspplb,2),", ",round(psppub,2),")"),round(psppmean1,2),paste0("(",round(pspplb1,2),", ",round(psppub1,2),")"),round(psppmean10,2),paste0("(",round(pspplb10,2),", ",round(psppub10,2),")")),c("AIPWT",round(aipwtmean,2),paste0("(",round(aipwtlb,2),", ",round(aipwtub,2),")"),round(aipwtmean1,2),paste0("(",round(aipwtlb1,2),", ",round(aipwtub1,2),")"),round(aipwtmean10,2),paste0("(",round(aipwtlb10,2),", ",round(aipwtub10,2),")")),c("PSBPP",round(psbppmean,2),paste0("(",round(psbpplb,2),", ",round(psbppub,2),")"),round(psbppmean1,2),paste0("(",round(psbpplb1,2),", ",round(psbppub1,2),")"),round(psbppmean10,2),paste0("(",round(psbpplb10,2),", ",round(psbppub10,2),")")),c("BARTps",round(bartpsmean,2),paste0("(",round(bartpslb,2),", ",round(bartpsub,2),")"),round(bartpsmean1,2),paste0("(",round(bartpslb1,2),", ",round(bartpsub1,2),")"),round(bartpsmean10,2),paste0("(",round(bartpslb10,2),", ",round(bartpsub10,2),")"))))
