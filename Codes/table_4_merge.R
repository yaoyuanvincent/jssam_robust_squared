#Takes in results from table_4.R to produce Table 4.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

nsim=250

vec_CC_mean=numeric(nsim)
vec_MI_mean=numeric(nsim)
vec_PSPP_mean=numeric(nsim)
vec_AIPWT_mean=numeric(nsim)
vec_PSBPP_mean=numeric(nsim)
vec_BARTps_mean=numeric(nsim)
vec_CC_OR_IS_any1=numeric(nsim)
vec_MI_OR_IS_any1=numeric(nsim)
vec_PSPP_OR_IS_any1=numeric(nsim)
vec_AIPWT_OR_IS_any1=numeric(nsim)
vec_PSBPP_OR_IS_any1=numeric(nsim)
vec_BARTps_OR_IS_any1=numeric(nsim)
vec_CC_OR_IS_sev1=numeric(nsim)
vec_MI_OR_IS_sev1=numeric(nsim)
vec_PSPP_OR_IS_sev1=numeric(nsim)
vec_AIPWT_OR_IS_sev1=numeric(nsim)
vec_PSBPP_OR_IS_sev1=numeric(nsim)
vec_BARTps_OR_IS_sev1=numeric(nsim)
vec_CC_OR_IS_any2=numeric(nsim)
vec_MI_OR_IS_any2=numeric(nsim)
vec_PSPP_OR_IS_any2=numeric(nsim)
vec_AIPWT_OR_IS_any2=numeric(nsim)
vec_PSBPP_OR_IS_any2=numeric(nsim)
vec_BARTps_OR_IS_any2=numeric(nsim)
vec_CC_OR_IS_sev2=numeric(nsim)
vec_MI_OR_IS_sev2=numeric(nsim)
vec_PSPP_OR_IS_sev2=numeric(nsim)
vec_AIPWT_OR_IS_sev2=numeric(nsim)
vec_PSBPP_OR_IS_sev2=numeric(nsim)
vec_BARTps_OR_IS_sev2=numeric(nsim)

for(loop in 1:nsim){
	load(paste0("table_4_",loop,".RData"))
	vec_CC_mean[loop]=CC_mean
	vec_MI_mean[loop]=MI_mean
	vec_PSPP_mean[loop]=PSPP_mean
	vec_AIPWT_mean[loop]=AIPWT_mean
	vec_PSBPP_mean[loop]=PSBPP_mean
	vec_BARTps_mean[loop]=BARTps_mean
	vec_CC_OR_IS_any1[loop]=CC_OR_IS_any1
	vec_MI_OR_IS_any1[loop]=MI_OR_IS_any1
	vec_PSPP_OR_IS_any1[loop]=PSPP_OR_IS_any1
	vec_AIPWT_OR_IS_any1[loop]=AIPWT_OR_IS_any1
	vec_PSBPP_OR_IS_any1[loop]=PSBPP_OR_IS_any1
	vec_BARTps_OR_IS_any1[loop]=BARTps_OR_IS_any1
	vec_CC_OR_IS_sev1[loop]=CC_OR_IS_sev1
	vec_MI_OR_IS_sev1[loop]=MI_OR_IS_sev1
	vec_PSPP_OR_IS_sev1[loop]=PSPP_OR_IS_sev1
	vec_AIPWT_OR_IS_sev1[loop]=AIPWT_OR_IS_sev1
	vec_PSBPP_OR_IS_sev1[loop]=PSBPP_OR_IS_sev1
	vec_BARTps_OR_IS_sev1[loop]=BARTps_OR_IS_sev1
	vec_CC_OR_IS_any2[loop]=CC_OR_IS_any2
	vec_MI_OR_IS_any2[loop]=MI_OR_IS_any2
	vec_PSPP_OR_IS_any2[loop]=PSPP_OR_IS_any2
	vec_AIPWT_OR_IS_any2[loop]=AIPWT_OR_IS_any2
	vec_PSBPP_OR_IS_any2[loop]=PSBPP_OR_IS_any2
	vec_BARTps_OR_IS_any2[loop]=BARTps_OR_IS_any2
	vec_CC_OR_IS_sev2[loop]=CC_OR_IS_sev2
	vec_MI_OR_IS_sev2[loop]=MI_OR_IS_sev2
	vec_PSPP_OR_IS_sev2[loop]=PSPP_OR_IS_sev2
	vec_AIPWT_OR_IS_sev2[loop]=AIPWT_OR_IS_sev2
	vec_PSBPP_OR_IS_sev2[loop]=PSBPP_OR_IS_sev2
	vec_BARTps_OR_IS_sev2[loop]=BARTps_OR_IS_sev2
}

Bt1=50
sub_CC_mean=numeric(Bt1)
sub_MI_mean=numeric(Bt1)
sub_PSPP_mean=numeric(Bt1)
sub_AIPWT_mean=numeric(Bt1)
sub_PSBPP_mean=numeric(Bt1)
sub_BARTps_mean=numeric(Bt1)
sub_CC_OR_IS_any1=numeric(Bt1)
sub_MI_OR_IS_any1=numeric(Bt1)
sub_PSPP_OR_IS_any1=numeric(Bt1)
sub_AIPWT_OR_IS_any1=numeric(Bt1)
sub_PSBPP_OR_IS_any1=numeric(Bt1)
sub_BARTps_OR_IS_any1=numeric(Bt1)
sub_CC_OR_IS_sev1=numeric(Bt1)
sub_MI_OR_IS_sev1=numeric(Bt1)
sub_PSPP_OR_IS_sev1=numeric(Bt1)
sub_AIPWT_OR_IS_sev1=numeric(Bt1)
sub_PSBPP_OR_IS_sev1=numeric(Bt1)
sub_BARTps_OR_IS_sev1=numeric(Bt1)
sub_CC_OR_IS_any2=numeric(Bt1)
sub_MI_OR_IS_any2=numeric(Bt1)
sub_PSPP_OR_IS_any2=numeric(Bt1)
sub_AIPWT_OR_IS_any2=numeric(Bt1)
sub_PSBPP_OR_IS_any2=numeric(Bt1)
sub_BARTps_OR_IS_any2=numeric(Bt1)
sub_CC_OR_IS_sev2=numeric(Bt1)
sub_MI_OR_IS_sev2=numeric(Bt1)
sub_PSPP_OR_IS_sev2=numeric(Bt1)
sub_AIPWT_OR_IS_sev2=numeric(Bt1)
sub_PSBPP_OR_IS_sev2=numeric(Bt1)
sub_BARTps_OR_IS_sev2=numeric(Bt1)

for(j in 1:Bt1){
	end=j*5
	start=end-4
	sub_CC_mean[j]=mean(vec_CC_mean[start:end])
	sub_MI_mean[j]=mean(vec_MI_mean[start:end])
	sub_PSPP_mean[j]=mean(vec_PSPP_mean[start:end])
	sub_AIPWT_mean[j]=mean(vec_AIPWT_mean[start:end])
	sub_PSBPP_mean[j]=mean(vec_PSBPP_mean[start:end])
	sub_BARTps_mean[j]=mean(vec_BARTps_mean[start:end])
	sub_CC_OR_IS_any1[j]=mean(vec_CC_OR_IS_any1[start:end])
	sub_MI_OR_IS_any1[j]=mean(vec_MI_OR_IS_any1[start:end])
	sub_PSPP_OR_IS_any1[j]=mean(vec_PSPP_OR_IS_any1[start:end])
	sub_AIPWT_OR_IS_any1[j]=mean(vec_AIPWT_OR_IS_any1[start:end])
	sub_PSBPP_OR_IS_any1[j]=mean(vec_PSBPP_OR_IS_any1[start:end])
	sub_BARTps_OR_IS_any1[j]=mean(vec_BARTps_OR_IS_any1[start:end])
	sub_CC_OR_IS_sev1[j]=mean(vec_CC_OR_IS_sev1[start:end])
	sub_MI_OR_IS_sev1[j]=mean(vec_MI_OR_IS_sev1[start:end])
	sub_PSPP_OR_IS_sev1[j]=mean(vec_PSPP_OR_IS_sev1[start:end])
	sub_AIPWT_OR_IS_sev1[j]=mean(vec_AIPWT_OR_IS_sev1[start:end])
	sub_PSBPP_OR_IS_sev1[j]=mean(vec_PSBPP_OR_IS_sev1[start:end])
	sub_BARTps_OR_IS_sev1[j]=mean(vec_BARTps_OR_IS_sev1[start:end])
	sub_CC_OR_IS_any2[j]=mean(vec_CC_OR_IS_any2[start:end])
	sub_MI_OR_IS_any2[j]=mean(vec_MI_OR_IS_any2[start:end])
	sub_PSPP_OR_IS_any2[j]=mean(vec_PSPP_OR_IS_any2[start:end])
	sub_AIPWT_OR_IS_any2[j]=mean(vec_AIPWT_OR_IS_any2[start:end])
	sub_PSBPP_OR_IS_any2[j]=mean(vec_PSBPP_OR_IS_any2[start:end])
	sub_BARTps_OR_IS_any2[j]=mean(vec_BARTps_OR_IS_any2[start:end])
	sub_CC_OR_IS_sev2[j]=mean(vec_CC_OR_IS_sev2[start:end])
	sub_MI_OR_IS_sev2[j]=mean(vec_MI_OR_IS_sev2[start:end])
	sub_PSPP_OR_IS_sev2[j]=mean(vec_PSPP_OR_IS_sev2[start:end])
	sub_AIPWT_OR_IS_sev2[j]=mean(vec_AIPWT_OR_IS_sev2[start:end])
	sub_PSBPP_OR_IS_sev2[j]=mean(vec_PSBPP_OR_IS_sev2[start:end])
	sub_BARTps_OR_IS_sev2[j]=mean(vec_BARTps_OR_IS_sev2[start:end])
}

res1=c(round(mean(sub_CC_mean),2),round(mean(sub_CC_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_mean)),2),round(mean(sub_CC_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_mean)),2),round(mean(sub_CC_OR_IS_any1),2),round(mean(sub_CC_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_any1)),2),round(mean(sub_CC_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_any1)),2),round(mean(sub_CC_OR_IS_any2),2),round(mean(sub_CC_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_any2)),2),round(mean(sub_CC_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_any2)),2),round(mean(sub_CC_OR_IS_sev1),2),round(mean(sub_CC_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_sev1)),2),round(mean(sub_CC_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_sev1)),2),round(mean(sub_CC_OR_IS_sev2),2),round(mean(sub_CC_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_sev2)),2),round(mean(sub_CC_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_CC_OR_IS_sev2)),2))
res2=c(round(mean(sub_MI_mean),2),round(mean(sub_MI_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_mean)),2),round(mean(sub_MI_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_mean)),2),round(mean(sub_MI_OR_IS_any1),2),round(mean(sub_MI_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_any1)),2),round(mean(sub_MI_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_any1)),2),round(mean(sub_MI_OR_IS_any2),2),round(mean(sub_MI_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_any2)),2),round(mean(sub_MI_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_any2)),2),round(mean(sub_MI_OR_IS_sev1),2),round(mean(sub_MI_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_sev1)),2),round(mean(sub_MI_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_sev1)),2),round(mean(sub_MI_OR_IS_sev2),2),round(mean(sub_MI_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_sev2)),2),round(mean(sub_MI_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_MI_OR_IS_sev2)),2))
res3=c(round(mean(sub_PSPP_mean),2),round(mean(sub_PSPP_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_mean)),2),round(mean(sub_PSPP_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_mean)),2),round(mean(sub_PSPP_OR_IS_any1),2),round(mean(sub_PSPP_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_any1)),2),round(mean(sub_PSPP_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_any1)),2),round(mean(sub_PSPP_OR_IS_any2),2),round(mean(sub_PSPP_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_any2)),2),round(mean(sub_PSPP_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_any2)),2),round(mean(sub_PSPP_OR_IS_sev1),2),round(mean(sub_PSPP_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_sev1)),2),round(mean(sub_PSPP_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_sev1)),2),round(mean(sub_PSPP_OR_IS_sev2),2),round(mean(sub_PSPP_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_sev2)),2),round(mean(sub_PSPP_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSPP_OR_IS_sev2)),2))
res4=c(round(mean(sub_AIPWT_mean),2),round(mean(sub_AIPWT_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_mean)),2),round(mean(sub_AIPWT_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_mean)),2),round(mean(sub_AIPWT_OR_IS_any1),2),round(mean(sub_AIPWT_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_any1)),2),round(mean(sub_AIPWT_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_any1)),2),round(mean(sub_AIPWT_OR_IS_any2),2),round(mean(sub_AIPWT_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_any2)),2),round(mean(sub_AIPWT_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_any2)),2),round(mean(sub_AIPWT_OR_IS_sev1),2),round(mean(sub_AIPWT_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_sev1)),2),round(mean(sub_AIPWT_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_sev1)),2),round(mean(sub_AIPWT_OR_IS_sev2),2),round(mean(sub_AIPWT_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_sev2)),2),round(mean(sub_AIPWT_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_AIPWT_OR_IS_sev2)),2))
res5=c(round(mean(sub_PSBPP_mean),2),round(mean(sub_PSBPP_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_mean)),2),round(mean(sub_PSBPP_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_mean)),2),round(mean(sub_PSBPP_OR_IS_any1),2),round(mean(sub_PSBPP_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_any1)),2),round(mean(sub_PSBPP_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_any1)),2),round(mean(sub_PSBPP_OR_IS_any2),2),round(mean(sub_PSBPP_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_any2)),2),round(mean(sub_PSBPP_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_any2)),2),round(mean(sub_PSBPP_OR_IS_sev1),2),round(mean(sub_PSBPP_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_sev1)),2),round(mean(sub_PSBPP_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_sev1)),2),round(mean(sub_PSBPP_OR_IS_sev2),2),round(mean(sub_PSBPP_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_sev2)),2),round(mean(sub_PSBPP_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_PSBPP_OR_IS_sev2)),2))
res6=c(round(mean(sub_BARTps_mean),2),round(mean(sub_BARTps_mean)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_mean)),2),round(mean(sub_BARTps_mean)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_mean)),2),round(mean(sub_BARTps_OR_IS_any1),2),round(mean(sub_BARTps_OR_IS_any1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_any1)),2),round(mean(sub_BARTps_OR_IS_any1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_any1)),2),round(mean(sub_BARTps_OR_IS_any2),2),round(mean(sub_BARTps_OR_IS_any2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_any2)),2),round(mean(sub_BARTps_OR_IS_any2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_any2)),2),round(mean(sub_BARTps_OR_IS_sev1),2),round(mean(sub_BARTps_OR_IS_sev1)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_sev1)),2),round(mean(sub_BARTps_OR_IS_sev1)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_sev1)),2),round(mean(sub_BARTps_OR_IS_sev2),2),round(mean(sub_BARTps_OR_IS_sev2)-qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_sev2)),2),round(mean(sub_BARTps_OR_IS_sev2)+qt(0.975,Bt1-1)*sqrt((1+1/Bt1)*var(sub_BARTps_OR_IS_sev2)),2))

write.csv(rbind(res1,res2,res3,res4,res5,res6),file="deltav.csv")