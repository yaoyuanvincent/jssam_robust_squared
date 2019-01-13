#Takes in results from table_1.R, table_2.R, or table_3.R to produce Tables 1, 2, or 3.

#Author: Yaoyuan Vincent, Tan
#Date created: 11 Jan 2019
#Date updated: 11 Jan 2019

rm(list=ls())

nsim=500

bd_param=matrix(NA,nrow=nsim,ncol=4)
cc_param=matrix(NA,nrow=nsim,ncol=4)
mi_param=matrix(NA,nrow=nsim,ncol=4)
np_param=matrix(NA,nrow=nsim,ncol=4)
pspp_param=matrix(NA,nrow=nsim,ncol=4)
psbpp_param=matrix(NA,nrow=nsim,ncol=4)
rrz_param=matrix(NA,nrow=nsim,ncol=4)
rrz_bart_param=matrix(NA,nrow=nsim,ncol=4)
bart_param=matrix(NA,nrow=nsim,ncol=4)
bartps_param=matrix(NA,nrow=nsim,ncol=4)
mr_param=matrix(NA,nrow=nsim,ncol=4)

bd_cover=matrix(NA,nrow=nsim,ncol=4)
cc_cover=matrix(NA,nrow=nsim,ncol=4)
mi_cover=matrix(NA,nrow=nsim,ncol=4)
np_cover=matrix(NA,nrow=nsim,ncol=4)
pspp_cover=matrix(NA,nrow=nsim,ncol=4)
psbpp_cover=matrix(NA,nrow=nsim,ncol=4)
rrz_cover=matrix(NA,nrow=nsim,ncol=4)
rrz_bart_cover=matrix(NA,nrow=nsim,ncol=4)
bart_cover=matrix(NA,nrow=nsim,ncol=4)
bartps_cover=matrix(NA,nrow=nsim,ncol=4)
mr_cover=matrix(NA,nrow=nsim,ncol=4)

bd_length=matrix(NA,nrow=nsim,ncol=4)
cc_length=matrix(NA,nrow=nsim,ncol=4)
mi_length=matrix(NA,nrow=nsim,ncol=4)
np_length=matrix(NA,nrow=nsim,ncol=4)
pspp_length=matrix(NA,nrow=nsim,ncol=4)
psbpp_length=matrix(NA,nrow=nsim,ncol=4)
rrz_length=matrix(NA,nrow=nsim,ncol=4)
rrz_bart_length=matrix(NA,nrow=nsim,ncol=4)
bart_length=matrix(NA,nrow=nsim,ncol=4)
bartps_length=matrix(NA,nrow=nsim,ncol=4)
mr_length=matrix(NA,nrow=nsim,ncol=4)

for(i in 1:nsim){
	load(paste0("table_1_",i,".RData"))
	#load(paste0("table_2_",i,".RData"))
	#load(paste0("table_3_",i,".RData"))

	param=10 #For Tables 1 and 2
	#param=210 #For Table 3

	bd_param[i,]=est_mean_bd
	cc_param[i,]=est_mean_cc
	mi_param[i,]=est_mean_mi
	np_param[i,c(1,3)]=mean_np[1]
	np_param[i,c(2,4)]=mean_np[2]	
	pspp_param[i,]=est_mean_pspp
	psbpp_param[i,]=est_mean_psbpp
	rrz_param[i,]=est_mean_rrz
	rrz_bart_param[i,]=est_mean_rrz_bart
	bart_param[i,]=est_mean_bart
	bartps_param[i,]=est_mean_bartps
	mr_param[i,]=colMeans(est_mean_mr)	

	bd_cover[i,]=ifelse(lower_bd<param&param<upper_bd,1,0)
	cc_cover[i,]=ifelse(lower_cc<param&param<upper_cc,1,0)
	mi_cover[i,]=ifelse(lower_mi<param&param<upper_mi,1,0)
	np_cover[i,c(1,3)]=ifelse(lower_np[1]<param&param<upper_np[1],1,0)
	np_cover[i,c(2,4)]=ifelse(lower_np[2]<param&param<upper_np[2],1,0)
	pspp_cover[i,]=ifelse(lower_pspp<param&param<upper_pspp,1,0)
	psbpp_cover[i,]=ifelse(lower_psbpp<param&param<upper_psbpp,1,0)
	rrz_cover[i,]=ifelse(lower_rrz<param&param<upper_rrz,1,0)
	rrz_bart_cover[i,]=ifelse(lower_rrz_bart<param&param<upper_rrz_bart,1,0)
	bart_cover[i,]=ifelse(lower_bart<param&param<upper_bart,1,0)
	bartps_cover[i,]=ifelse(lower_bartps<param&param<upper_bartps,1,0)	
	mr_cover[i,]=ifelse(lower_mr<param&param<upper_mr,1,0)

	bd_length[i,]=upper_bd-lower_bd
	cc_length[i,]=upper_cc-lower_cc
	mi_length[i,]=upper_mi-lower_mi
	np_length[i,c(1,3)]=upper_np[1]-lower_np[1]
	np_length[i,c(2,4)]=upper_np[2]-lower_np[2]
	pspp_length[i,]=upper_pspp-lower_pspp
	psbpp_length[i,]=upper_psbpp-lower_psbpp
	rrz_length[i,]=upper_rrz-lower_rrz
	rrz_bart_length[i,]=upper_rrz_bart-lower_rrz_bart
	bart_length[i,]=upper_bart-lower_bart
	bartps_length[i,]=upper_bartps-lower_bartps
	mr_length[i,]=upper_mr-lower_mr

}

bias_bd=numeric(4)
bias_cc=numeric(4)
bias_mi=numeric(4)
bias_np=numeric(4)
bias_pspp=numeric(4)
bias_psbpp=numeric(4)
bias_rrz=numeric(4)
bias_rrz_bart=numeric(4)
bias_bart=numeric(4)
bias_bartps=numeric(4)
bias_mr=numeric(4)

rmse_bd=numeric(4)
rmse_cc=numeric(4)
rmse_mi=numeric(4)
rmse_np=numeric(4)
rmse_pspp=numeric(4)
rmse_psbpp=numeric(4)
rmse_rrz=numeric(4)
rmse_rrz_bart=numeric(4)
rmse_bart=numeric(4)
rmse_bartps=numeric(4)
rmse_mr=numeric(4)

cover95_bd=numeric(4)
cover95_cc=numeric(4)
cover95_mi=numeric(4)
cover95_np=numeric(4)
cover95_pspp=numeric(4)
cover95_psbpp=numeric(4)
cover95_rrz=numeric(4)
cover95_rrz_bart=numeric(4)
cover95_bart=numeric(4)
cover95_bartps=numeric(4)
cover95_mr=numeric(4)

avglen_bd=numeric(4)
avglen_cc=numeric(4)
avglen_mi=numeric(4)
avglen_np=numeric(4)
avglen_pspp=numeric(4)
avglen_psbpp=numeric(4)
avglen_rrz=numeric(4)
avglen_rrz_bart=numeric(4)
avglen_bart=numeric(4)
avglen_bartps=numeric(4)
avglen_mr=numeric(4)

for(i in 1:4){
	bias_bd[i]=mean(bd_param[,i])-param
	bias_cc[i]=mean(cc_param[,i])-param
	bias_mi[i]=mean(mi_param[,i])-param
	bias_np[i]=mean(np_param[,i])-param
	bias_pspp[i]=mean(pspp_param[,i])-param
	bias_psbpp[i]=mean(psbpp_param[,i])-param
	bias_rrz[i]=mean(rrz_param[,i])-param
	bias_rrz_bart[i]=mean(rrz_bart_param[,i])-param
	bias_bart[i]=mean(bart_param[,i])-param
	bias_bartps[i]=mean(bartps_param[,i])-param
	bias_mr[i]=mean(mr_param[,i])-param

	rmse_bd[i]=sqrt(sum((bd_param[,i]-param)^2)/nsim)
	rmse_cc[i]=sqrt(sum((cc_param[,i]-param)^2)/nsim)
	rmse_mi[i]=sqrt(sum((mi_param[,i]-param)^2)/nsim)
	rmse_np[i]=sqrt(sum((np_param[,i]-param)^2)/nsim)
	rmse_pspp[i]=sqrt(sum((pspp_param[,i]-param)^2)/nsim)
	rmse_psbpp[i]=sqrt(sum((psbpp_param[,i]-param)^2)/nsim)
	rmse_rrz[i]=sqrt(sum((rrz_param[,i]-param)^2)/nsim)
	rmse_rrz_bart[i]=sqrt(sum((rrz_bart_param[,i]-param)^2)/nsim)
	rmse_bart[i]=sqrt(sum((bart_param[,i]-param)^2)/nsim)
	rmse_bartps[i]=sqrt(sum((bartps_param[,i]-param)^2)/nsim)	
	rmse_mr[i]=sqrt(sum((mr_param[,i]-param)^2)/nsim)

	cover95_bd[i]=mean(bd_cover[,i])*100
	cover95_cc[i]=mean(cc_cover[,i])*100
	cover95_mi[i]=mean(mi_cover[,i])*100
	cover95_np[i]=mean(np_cover[,i])*100
	cover95_pspp[i]=mean(pspp_cover[,i])*100
	cover95_psbpp[i]=mean(psbpp_cover[,i])*100
	cover95_rrz[i]=mean(rrz_cover[,i])*100
	cover95_rrz_bart[i]=mean(rrz_bart_cover[,i])*100
	cover95_bart[i]=mean(bart_cover[,i])*100
	cover95_bartps[i]=mean(bartps_cover[,i])*100
	cover95_mr[i]=mean(mr_cover[,i])*100

	avglen_bd[i]=mean(bd_length[,i])
	avglen_cc[i]=mean(cc_length[,i])
	avglen_mi[i]=mean(mi_length[,i])
	avglen_np[i]=mean(np_length[,i])
	avglen_pspp[i]=mean(pspp_length[,i])
	avglen_psbpp[i]=mean(psbpp_length[,i])
	avglen_rrz[i]=mean(rrz_length[,i])
	avglen_rrz_bart[i]=mean(rrz_bart_length[,i])
	avglen_bart[i]=mean(bart_length[,i])
	avglen_bartps[i]=mean(bartps_length[,i])
	avglen_mr[i]=mean(mr_length[,i])

}

res_bd=c(round(bias_bd[1],2),round(rmse_bd[1],2),round(cover95_bd[1],2),round(avglen_bd[1],2),round(bias_bd[2],2),round(rmse_bd[2],2),round(cover95_bd[2],2),round(avglen_bd[2],2),round(bias_bd[3],2),round(rmse_bd[3],2),round(cover95_bd[3],2),round(avglen_bd[3],2),round(bias_bd[4],2),round(rmse_bd[4],2),round(cover95_bd[4],2),round(avglen_bd[4],2))
res_cc=c(round(bias_cc[1],2),round(rmse_cc[1],2),round(cover95_cc[1],2),round(avglen_cc[1],2),round(bias_cc[2],2),round(rmse_cc[2],2),round(cover95_cc[2],2),round(avglen_cc[2],2),round(bias_cc[3],2),round(rmse_cc[3],2),round(cover95_cc[3],2),round(avglen_cc[3],2),round(bias_cc[4],2),round(rmse_cc[4],2),round(cover95_cc[4],2),round(avglen_cc[4],2))
res_mi=c(round(bias_mi[1],2),round(rmse_mi[1],2),round(cover95_mi[1],2),round(avglen_mi[1],2),round(bias_mi[2],2),round(rmse_mi[2],2),round(cover95_mi[2],2),round(avglen_mi[2],2),round(bias_mi[3],2),round(rmse_mi[3],2),round(cover95_mi[3],2),round(avglen_mi[3],2),round(bias_mi[4],2),round(rmse_mi[4],2),round(cover95_mi[4],2),round(avglen_mi[4],2))
res_np=c(round(bias_np[1],2),round(rmse_np[1],2),round(cover95_np[1],2),round(avglen_np[1],2),round(bias_np[2],2),round(rmse_np[2],2),round(cover95_np[2],2),round(avglen_np[2],2),round(bias_np[3],2),round(rmse_np[3],2),round(cover95_np[3],2),round(avglen_np[3],2),round(bias_np[4],2),round(rmse_np[4],2),round(cover95_np[4],2),round(avglen_np[4],2))
res_pspp=c(round(bias_pspp[1],2),round(rmse_pspp[1],2),round(cover95_pspp[1],2),round(avglen_pspp[1],2),round(bias_pspp[2],2),round(rmse_pspp[2],2),round(cover95_pspp[2],2),round(avglen_pspp[2],2),round(bias_pspp[3],2),round(rmse_pspp[3],2),round(cover95_pspp[3],2),round(avglen_pspp[3],2),round(bias_pspp[4],2),round(rmse_pspp[4],2),round(cover95_pspp[4],2),round(avglen_pspp[4],2))
res_psbpp=c(round(bias_psbpp[1],2),round(rmse_psbpp[1],2),round(cover95_psbpp[1],2),round(avglen_psbpp[1],2),round(bias_psbpp[2],2),round(rmse_psbpp[2],2),round(cover95_psbpp[2],2),round(avglen_psbpp[2],2),round(bias_psbpp[3],2),round(rmse_psbpp[3],2),round(cover95_psbpp[3],2),round(avglen_psbpp[3],2),round(bias_psbpp[4],2),round(rmse_psbpp[4],2),round(cover95_psbpp[4],2),round(avglen_psbpp[4],2))
res_rrz=c(round(bias_rrz[1],2),round(rmse_rrz[1],2),round(cover95_rrz[1],2),round(avglen_rrz[1],2),round(bias_rrz[2],2),round(rmse_rrz[2],2),round(cover95_rrz[2],2),round(avglen_rrz[2],2),round(bias_rrz[3],2),round(rmse_rrz[3],2),round(cover95_rrz[3],2),round(avglen_rrz[3],2),round(bias_rrz[4],2),round(rmse_rrz[4],2),round(cover95_rrz[4],2),round(avglen_rrz[4],2))
res_rrz_bart=c(round(bias_rrz_bart[1],2),round(rmse_rrz_bart[1],2),round(cover95_rrz_bart[1],2),round(avglen_rrz_bart[1],2),round(bias_rrz_bart[2],2),round(rmse_rrz_bart[2],2),round(cover95_rrz_bart[2],2),round(avglen_rrz_bart[2],2),round(bias_rrz_bart[3],2),round(rmse_rrz_bart[3],2),round(cover95_rrz_bart[3],2),round(avglen_rrz_bart[3],2),round(bias_rrz_bart[4],2),round(rmse_rrz_bart[4],2),round(cover95_rrz_bart[4],2),round(avglen_rrz_bart[4],2))
res_bart=c(round(bias_bart[1],2),round(rmse_bart[1],2),round(cover95_bart[1],2),round(avglen_bart[1],2),round(bias_bart[2],2),round(rmse_bart[2],2),round(cover95_bart[2],2),round(avglen_bart[2],2),round(bias_bart[3],2),round(rmse_bart[3],2),round(cover95_bart[3],2),round(avglen_bart[3],2),round(bias_bart[4],2),round(rmse_bart[4],2),round(cover95_bart[4],2),round(avglen_bart[4],2))
res_bartps=c(round(bias_bartps[1],2),round(rmse_bartps[1],2),round(cover95_bartps[1],2),round(avglen_bartps[1],2),round(bias_bartps[2],2),round(rmse_bartps[2],2),round(cover95_bartps[2],2),round(avglen_bartps[2],2),round(bias_bartps[3],2),round(rmse_bartps[3],2),round(cover95_bartps[3],2),round(avglen_bartps[3],2),round(bias_bartps[4],2),round(rmse_bartps[4],2),round(cover95_bartps[4],2),round(avglen_bartps[4],2))
res_mr=c(round(bias_mr[1],2),round(rmse_mr[1],2),round(cover95_mr[1],2),round(avglen_mr[1],2),round(bias_mr[2],2),round(rmse_mr[2],2),round(cover95_mr[2],2),round(avglen_mr[2],2),round(bias_mr[3],2),round(rmse_mr[3],2),round(cover95_mr[3],2),round(avglen_mr[3],2),round(bias_mr[4],2),round(rmse_mr[4],2),round(cover95_mr[4],2),round(avglen_mr[4],2))

write.csv(rbind(res_bd,res_cc,res_mi,res_np,res_pspp,res_psbpp,res_rrz,res_rrz_bart,res_bart,res_bartps,res_mr),file="Table_1.csv")
#write.csv(rbind(res_bd,res_cc,res_mi,res_np,res_pspp,res_psbpp,res_rrz,res_rrz_bart,res_bart,res_bartps,res_mr),file="Table_2.csv")
#write.csv(rbind(res_bd,res_cc,res_mi,res_np,res_pspp,res_psbpp,res_rrz,res_rrz_bart,res_bart,res_bartps,res_mr),file="Table_3.csv")