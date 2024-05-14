# Supporting Material S2
# Analyses for Comparative life-history responses of lacewings to changes in temperature
# created by Maria Paniw
# last modified: 04-11-2023

set.seed(14052024)

#############################   Load necessary packages
library(MCMCglmm)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(hdrcde)
library(coda)

###########################   Data preparation 

df=read.csv("data/LH_Neuroptera_red.csv") %>%
  mutate(temp = as.numeric(temp),
         L_surv = as.numeric(L_surv),
         M_sex_rep_rate_M = as.numeric(M_sex_rep_rate_M))


sub_noNA=na.omit(df[,c('St_ID', 'sp.', 'temp', 'Latitude', "L_surv", "M_sex_rep_rate_M", "Types_of_Lit_Sources")])

# Values necessary to back-transform scaled temperatures for plotting
mean.temp=25.05
sd.temp=3.80

sub_noNA$temp=scale(as.numeric(sub_noNA$temp))
sub_noNA$temp2=sub_noNA$temp^2

sub_noNA$Latitude=scale(as.numeric(sub_noNA$Latitude))

sub_noNA$L_surv=asin(sqrt(sub_noNA$L_surv/100))
sub_noNA=sub_noNA[sub_noNA$M_sex_rep_rate_M>0,]

sub_noNA$M_sex_rep_rate_M=log(sub_noNA$M_sex_rep_rate_M)

nrow(sub_noNA)

###########################   MCMC analyses
prior = list(R = list(V = diag(2)/3, n = 2, nu=0.002),
             G = list(G1 = list(V = diag(2)/3, n = 2, nu=0.002),
                      G2 = list(V = diag(2)/3, n = 2, nu=0.002)))

m1=MCMCglmm(cbind(L_surv,M_sex_rep_rate_M)~trait+trait:temp+trait:temp2 + trait:Latitude + trait:Types_of_Lit_Sources,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 2), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)

m2=MCMCglmm(cbind(L_surv,M_sex_rep_rate_M)~trait+trait:temp+trait:temp2 + trait:Latitude + trait:Types_of_Lit_Sources,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 2), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)

m3=MCMCglmm(cbind(L_surv,M_sex_rep_rate_M)~trait+trait:temp+trait:temp2 + trait:Latitude + trait:Types_of_Lit_Sources,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 2), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)



param.coda.fixed=mcmc.list(list(mcmc(m1$Sol),mcmc(m2$Sol),mcmc(m3$Sol)))

summary(param.coda.fixed)
gelman.diag(param.coda.fixed,multivariate=F)

#Trace plots (to check if chains are well mixed)

par(mar=c(2,2,2,2))
plot(param.coda.fixed,smooth=F) # The different colors indicate different chains

### Plot response as function of temperature
out.mcmc=rbind(m1$Sol,m2$Sol,m3$Sol)

temp.pred=seq(min(sub_noNA$temp),max(sub_noNA$temp),length.out=20)
temp.pred2=temp.pred^2

new.data_1=data.frame(temp=temp.pred,
                      dev=sin(mean(out.mcmc[,1])+mean(out.mcmc[,3])*temp.pred+mean(out.mcmc[,5])*temp.pred2)^2,
                      LB=sin((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,3],0.025)*temp.pred+quantile(out.mcmc[,5],0.025)*temp.pred2))^2,
                      UB=sin((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,3],0.975)*temp.pred+quantile(out.mcmc[,5],0.975)*temp.pred2))^2,
                      stage="L_surv")


new.data_2=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,2])+mean(out.mcmc[,4])*temp.pred+mean(out.mcmc[,6])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,2],0.025)+quantile(out.mcmc[,4],0.025)*temp.pred+quantile(out.mcmc[,6],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,2],0.975)+quantile(out.mcmc[,4],0.975)*temp.pred+quantile(out.mcmc[,6],0.975)*temp.pred2)),
                      stage="M_sex_rep_rate_M")


pred.data=rbind(new.data_1,new.data_2)

pred.data$stage=factor(pred.data$stage,levels=c("L_surv","M_sex_rep_rate_M"))

levels(pred.data$stage) = c("S pupae-adult","# eggs/female")

#observed 

df.obs <- sub_noNA %>%
  pivot_longer(
    cols = "L_surv":"M_sex_rep_rate_M",
    names_to = "stage",
    values_to = "dev"
  )


df.obs$dev[!df.obs$stage%in%"M_sex_rep_rate_M"]=sin(df.obs$dev[!df.obs$stage%in%"M_sex_rep_rate_M"])^2
df.obs$dev[df.obs$stage%in%"M_sex_rep_rate_M"]=exp(df.obs$dev[df.obs$stage%in%"M_sex_rep_rate_M"])


df.obs$stage=factor(df.obs$stage,levels=c("L_surv","M_sex_rep_rate_M"))

levels(df.obs$stage)= c("S pupae-adult", "# eggs/female")

df.obs$sp.=gsub("_"," ", df.obs$sp.)
pred.data$temp=pred.data$temp*sd.temp+mean.temp
df.obs$temp=df.obs$temp*sd.temp+mean.temp

p.temp=ggplot(pred.data, aes(temp, dev))+
  facet_grid(stage~.,scales="free")+
  geom_line() +
  geom_point(data=df.obs,aes(temp, dev, col=sp.),size=2)+
  geom_ribbon(aes(ymin = LB, ymax = UB),alpha=0.1,col=NA) +
  xlab("Temperature (ºC)")+
  ylab("")+
  theme_bw(base_size = 18)+
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic",size=10))+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=16))


p.temp


ggsave(filename = "results/plot_temp_resp_surv_only_V2.pdf",plot=p.temp,width = 12,height = 6)


#covariance
colnames(m1$VCV)=gsub(".units","",colnames(m1$VCV))
colnames(m1$VCV)=gsub("trait","",colnames(m1$VCV))

colnames(m2$VCV)=gsub(".units","",colnames(m2$VCV))
colnames(m2$VCV)=gsub("trait","",colnames(m2$VCV))

colnames(m3$VCV)=gsub(".units","",colnames(m3$VCV))
colnames(m3$VCV)=gsub("trait","",colnames(m3$VCV))

m1.sub.res=m1$VCV[,c(5,6,8)]
m2.sub.res=m2$VCV[,c(5,6,8)]
m3.sub.res=m3$VCV[,c(5,6,8)]
colnames(m1.sub.res)=colnames(m2.sub.res)=colnames(m3.sub.res)=c("S PA - S PA","S PA - Repro","Repro - Repro")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.res),mcmc(m2.sub.res),mcmc(m3.sub.res)))

library(MCMCvis)


MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_surv_repro",wd="results/")

pdf("results/FigS2.6b.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")

dev.off()

#### Species level covariance

m1.sub.sp=m1$VCV[,c(1,2,4)]
m2.sub.sp=m2$VCV[,c(1,2,4)]
m3.sub.sp=m3$VCV[,c(1,2,4)]
colnames(m1.sub.sp)=colnames(m2.sub.sp)=colnames(m3.sub.sp)=c("S PA - S PA","S PA - Repro","Repro - Repro")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.sp),mcmc(m2.sub.sp),mcmc(m3.sub.sp)))

#### species covariance: 

MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_surv_repro_sp",wd="results/")

pdf("results/FigS.2.6a.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random species-specific covariance")

dev.off()


###### Percent variance expalin by species

m1.sub=abs(m1.sub.sp)/(abs(m1.sub.sp)+abs(m1.sub.res))
m2.sub=abs(m2.sub.sp)/(abs(m2.sub.sp)+abs(m2.sub.res))
m3.sub=abs(m3.sub.sp)/(abs(m3.sub.sp)+abs(m3.sub.res))

colnames(m1.sub)=colnames(m2.sub)=colnames(m3.sub)=c("S PA - S PA","S PA - Repro","Repro - Repro")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub),mcmc(m2.sub),mcmc(m3.sub)))

pdf("results/cov_Pexplained_surv_repro.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="% species-specific covariance",xlim = c(0,1))

dev.off()

