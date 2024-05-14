# Supporting Material S2
# Analyses for Comparative life-history responses of lacewings to changes in temperature
# created by Maria Paniw, adapted by Sanne Evers
# last modified: 08/05/2024

set.seed(14052024)

library(MCMCglmm)
library(tidyr)
library(dplyr)
library(coda)
library(ggplot2)

df=read.ggplot2df=read.csv("data/LH_Neuroptera_red.csv") %>%
  mutate(temp = as.numeric(temp),
         Dev_1st_inst = as.numeric(Dev_1st_inst),
         Dev_2nd_inst = as.numeric(Dev_2nd_inst),
         Dev_3rd_inst = as.numeric(Dev_3rd_inst),
         Dev_.P = as.numeric(Dev_.P)) 


mean.temp=24
sd.temp=3.84


sub_noNA=na.omit(df[,c('St_ID', 'sp.', 'vivo_situ', "Types_of_Lit_Sources", 'temp', 'Latitude', 'Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')])

sub_noNA$temp=as.numeric(scale(sub_noNA$temp))
sub_noNA$Latitude=as.numeric(scale(sub_noNA$Latitude))

sub_noNA=sub_noNA[sub_noNA$vivo_situ%in%"in_situ",-3]

sub_noNA=sub_noNA[sub_noNA$Dev_1st_inst>0,]
sub_noNA=sub_noNA[sub_noNA$Dev_.P>0,]
sub_noNA[,c('Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')]=log(sub_noNA[,c('Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')])

prior = list(R = list(V = diag(4)/5, n = 4, nu=0.002),
             G = list(G1 = list(V = diag(4)/5, n = 4, nu=0.002),
                      G2 = list(V = diag(2)/5, n = 4, nu=0.002)))

m1=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait:temp + trait:Latitude,
            random = ~ us(trait):sp. + us(1 + temp):St_ID, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 100000, burnin = 50000,
            pr=F,thin=25, data = sub_noNA)

m2=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait:temp + trait:Latitude,
            random = ~ us(trait):sp. + us(1 + temp):St_ID, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 100000, burnin = 50000,
            pr=F,thin=25, data = sub_noNA)

m3=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait:temp + trait:Latitude,
            random = ~ us(trait):sp. + us(1 + temp):St_ID, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 100000, burnin = 50000,
            pr=F,thin=25, data = sub_noNA)


# Plot main (fixed) effects

param.coda.fixed=mcmc.list(list(mcmc(m1$Sol),mcmc(m2$Sol),mcmc(m3$Sol)))

summary(param.coda.fixed)
gelman.diag(param.coda.fixed,multivariate=F)

#Trace plots (to check if chains are well mixed)

par(mar=c(2,2,2,2))
plot(param.coda.fixed,smooth=F) # The different colors indicate different chains

### Plot response as function of temperature
out.mcmc=rbind(m1$Sol,m2$Sol,m3$Sol)

temp.pred=seq(min(sub_noNA$temp),max(sub_noNA$temp),length.out=20)

new.data_1=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,2])*temp.pred),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,2],0.025)*temp.pred)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,2],0.975)*temp.pred)),
                      stage="Dev_1st_inst")
new.data_2=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,3])*temp.pred),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,3],0.025)*temp.pred)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,3],0.975)*temp.pred)),
                      stage="Dev_2nd_inst")

new.data_3=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,4])*temp.pred),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,4],0.025)*temp.pred)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,4],0.975)*temp.pred)),
                      stage="Dev_3rd_inst")

new.data_4=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,5])*temp.pred),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,5],0.025)*temp.pred)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,5],0.975)*temp.pred)),
                      stage="Dev_.P")

pred.data=rbind(new.data_1,new.data_2,new.data_3,new.data_4)

pred.data$stage=factor(pred.data$stage,levels=c("Dev_1st_inst","Dev_2nd_inst","Dev_3rd_inst","Dev_.P"))

levels(pred.data$stage) = c("1st instar","2nd instar", "3rd instar","Pupae")

#observed 

df.obs <- sub_noNA %>%
  pivot_longer(
    cols = "Dev_1st_inst":"Dev_.P",
    names_to = "stage",
    values_to = "dev"
  )

df.obs$dev=exp(df.obs$dev)
df.obs$stage=factor(df.obs$stage,levels=c("Dev_1st_inst","Dev_2nd_inst","Dev_3rd_inst","Dev_.P"))

levels(df.obs$stage)=c("1st instar","2nd instar", "3rd instar","Pupae")
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

ggsave(filename = "results/plot_dev_in_vivo.pdf",plot=p.temp,width = 12,height = 8)
ggsave(filename = "results/plot_dev_in_vivo.png",plot=p.temp,width = 12,height = 8, dpi = 600)


