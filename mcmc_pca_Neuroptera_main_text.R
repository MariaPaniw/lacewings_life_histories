# Supporting Material S2
# Analyses for Comparative life-history responses of lacewings to changes in temperature
# created by Maria Paniw
# last modified: 04-11-2023

#############################   Load necessary packages
library(MCMCglmm)
library(tidyr)
library(ggplot2)
library(viridis)
library(hdrcde)

###########################   Data perparation 

df=read.csv("LH_Neuroptera_red.csv")%>%
  mutate(temp = as.numeric(temp),Dev_1st_inst = as.numeric(Dev_1st_inst),
         Dev_2nd_inst = as.numeric(Dev_2nd_inst),
         Dev_3rd_inst = as.numeric(Dev_3rd_inst),
         Dev_.P = as.numeric(Dev_.P),
         L_surv = as.numeric(L_surv),
         M_sex_rep_rate_M = as.numeric(M_sex_rep_rate_M))

# subset to the columns of interest
sub_noNA=na.omit(df[,c(2,4,10,35:38,58,39)])

# Values necessary to back-transform scaled temperatures for plotting
mean.temp=25.1569
sd.temp=3.571334

sub_noNA$temp=as.numeric(scale(sub_noNA$temp))

sub_noNA=sub_noNA[sub_noNA$Dev_1st_inst>0,]
sub_noNA=sub_noNA[sub_noNA$Dev_.P>0,]
sub_noNA=sub_noNA[sub_noNA$M_sex_rep_rate_M>0,]

sub_noNA$L_surv=asin(sqrt(sub_noNA$L_surv/100)) # arcsin transformation of survival to normalize data

nrow(sub_noNA) # final sample size

sub_noNA[,4:8]=log(sub_noNA[,4:8])

sub_noNA$temp2=sub_noNA$temp^2

###########################   MCMC analyses

prior = list(R = list(V = diag(6)/7, n = 6, nu=0.002),
             G = list(G1 = list(V = diag(6)/7, n = 6, nu=0.002),
                      G2 = list(V = diag(2)/7, n = 6, nu=0.002)))

m1=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P,M_sex_rep_rate_M,L_surv)~trait+trait:temp+trait:temp2,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 6), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)

m2=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P,M_sex_rep_rate_M,L_surv)~trait+trait:temp+trait:temp2,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 6), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)

m3=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P,M_sex_rep_rate_M,L_surv)~trait+trait:temp+trait:temp2,
            random = ~ us(trait):sp. + us(1 + temp):St_ID,rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 6), nitt = 60000, burnin = 10000,
            pr=F,thin=25, data = sub_noNA)

library(coda)

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
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,7])*temp.pred+mean(out.mcmc[,13])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,7],0.025)*temp.pred+quantile(out.mcmc[,13],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,7],0.975)*temp.pred+quantile(out.mcmc[,13],0.975)*temp.pred2)),
                      stage="Dev_1st_inst")
new.data_2=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,2])+mean(out.mcmc[,8])*temp.pred+mean(out.mcmc[,14])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,2],0.025)+quantile(out.mcmc[,8],0.025)*temp.pred+quantile(out.mcmc[,14],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,2],0.975)+quantile(out.mcmc[,8],0.975)*temp.pred+quantile(out.mcmc[,14],0.975)*temp.pred2)),
                      stage="Dev_2nd_inst")

new.data_3=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,3])+mean(out.mcmc[,9])*temp.pred+mean(out.mcmc[,15])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,3],0.025)+quantile(out.mcmc[,9],0.025)*temp.pred+quantile(out.mcmc[,15],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,3],0.975)+quantile(out.mcmc[,9],0.975)*temp.pred+quantile(out.mcmc[,15],0.975)*temp.pred2)),
                      stage="Dev_3rd_inst")

new.data_4=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,4])+mean(out.mcmc[,10])*temp.pred+mean(out.mcmc[,16])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,4],0.025)+quantile(out.mcmc[,10],0.025)*temp.pred+quantile(out.mcmc[,16],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,4],0.975)+quantile(out.mcmc[,10],0.975)*temp.pred+quantile(out.mcmc[,16],0.975)*temp.pred2)),
                      stage="Dev_.P")

new.data_5=data.frame(temp=temp.pred,
                      dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,5])+mean(out.mcmc[,11])*temp.pred+mean(out.mcmc[,17])*temp.pred2),
                      LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,5],0.025)+quantile(out.mcmc[,11],0.025)*temp.pred+quantile(out.mcmc[,17],0.025)*temp.pred2)),
                      UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,5],0.975)+quantile(out.mcmc[,11],0.975)*temp.pred+quantile(out.mcmc[,17],0.975)*temp.pred2)),
                      stage="M_sex_rep_rate_M")

new.data_6=data.frame(temp=temp.pred,
                      dev=sin(mean(out.mcmc[,1])+mean(out.mcmc[,6])+mean(out.mcmc[,12])*temp.pred+mean(out.mcmc[,18])*temp.pred2)^2,
                      LB=sin((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,6],0.025)+quantile(out.mcmc[,12],0.025)*temp.pred+quantile(out.mcmc[,18],0.025)*temp.pred2))^2,
                      UB=sin((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,6],0.975)+quantile(out.mcmc[,12],0.975)*temp.pred+quantile(out.mcmc[,18],0.975)*temp.pred2))^2,
                      stage="L_surv")

pred.data=rbind(new.data_1,new.data_2,new.data_3,new.data_4,new.data_5,new.data_6)

pred.data$stage=factor(pred.data$stage,levels=c("Dev_1st_inst","Dev_2nd_inst","Dev_3rd_inst","Dev_.P","M_sex_rep_rate_M","L_surv"))

levels(pred.data$stage) = c("D 1st instar","D 2nd instar", "D 3rd instar","D Pupae","# eggs/female", "S pupae-adult")

#observed 

df.obs <- sub_noNA %>%
  pivot_longer(
    cols = "Dev_1st_inst":"L_surv",
    names_to = "stage",
    values_to = "dev"
  )

df.obs$dev[!df.obs$stage%in%"L_surv"]=exp(df.obs$dev[!df.obs$stage%in%"L_surv"])
df.obs$dev[df.obs$stage%in%"L_surv"]=sin(df.obs$dev[df.obs$stage%in%"L_surv"])^2


df.obs$stage=factor(df.obs$stage,levels=c("Dev_1st_inst","Dev_2nd_inst","Dev_3rd_inst","Dev_.P","M_sex_rep_rate_M","L_surv"))

levels(df.obs$stage)=c("D 1st instar","D 2nd instar", "D 3rd instar","D Pupae","# eggs/female","S pupae-adult")

df.obs$sp.=gsub("_"," ", df.obs$sp.)

#Back-transform temperature

pred.data$temp=pred.data$temp*sd.temp+mean.temp
df.obs$temp=df.obs$temp*sd.temp+mean.temp

# Figure 2 main text
p.temp=ggplot(pred.data, aes(temp, dev))+
  facet_grid(stage~.,scales="free")+
  geom_line() +
  geom_point(data=df.obs,aes(temp, dev, col=sp.),size=2)+
  geom_ribbon(aes(ymin = LB, ymax = UB),alpha=0.1,col=NA) +
  xlab("Temperature (ºC)")+
  ylab("")+
  theme_bw(base_size = 18)+
  theme(legend.position = "none",
        legend.text = element_text(face = "italic",size=10))+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=16))

p.temp


#covariance
colnames(m1$VCV)=gsub(".units","",colnames(m1$VCV))
colnames(m1$VCV)=gsub("trait","",colnames(m1$VCV))

colnames(m2$VCV)=gsub(".units","",colnames(m2$VCV))
colnames(m2$VCV)=gsub("trait","",colnames(m2$VCV))

colnames(m3$VCV)=gsub(".units","",colnames(m3$VCV))
colnames(m3$VCV)=gsub("trait","",colnames(m3$VCV))

## Residual variance
m1.sub.res=m1$VCV[,c(37:42,44:48,51:54,58:60,65,66,72)]
m2.sub.res=m2$VCV[,c(37:42,44:48,51:54,58:60,65,66,72)]
m3.sub.res=m3$VCV[,c(37:42,44:48,51:54,58:60,65,66,72)]

colnames(m1.sub.res)=colnames(m2.sub.res)=colnames(m3.sub.res)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                   "1st instar - pupae","1st instar - repro","1st instar - surv","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae","2nd instar - repro","2nd instar - surv",
                   "3d - 3rd instar","3rd instar - pupae","3rd instar - repro","3rd instar - surv","pupae - pupae","pupae - repro","pupae - surv","repro - repro","repro - surv","surv - surv")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.res),mcmc(m2.sub.res),mcmc(m3.sub.res)))

library(MCMCvis)

MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_surv+",wd="~/plots/")

pdf("~/Fig3b_main_text.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")

dev.off()


#### Species level covariance

m1.sub.sp=m1$VCV[,c(1:6,8:12,15:18,22:24,29,30,36)]
m2.sub.sp=m2$VCV[,c(1:6,8:12,15:18,22:24,29,30,36)]
m3.sub.sp=m3$VCV[,c(1:6,8:12,15:18,22:24,29,30,36)]

colnames(m1.sub.sp)=colnames(m2.sub.sp)=colnames(m3.sub.sp)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                                                                 "1st instar - pupae","1st instar - repro","1st instar - surv","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae","2nd instar - repro","2nd instar - surv",
                                                                 "3d - 3rd instar","3rd instar - pupae","3rd instar - repro","3rd instar - surv","pupae - pupae","pupae - repro","pupae - surv","repro - repro","repro - surv","surv - surv")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.sp),mcmc(m2.sub.sp),mcmc(m3.sub.sp)))

library(MCMCvis)


MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_surv+_sp",wd="~/plots/")

pdf("~/Fig3a_main_text.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random species-specific covariance")

dev.off()

###### Percent variance explained by species

m1.sub=abs(m1.sub.sp)/(abs(m1.sub.sp)+abs(m1.sub.res))
m2.sub=abs(m2.sub.sp)/(abs(m2.sub.sp)+abs(m2.sub.res))
m3.sub=abs(m3.sub.sp)/(abs(m3.sub.sp)+abs(m3.sub.res))

colnames(m1.sub)=colnames(m2.sub)=colnames(m3.sub)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                                                     "1st instar - pupae","1st instar - repro","1st instar - surv","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae","2nd instar - repro","2nd instar - surv",
                                                     "3d - 3rd instar","3rd instar - pupae","3rd instar - repro","3rd instar - surv","pupae - pupae","pupae - repro","pupae - surv","repro - repro","repro - surv","surv - surv")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub),mcmc(m2.sub),mcmc(m3.sub)))

pdf("~/cov_Pexplained.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="% species-specific covariance",xlim = c(0,1))

dev.off()


######### PCA ANALYSES

### Get residuals

temp.pred=sub_noNA$temp
temp.pred2=sub_noNA$temp2

new.data_1=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,7])$mode*temp.pred+hdr(out.mcmc[,13])$mode*temp.pred2,
                      stage="Dev_1st_inst")
new.data_2=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,2])$mode+hdr(out.mcmc[,8])$mode*temp.pred+hdr(out.mcmc[,14])$mode*temp.pred2,
                      stage="Dev_2nd_inst")

new.data_3=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,3])$mode+hdr(out.mcmc[,9])$mode*temp.pred+hdr(out.mcmc[,15])$mode*temp.pred2,
                      stage="Dev_3rd_inst")

new.data_4=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,4])$mode+hdr(out.mcmc[,10])$mode*temp.pred+hdr(out.mcmc[,16])$mode*temp.pred2,
                      stage="Dev_.P")

new.data_5=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,5])$mode+hdr(out.mcmc[,11])$mode*temp.pred+hdr(out.mcmc[,17])$mode*temp.pred2,
                      stage="M_sex_rep_rate_M")

new.data_6=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,6])$mode+hdr(out.mcmc[,12])$mode*temp.pred+hdr(out.mcmc[,18])$mode*temp.pred2,
                      stage="L_surv")

#Main PCA dataset
predicted_data <- data.frame(
  Dev_1st_inst = sub_noNA$Dev_1st_inst-new.data_1$dev,
  Dev_2nd_inst = sub_noNA$Dev_2nd_inst-new.data_2$dev,
  Dev_3rd_inst = sub_noNA$Dev_3rd_inst-new.data_3$dev,
  Dev_.P = sub_noNA$Dev_.P-new.data_4$dev,
  Rep = sub_noNA$M_sex_rep_rate_M-new.data_5$dev,
  Surv = sub_noNA$L_surv-new.data_6$dev
)
#### PCA 
pca_res <- prcomp(predicted_data,center=T,scale. = T)
pca_res

results <- pca_res
summary(results)$importance # information on the importance of PCA axes (how much variance in data they describe)

ncomp=2 # based on Kaiser criterion

# Varimax transforamtion

rawLoadings     <- results$rotation[,1:ncomp] %*% diag(results$sdev, ncomp, ncomp)
rotatedLoadings <- varimax(rawLoadings)$loadings
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- scale(predicted_data) %*% invLoadings


# This creates a new objec, x, that we will use to make the plot
x <- list() 
x$scores <- scores
colnames(x$scores)=c("PC1", "PC2")
x$loadings <- rotatedLoadings[,1:2]
colnames(x$loadings)=c("PC1", "PC2")

name1=0.539
name2=0.225

# PLOT

library(ggrepel)

PCA_S=data.frame(PC1=x$scores[,1],PC2=x$scores[,2])

LHTexpr <- c("D 1st instar","D 2nd instar", "D 3rd instar","D Pupae","# eggs/female","S pupae-adult")

arrow.dat=as.data.frame(x$scores)

PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  
  plot <- ggplot(PCA_S, aes_string(x=x, y=y)) 
  plot <- plot + geom_point(color="grey",alpha=0.5,size=3)
  
  datapc <- data.frame(varnames=LHTexpr, PC$loadings)
  mult <- min(
    (max(PCA_S[,y]) - min(PCA_S[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(PCA_S[,x]) - min(PCA_S[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .3 * mult* (get(x)),
                      v2 = .3 * mult* (get(y))
  )
  
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.5,"cm")), size=1.2, color="black") 
  plot<-plot + geom_text_repel(data=datapc, aes(x=v1, y=v2, label=LHTexpr), size = 5,  parse=F,color="black", point.padding = unit(1, 'lines'))
  
  plot <- plot+theme_bw()+theme(panel.grid = element_blank())+ylab(paste("PC 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PC 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=18,colour = "black"))+theme(axis.title = element_text(size=20))+theme(legend.text = element_text(size=20,colour = "black"),
                                                                                                                        legend.title=element_blank(),
                                                                                                                        legend.key.size = unit(1, "lines"),
                                                                                                                        legend.position = "none")
  
  plot
}


# Fgiure 4 main text
comparison = PCbiplot(x)
comparison


