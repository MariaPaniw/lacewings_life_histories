# Supporting Material S2
# Analyses for Comparative life-history responses of lacewings to changes in temperature
# created by Maria Paniw
# last modified: 04-11-2023

set.seed(14052024)

#############################   Load necessary packages
library(MCMCglmm)
library(tidyr)
library(ggplot2)
library(viridis)
library(hdrcde)

###########################   Data perparation 

df=read.csv("LH_Neuroptera_red.csv")

# Values necessary to back-transform scaled temperatures for plotting
mean.temp=24
sd.temp=3.84

sub_noNA=na.omit(df[,c('St_ID', 'sp.', 'vivo_situ', "Types_of_Lit_Sources", 'temp', 'Latitude', 'Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')])

sub_noNA$temp=as.numeric(scale(sub_noNA$temp))

sub_noNA=sub_noNA[sub_noNA$Dev_1st_inst>0,]
sub_noNA=sub_noNA[sub_noNA$Dev_.P>0,]
sub_noNA[,c('Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')]=log(sub_noNA[,c('Dev_1st_inst', 'Dev_2nd_inst', 'Dev_3rd_inst','Dev_.P')])

sub_noNA$species_pub=as.factor(paste(sub_noNA$sp.,sub_noNA$St_ID))

###########################   MCMC analyses

prior = list(R = list(V = diag(4)/5, n = 4, nu=0.002),
             G = list(G1 = list(V = diag(4)/5, n = 4, nu=0.002),
                      G2 = list(V = diag(2)/5, n = 4, nu=0.002)))

m1=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait+trait:temp + trait:Latitude,
            random = ~ us(trait):species_pub + us(1 + temp):species_pub, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 60000, burnin = 10000,
            pr=T,thin=25, data = sub_noNA)

m2=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait+trait:temp + trait:Latitude,
            random = ~ us(trait):species_pub + us(1 + temp):species_pub, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 60000, burnin = 10000,
            pr=T,thin=25, data = sub_noNA)

m3=MCMCglmm(cbind(Dev_1st_inst,Dev_2nd_inst,Dev_3rd_inst,Dev_.P)~trait+trait:temp + trait:Latitude,
            random =  ~ us(trait):species_pub + us(1 + temp):species_pub, rcov = ~us(trait):units,prior = prior, family = rep("gaussian", 4), nitt = 60000, burnin = 10000,
            pr=T,thin=25, data = sub_noNA)

library(coda)

param.coda.fixed=mcmc.list(list(mcmc(m1$Sol[,1:12]),mcmc(m2$Sol[,1:12]),mcmc(m3$Sol[,1:12])))

summary(param.coda.fixed)
gelman.diag(param.coda.fixed,multivariate=F)

#Trace plots (to check if chains are well mixed)

par(mar=c(2,2,2,2))
plot(param.coda.fixed,smooth=F) # The different colors indicate different chains

### Plot response as function of temperature
out.mcmc=rbind(m1$Sol,m2$Sol,m3$Sol)

temp.pred=seq(min(sub_noNA$temp),max(sub_noNA$temp),length.out=20)

new.data_1=expand.grid(temp=temp.pred,
                       species = unique(sub_noNA$species_pub)) %>%
  rowwise() %>%
  mutate(dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,5])*temp + 
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][13:37], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)]) * temp),
         LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,5],0.025)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][13:37], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.025) * temp)),
         UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,5],0.975)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][13:37], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.975) * temp)),
         stage="Dev_1st_inst")
  
new.data_2=expand.grid(temp=temp.pred,
                       species = unique(sub_noNA$species_pub)) %>%
  rowwise() %>%
  mutate(dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,2])+mean(out.mcmc[,6])*temp + 
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][38:62], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)]) * temp),
         LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,2],0.025)+quantile(out.mcmc[,6],0.025)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][38:62], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.025) * temp)),
         UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,2],0.025)+quantile(out.mcmc[,6],0.975)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][38:62], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.975) * temp)),
         stage="Dev_2nd_inst")

new.data_3=expand.grid(temp=temp.pred,
                       species = unique(sub_noNA$species_pub)) %>%
  rowwise() %>%
  mutate(dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,3])+mean(out.mcmc[,7])*temp + 
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][63:87], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)]) * temp),
         LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,3],0.025)+quantile(out.mcmc[,7],0.025)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][63:87], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.025) * temp)),
         UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,3],0.025)+quantile(out.mcmc[,7],0.975)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][63:87], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.975) * temp)),
         stage="Dev_3rd_inst")


new.data_4=expand.grid(temp=temp.pred,
                       species = unique(sub_noNA$species_pub)) %>%
  rowwise() %>%
  mutate(dev=exp(mean(out.mcmc[,1])+mean(out.mcmc[,4])+mean(out.mcmc[,8])*temp + 
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][88:112], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)]) +
                   mean(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)]) * temp),
         LB=exp((quantile(out.mcmc[,1],0.025)+quantile(out.mcmc[,4],0.025)+quantile(out.mcmc[,8],0.025)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][88:112], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.025) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.025) * temp)),
         UB=exp((quantile(out.mcmc[,1],0.975)+quantile(out.mcmc[,4],0.025)+quantile(out.mcmc[,8],0.975)*temp+
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][88:112], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][113:137], value = T)],0.975) +
                   quantile(out.mcmc[, grep(pattern = species, attr(out.mcmc, which = "dimnames")[[2]][138:162], value = T)],0.975) * temp)),
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

df.obs$species=as.factor(paste(df.obs$sp.,df.obs$St_ID))

levels(df.obs$stage)=c("1st instar","2nd instar", "3rd instar","Pupae")
df.obs$sp.=gsub("_"," ", df.obs$sp.)


pred.data$temp=pred.data$temp*sd.temp+mean.temp
df.obs$temp=df.obs$temp*sd.temp+mean.temp


# get only observed temp ranges

range_df <- df.obs %>%
  group_by(species) %>%
  summarise(min_temp = min(temp) * 0.9,
            max_temp = max(temp) * 1.1)

pred.data1 <- pred.data %>% 
  left_join(., range_df, by = "species") %>%
  rowwise() %>%
  filter(temp >= min_temp & temp <= max_temp)




# Figure S3.3

p.temp=ggplot(pred.data1, aes(temp, dev))+
  facet_grid(stage~.,scales="free")+
  geom_point(data=df.obs,aes(temp, dev, col=species),size=2, show.legend = FALSE)+
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = species),alpha=0.2) +
  geom_line(aes(colour = species), show.legend = FALSE) +
  xlab("Temperature (ºC)")+
  ylab("")+
  theme_bw(base_size = 18)+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5),
        legend.text = element_text(face = "italic",size=10))+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(strip.background =element_blank(),
        strip.text =element_text(size=16))+
  guides(fill = guide_legend(title = "Species & Study ID", override.aes = list(alpha = 0.8)))

p.temp


ggsave(filename = "results/plot_dev_FigS3.3.pdf",plot=p.temp,width = 13,height = 10)
ggsave(filename = "results/plot_dev_FigS3.3.png",plot=p.temp,width = 13,height = 10, dpi = 600)



#covariance
colnames(m1$VCV)=gsub(".units","",colnames(m1$VCV))
colnames(m1$VCV)=gsub("trait","",colnames(m1$VCV))

colnames(m2$VCV)=gsub(".units","",colnames(m2$VCV))
colnames(m2$VCV)=gsub("trait","",colnames(m2$VCV))

colnames(m3$VCV)=gsub(".units","",colnames(m3$VCV))
colnames(m3$VCV)=gsub("trait","",colnames(m3$VCV))

### Residual variance
m1.sub.res=m1$VCV[,c(17:20,22:24,27,28,32)]
m2.sub.res=m2$VCV[,c(17:20,22:24,27,28,32)]
m3.sub.res=m3$VCV[,c(17:20,22:24,27,28,32)]
colnames(m1.sub.res)=colnames(m2.sub.res)=colnames(m3.sub.res)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                   "1st instar - pupae","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae",
                   "3d - 3rd instar","3rd instar - pupae","pupae - pupae")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.res),mcmc(m2.sub.res),mcmc(m3.sub.res)))

library(MCMCvis)


MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_dev",wd="results")

pdf("results/FigS2.3b.pdf",width=6,height=7)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")
dev.off()

png("results/FigS2.3b.png", width = 6, height = 7, units = "in", res = 400)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Residual covariance")
dev.off()


#### Species level covariance

m1.sub.sp=m1$VCV[,c(1:4,6:8,11,12,16)]
m2.sub.sp=m2$VCV[,c(1:4,6:8,11,12,16)]
m3.sub.sp=m3$VCV[,c(1:4,6:8,11,12,16)]
colnames(m1.sub.sp)=colnames(m2.sub.sp)=colnames(m3.sub.sp)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                                                                 "1st instar - pupae","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae",
                                                                 "3d - 3rd instar","3rd instar - pupae","pupae - pupae")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub.sp),mcmc(m2.sub.sp),mcmc(m3.sub.sp)))

library(MCMCvis)

MCMCtrace(param.coda.vcv,pdf=T,filename="Fig.covariance_dev_sp",wd="results/")

pdf("results/FigS.2.3a.pdf",width=6,height=7)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random species-specific covariance")
dev.off()

png("results/FigS2.3a.png", width = 6, height = 7, units = "in", res = 400)
MCMCplot(param.coda.vcv,ref_ovl = T,xlab="Random species-specific covariance")
dev.off()


###### Percent variance explained by species

m1.sub=abs(m1.sub.sp)/(abs(m1.sub.sp)+abs(m1.sub.res))
m2.sub=abs(m2.sub.sp)/(abs(m2.sub.sp)+abs(m2.sub.res))
m3.sub=abs(m3.sub.sp)/(abs(m3.sub.sp)+abs(m3.sub.res))

colnames(m1.sub)=colnames(m2.sub)=colnames(m3.sub)=c("1st - 1st instar","1st - 2nd instar","1st - 3rd instar",
                                                     "1st instar - pupae","2nd - 2nd instar","2nd - 3rd instar","2nd instar - pupae",
                                                     "3d - 3rd instar","3rd instar - pupae","pupae - pupae")

param.coda.vcv=mcmc.list(list(mcmc(m1.sub),mcmc(m2.sub),mcmc(m3.sub)))

pdf("results/cov_Pexplained_dev.pdf",width=6,height=7)

MCMCplot(param.coda.vcv,ref_ovl = T,xlab="% species-specific covariance",xlim = c(0,1))

dev.off()


######### PCA ANALYSES

### Get residuals

temp.pred=sub_noNA$temp

new.data_1=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,5])$mode*temp.pred,
                      stage="Dev_1st_inst")
new.data_2=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,2])$mode+hdr(out.mcmc[,6])$mode*temp.pred,
                      stage="Dev_2nd_inst")

new.data_3=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,3])$mode+hdr(out.mcmc[,7])$mode*temp.pred,
                      stage="Dev_3rd_inst")

new.data_4=data.frame(temp=temp.pred,
                      dev=hdr(out.mcmc[,1])$mode+hdr(out.mcmc[,4])$mode+hdr(out.mcmc[,8])$mode*temp.pred,
                      stage="Dev_.P")
#Main PCA dataset
predicted_data <- data.frame(
  Dev_1st_inst = sub_noNA$Dev_1st_inst-new.data_1$dev,
  Dev_2nd_inst = sub_noNA$Dev_2nd_inst-new.data_2$dev,
  Dev_3rd_inst = sub_noNA$Dev_3rd_inst-new.data_3$dev,
  Dev_.P = sub_noNA$Dev_.P-new.data_4$dev
)

# This performs the PCA
pca_res <- prcomp(predicted_data,center=T,scale. = T)
pca_res

results <- pca_res
summary(results)$importance 

ncomp=2
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


name1=0.61
name2=0.27

# PLOT

library(ggrepel)

PCA_S=data.frame(PC1=x$scores[,1],PC2=x$scores[,2])

LHTexpr <- c("D 1st instar","D 2nd instar", "D 3rd instar","D Pupae")

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
  
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.5,"cm")), size=1.2, color="grey50") 
  
  # plot<-plot + geom_text_repel(data=datapc, aes(x=v1, y=v2, label=LHTexpr), size = 5,  parse=F,color="black", point.padding = unit(1, 'lines'))
  
  plot <- plot+theme_bw()+theme(panel.grid = element_blank())+ylab(paste("PC 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PC 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=18,colour = "black"))+theme(axis.title = element_text(size=20))+theme(legend.text = element_text(size=20,colour = "black"),
                                                                                                                        legend.title=element_blank(),
                                                                                                                        legend.key.size = unit(1, "lines"),
                                                                                                                        legend.position = "none")
  
  plot
}



comparison = PCbiplot(x)
comparison # Fig. S2.4


ggsave(filename = "results/PCA_dev_only.pdf",plot=comparison,width = 8,height = 10)
ggsave(filename = "results/PCA_dev_only.png",plot=comparison, width = 15.2, height = 15.2, units = "cm", dpi = 600)

