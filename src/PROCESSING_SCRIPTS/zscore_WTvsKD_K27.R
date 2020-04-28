# Author: Heurteau Alexandre
# Date: Mon Nov 20 10:12:13 2017

#####################################################################################--
#
#          DESCRIPTION
#
#####################################################################################--

#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          PATH  ----
#####################################################################################-
SPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
OUT.d <- create(paste0(SPATH,"/"))

#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source(paste0(OUT.d,"/LIBS/lrc_lib.R"))



#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

WT1 <- read.table(paste0(OUT.d,'DATA/PROCESSED/WT4KBgene.prof'),row.names=1)
KD1 <-read.table(paste0(OUT.d,'DATA/PROCESSED/KD4KBgene.prof'),row.names=1)

#####################################################################################-
#
#          RUN
#
#####################################################################################-

# Normalize reads number
KD1_norm=KD1*(sum(WT1)/sum(KD1))
plot(seq(-2000,2000,l=400),apply(WT1,2,mean),type='l',xlab='distance from TSS',ylab='average K27 level')
lines(seq(-2000,2000,l=400),apply(KD1_norm,2,mean),col='red')
abline(v=c(-500,0),lty=2)
# Seems ok

# Gènes avec que des 0 sur [-500,0]
sum(apply(WT1[,151:200],1,max)==0)
#~ [1] 812
sum(apply(KD1_norm[,151:200],1,max)==0)
#~ [1] 607
sum(apply(KD1_norm[,151:200],1,max)==0 & apply(WT1[,151:200],1,max)==0)
#~ [1] 420


# Remove gènes sans valeur de reads dans l'intervalle -500;0 #14165
WT1_no0=WT1[-which(apply(KD1_norm[,151:200],1,max)==0 & apply(WT1[,151:200],1,max)==0),]
KD1_no0=KD1_norm[-which(apply(KD1_norm[,151:200],1,max)==0 & apply(WT1[,151:200],1,max)==0),]


# COMPUTE ZSCORES ----
WT=WT1_no0
KD=KD1_no0

#calcul "zscore" sur (-500;0)
KD_sum=apply(KD[151:200],1,sum)
WT_sum=apply(WT[151:200],1,sum)
zsm500bp=(KD_sum-WT_sum)/sqrt((KD_sum+WT_sum)/2)

#calcul "zscore" sur +- 1Kbp
KD_sum=apply(KD[101:300],1,sum)
WT_sum=apply(WT[101:300],1,sum)
zspm1kb=(KD_sum-WT_sum)/sqrt((KD_sum+WT_sum)/2)

#calcul "zscore" sur +- 2Kbp
KD_sum=apply(KD[1:400],1,sum)
WT_sum=apply(WT[1:400],1,sum)
zspm2kb=(KD_sum-WT_sum)/sqrt((KD_sum+WT_sum)/2)



# summary(zs)
# hist(zs)
# length(zs[zs>0])
# #~ [1] 7332
# shapiro.test(sample(zs,4000))
# #Distribution des z-scores
# plot(density(zs),type="n",main="z-scores density",ylim=c(0,0.4))
# lines(density(zs),col="blue")
# lines(density(rnorm(10000)),col="red")
# abline(v=c(0,mean(zs)),col=c('red','blue'),lty=2)
# Ok follow normal distribution
# ----



# Adjust PVAL
pvm500 <- 2*pnorm(abs(zsm500bp),lower.tail = F)
pv_adj_m500 <- p.adjust(pvm500,method="BH")

pvpm1kb <- 2*pnorm(abs(zspm1kb),lower.tail = F)
pv_adj_pm1kb <- p.adjust(pvpm1kb,method="BH")

pvpm2kb <- 2*pnorm(abs(zspm2kb),lower.tail = F)
pv_adj_pm2kb <- p.adjust(pvpm2kb,method="BH")


mydf <- data.frame(genes=rownames(WT1_no0),zscore=zsm500bp,pv=pvm500,pv_adj=pv_adj_m500)
saveRDS(mydf,paste0(OUT.d,"/DATA/PROCESSED/ZSCORE/K27_wt_kd_zscoreM500_matrecap.RData"))

mydf <- data.frame(genes=rownames(WT1_no0),zscore=zspm1kb,pv=pvpm1kb,pv_adj=pv_adj_pm1kb)
saveRDS(mydf,paste0(OUT.d,"/DATA/PROCESSED/ZSCORE//K27_wt_kd_zscorePM1KB_matrecap.RData"))

mydf <- data.frame(genes=rownames(WT1_no0),zscore=zspm2kb,pv=pvpm2kb,pv_adj=pv_adj_pm2kb)
saveRDS(mydf,paste0(OUT.d,"/DATA/PROCESSED/ZSCORE/K27_wt_kd_zscorePM2KB_matrecap.RData"))
