
#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          PATH  ----
#####################################################################################-
SPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
OUT.d <- paste0(SPATH,"/")
#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source(paste0(OUT.d,"/LIBS/lrc_lib.R"))
source(paste0(OUT.d,"/LIBS/data_lib.R")) # load data
source(paste0(OUT.d,"LIBS/mrc_lib.R"))
TMP.d <- create(paste0(SPATH,"/TMP/"))


#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-
# Figure 2A -- HISTOGRAM OF NORMR NO HMM 15 KB AWAY ------------------------------------

# HISTOGRAM OF NORMR NO HMM 15 KB AWAY
# Limit the size of the domains to c [100,5000]
# Distance to K27 dom at least 10 Kb
SUB <- "FIGURES/Figure_2/Figure_2A/"
bSup <- 5000
bInf <- 80
DTK27 <- 10e3

normRmicroTMP.gr <- subset(mDOM.dm3.gr,width > bInf & width <= bSup  & dtk27 > DTK27)
OUTF.d <- create(paste0(OUT.d,SUB))
pdf(paste0(OUTF.d,"HIST_mDOM_bsup",bSup,"_binf",bInf,"_DTk27_",DTK27,".pdf"),width=10)
hist(width(normRmicroTMP.gr),prob=TRUE,breaks=100,ylim=c(0,2e-3),xlim=c(0,5000),
     xlab = 'normR domains width',main='Histogram of normR domains width')
lines(density(width(normRmicroTMP.gr),bw=30), col="orange", lwd=1.1)
legend("topright", c("Domain size fits with nucleosome and linker", "Domain size is shifted from nucleosome and linker"), col=c("red", "grey"), lty = 16,lwd=1.5,cex = 0.7)
abline(v=c(147,(seq(1,20)*147)+(seq(1,20)-1)*60), col=c(rep("red",10),rep("grey",11)),lty=16,lwd=0.8,ylim=c(0.001,0.002))
dev.off()

#   -----------------------------------------------------------------------

# Figure 2B - Boxplot LFC of K27 level upon micro domains by size ---------
SUB <- "FIGURES/Figure_2/Figure_2B/"
bSup <- 5000
bInf <- 40
DTK27 <- 10e3
ylimT <- list(SQD=c(-1,1))

# MICRO DOMAINS
# Filter mDOM to be less thant 10kb in width
mDOM.sub.dm3.gr <- subset(mDOM.dm3.gr,width < 10e3)
microDomains.grl <- list()

normRmicroTMP_SQD.gr <- mDOM.sub.dm3.gr
for(my.bwp in ls(pattern ="h3k27(.*)norm.bwp")){
  NME <- ifelse(grepl("h3k27kd*",my.bwp),"k27KD","k27WT")
  normR_noHMM.cov <- importData(get(my.bwp),"BigWig",normRmicroTMP_SQD.gr)
  normR_noHMM.cov <- normR_noHMM.cov[[1]]
  normR_noHMM_mat.cov <- as.data.frame(as.matrix(normR_noHMM.cov[normRmicroTMP_SQD.gr]))
  mcols(normRmicroTMP_SQD.gr)[NME] <- rowSums(normR_noHMM_mat.cov,na.rm=T)/width(normRmicroTMP_SQD.gr)
}
normRmicroTMP_SQD.gr$zscore <- with(normRmicroTMP_SQD.gr,k27KD-k27WT/sqrt(abs(rowMeans(cbind(k27WT,k27KD)))))
normRmicroTMP_SQD.gr$lfc <- with(normRmicroTMP_SQD.gr,log2((k27KD + 1e-6) /(k27WT + 1e-6)))
microDomains.grl[["SQD"]] <- normRmicroTMP_SQD.gr


OUTF.d <- create(paste0(OUT.d,SUB))
for(type in names(microDomains.grl)){
    mydf <- NULL
    normR_nohmmFit.gr <- microDomains.grl[[type]]
    normR_nohmmFit.gr$lfc <- ifelse(is.nan(normR_nohmmFit.gr$lfc),NA,normR_nohmmFit.gr$lfc)
    supp15kb <- mcols(distanceToNearest(normR_nohmmFit.gr,Het.gr))$distance >= DTK27
    normR_nohmmFitsup15kb.gr <- normR_nohmmFit.gr[supp15kb]
    normR_nohmmFitsup15kb.gr$wdth <- width(normR_nohmmFitsup15kb.gr)
    nomrR_binf_bsup.gr <- subset(normR_nohmmFitsup15kb.gr,wdth < bSup & wdth > bInf) #& k27WT > k27KD)
    
    
    mydf <- data.frame(wdth=width(nomrR_binf_bsup.gr),lfc=nomrR_binf_bsup.gr$lfc)
    out_b <- create(paste0(OUTF.d,"/BOXPLOT/LFCK27",type,"_by_normR_DOMAINSSIZE/"))
    tiff(paste0(out_b,"BOXPLOT_LFC_DOMAINSIZE_normR_bsup",bSup,"_binf",bInf,"_DTk27_",DTK27,".tiff"), width = 8, height = 8, units = 'in', res = 300)
    p <- ggline(mydf,x="wdth",y="lfc",add="mean_se",plot_type = "p") +
      geom_vline(xintercept=147,col="red",linetype="dashed",alpha=0.6,size=1.4) + 
      theme_tiff
    print(ggpar(p, ylim=ylimT[[type]]))
    dev.off()
    
    pdf(paste0(out_b,"BOXPLOT_LFC_DOMAINSIZE_normR_bsup",bSup,"_binf",bInf,"_DTk27_",DTK27,".pdf"),width=30)
    p <- ggline(mydf,x="wdth",y="lfc",add="mean_se",plot_type = "p") +
      geom_vline(xintercept=147,col="red",linetype="dashed",alpha=0.6,size=1.4)  +
      geom_hline(yintercept = 0,col="grey",linetype="dashed",alpha=0.6,size=1.4)  +
      geom_vline(xintercept=c(147,(seq(1,25)*147)+(seq(1,25)-1)*60),col="grey",linetype="dashed",alpha=0.3)
    print(ggpar(p, ylim=ylimT[[type]]))
    dev.off()
}


#   -----------------------------------------------------------------------

# Figure 2C - Profile K27 WT/KD on mDOM ---------------------------------------------------------------
SUB <- "FIGURES/Figure_2/Figure_2C/"
TYPE <- "RPKM"
DIFF_SCORE <- "ZS"
NTILE <- 10
SEP_NTILE <- 2
DST_BORDER <- 5000

WDTH_MAX <- 1500
WDTH_MIN <- 40


ctl.bwp <- bw.ll[[TYPE]][["sig"]][1]
trt.bwp <- bw.ll[[TYPE]][["sig"]][2]

# microDOMAINS
mDOM.sub.dm3.gr <- subset(mDOM.dm3.gr,width < WDTH_MAX & width >= WDTH_MIN)
seqlevelsStyle(mDOM.sub.dm3.gr) <- "UCSC"

mDOM.sub.dm3.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = mDOM.sub.dm3.gr)
mDOM.sub.dm3.gr <- lfc_zs_dif_decile(anc.gr = mDOM.sub.dm3.gr,NTILE=NTILE,ASYM=T)

microWT_KD_EUC.gr <- subset(mDOM.sub.dm3.gr,dtk27 > DST_BORDER)
microWT_KD_HET.gr <- subset(mDOM.sub.dm3.gr,dteuc > DST_BORDER)

microWT_KD_DOWN_EUC.gr <- GRanges(data.frame(microWT_KD_EUC.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) > NTILE + SEP_NTILE)) 
microWT_KD_UP_EUC.gr <- GRanges(data.frame(microWT_KD_EUC.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) <= NTILE - SEP_NTILE)) 
microWT_KD_EQ_EUC.gr <- GRanges(data.frame(microWT_KD_EUC.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) <= NTILE + SEP_NTILE & UQ(as.name(paste0("dec",DIFF_SCORE))) > NTILE - SEP_NTILE)) 

microWT_KD_DOWN_HET.gr <- GRanges(data.frame(microWT_KD_HET.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) > NTILE + 7)) 
microWT_KD_UP_HET.gr <- GRanges(data.frame(microWT_KD_HET.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) < NTILE - 7)) 


microKD_WT_top50.gr <- microWT_KD_DOWN_EUC.gr
microKD_WT_min50.gr <- microWT_KD_EQ_EUC.gr
SD <- 2

OUT_PROF.d <- create(OUT.d,SUB)

pdf(paste0(OUT_PROF.d,"PROFILE_K27_WT_KD_TOP50perc_LOSS_K27_KD.pdf"))
bw.l <- c(h3k27kd_RPKM.bwp,h3k27wt_RPKM.bwp)
seqPlotSDoutliers(bw.l,tmp,c("microKD_WT_top50.gr"),ylim=c(0,35),xlim=c(1300,1300),bin=20L,sd=SD,err = F,smooth = T,type="pf",gnme = "dm3",ignore.strand = F,COLVEC = c("firebrick","dodgerblue"))
dev.off()

# WITH BOXPLOT
OUT_BOX.d <- create(OUT.d,SUB)
microKD_WT_top50.gr <- resize(microKD_WT_top50.gr,200,"start") # see profile
seqlevelsStyle(microKD_WT_top50.gr) <- "UCSC"
microKD_WT_top50.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microKD_WT_top50.gr)
all.df <- data.frame(microKD_WT_top50.gr)
melt.df <- melt(all.df,measure.vars = c("trtSIGNAL","ctlSIGNAL"))
melt.df$variable <- factor(melt.df$variable, levels = levels(melt.df$variable)[c(2,1)])
# ypos <- mean(Cmelt.df$value) - sd(melt.df$value)
ypos <- 5
pdf(paste0(OUT_BOX.d,"BOXPLOT_K27_WT_KD_SIG_AT_microDOMAINS_top50_DOWN_K127_NORM_RPKM_FROM_0BP_to_200BP.pdf"),width=5)
p <- ggplot(melt.df, aes(x = variable, y = value,color=variable),
            size=2,width=0.2) +
  scale_color_manual(values=c("dodgerblue", "firebrick"))+
  coord_cartesian(ylim = c(0,40)) +
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=0.2,width = 0.1) + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "linerange", size=0.2,linetype="dashed") + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=0.2,width = 0.1) +
  geom_boxplot(width=0.1,outlier.size = 0, coef = 0,size=0.2,outlier.shape=NA,position=position_dodge(0.6)) +
  stat_compare_means(aes(group = variable),paired=F,size=3,label.y = ypos,label.x = 1.5,
                     # method.args = list(alternative = "greater"),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     label = "p.signif") +
  theme_classic() + theme(plot.title = element_text(size = 10, face = "bold")) +
  ggtitle(paste0("Boxplot - Wilcoxon TEST \n WT vs KD BGbord or noBGbord"))
print(p)
dev.off()


#   -----------------------------------------------------------------------

# Figure 2D - Distribution plot K27 level in EUC ctl vs HET ctl vs mDOM---------------------------------------------------------------

# PARAMETERS
SUB <- "FIGURES/Figure_2/Figure_2D/"
WDTH_MAX <- 1500
WDTH_MIN <- 80
TYPE <- "SQDNORM"
DIFF_SCORE <- "ZS"
NTILE <- 10
SEP_NTILE <- 5


# MICRO DOMAINS
mDOM.sub.dm3.gr <- subset(mDOM.dm3.gr,width < WDTH_MAX & width >= WDTH_MIN)


# CONTROL microDOM
Euc_2KBK27.gr <- subset(Euc.gr,width > 5000)-2000
nmicroDOM <- length(subset(mDOM.sub.dm3.gr,dtk27 > 2000))
microEUC.gr <- subset(mDOM.sub.dm3.gr,dtk27 > 2000)
ctrl_microDOM.gr <- sapply(1:length(microEUC.gr),
                           function(X){
                             wdth <- width(microEUC.gr[X])
                             tileEUC.gr <- sample(Euc_2KBK27.gr,1)
                             chr <- seqnames(tileEUC.gr)
                             start <- sample(start(tileEUC.gr):end(tileEUC.gr),1)
                             end <- start + wdth
                             GRanges(seqnames=chr,IRanges(start=start,end=end))
                           }
)
ctrlEUC_microDOM.gr <- do.call("c",ctrl_microDOM.gr)

# CONTROL HET mDOM
Het_2KBK27.gr <- subset(Het.gr,width > 5000)-2000
nmicroDOM <- length(subset(mDOM.sub.dm3.gr,dtk27 > 2000))
microEUC.gr <- subset(mDOM.sub.dm3.gr,dtk27 > 2000)
ctrl_microDOM.gr <- sapply(1:length(microEUC.gr),
                           function(X){
                             wdth <- width(microEUC.gr[X])
                             tileHET.gr <- sample(Het_2KBK27.gr,1)
                             chr <- seqnames(tileHET.gr)
                             start <- sample(start(tileHET.gr):end(tileHET.gr),1)
                             end <- start + wdth
                             GRanges(seqnames=chr,IRanges(start=start,end=end))
                           }
)
ctrlHET_microDOM.gr <- do.call("c",ctrl_microDOM.gr)


# COMPUTE DISTRIBUTION OF K27
ctl.bwp <- bw.ll[[TYPE]][["sig"]][1]
trt.bwp <- bw.ll[[TYPE]][["sig"]][2]

seqlevelsStyle(microEUC.gr) <- "Ensembl"
ctrlHET_microDOM.gr <- bw_signal_simple(ctl.bwp = ctl.bwp,anc.gr = ctrlHET_microDOM.gr,name="K27_WT")
ctrlEUC_microDOM.gr <- bw_signal_simple(ctl.bwp = ctl.bwp,anc.gr = ctrlEUC_microDOM.gr,name="K27_WT")
microEUC.gr <- bw_signal_simple(ctl.bwp = ctl.bwp,anc.gr = microEUC.gr,name="K27_WT")

myDF <- rbind.data.frame(
  cbind.data.frame(NM="CTL_HET",VAL=ctrlHET_microDOM.gr$K27_WT),
  cbind.data.frame(NM="CTL_EUC",VAL=ctrlEUC_microDOM.gr$K27_WT),
  cbind.data.frame(NM="mDOM",VAL=microEUC.gr$K27_WT)
)

OUTF.d <- create(OUT.d,SUB)

pdf(paste0(OUTF.d,"/distribution_plot_mDOM_WDTHMAX_",WDTH_MAX,"_WDTHMIN_",WDTH_MIN,".pdf"))
gghistogram(myDF,x="VAL",y="..density..",add = "mean", col="NM", fill="NM",alpha = 0.6,palette="jco",bins=500) + coord_cartesian(xlim = c(0,0.05))
ggdensity(myDF,x="VAL",add = "mean", col="NM", fill="NM",palette="jco",bins=100) + scale_x_continuous(limits = c(0,0.05))
ggdensity(myDF,x="VAL",add = "mean", col="NM", fill="NM",palette="jco",bins=100) + scale_x_continuous(limits = c(0,0.05))+ theme_tiff
dev.off()

#   -----------------------------------------------------------------------

# Figure 2E - Distrib GAF /CTCF down K27 ---------------------------------------------------------------
SUB <-  "FIGURES/Figure_2/Figure_2E/"
OUTFD.d <- create(OUT.d,SUB)

gafnob <- gaf_kc  & !bg & cp & bs

myL <- NULL
mybait <- c("gafnob")
myfilt <- c("dwnK27bM500")
myanchor <- c("bgbord","nobgbord")
for(anchor in myanchor){
  for(bait in mybait){
    for(filt in myfilt){
      myL[[bait]][[anchor]][[filt]] <- fisherDistPlotBothSided(anchor,bait,filt,60000,2000,10000)
    }
  }
}



bg_gafnobg.df <- myL$gafnob$bgbord$dwnK27bM500
bg_gafnobg.df <- cbind(bg_gafnobg.df,NM="BG_GAFnoBG")
nobg_gafnobg.df <- myL$gafnob$nobgbord$dwnK27bM500
nobg_gafnobg.df <- cbind(nobg_gafnobg.df,NM="noBG_GAFnoBG")

perc.df <- rbind(bg_gafnobg.df,nobg_gafnobg.df) 
wdth <- 25000/(60000/2000)


perc.df$sd

pdf(paste0(OUTFD.d,"Comparison_BGbord_GAF_VS_noBGbord_GAF_N.pdf"),width=20)
ggplot(perc.df, aes(x=x , y=mean,fill=NM)) + geom_bar(position=position_dodge(1200),stat="identity",width = wdth)   +
  scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + scale_x_continuous(limits = c(0,30000)) +
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(1200),size=0.5,width=.15) + scale_fill_manual(values = c("dodgerblue","firebrick"))
dev.off()

#   -----------------------------------------------------------------------

# Figure 2F - Beaf K27 and mDOM distribution / normalized TAD borders ---------------------------------------------------------------
SUB <-  "FIGURES/Figure_2/Figure_2F/"
OUTFD.d <- create(OUT.d,SUB)
# Distance to Heterochromatin to filter out Active mDOM
DT_K27 <- 0
mDOM.sub.dm3.gr <- subset(mDOM.dm3.gr,dtk27 > DT_K27)
seqlevelsStyle(mDOM.sub.dm3.gr) <- "Ensembl"
mDOM.subD.dm3.gr <- subset(mDOM.sub.dm3.gr,decLFC > 10)
mDOM.subU.dm3.gr <- subset(mDOM.sub.dm3.gr,decLFC <= 10)
upk27.dm3.gr <- subset(mrc_tss.gr,gene_id %in% matrecap[matrecap$upK27bM500==1,"Fbgnid"])
TAD.l <- c("wang.tad.gr","ramirez.act.tad.gr","ramirez.ina.tad.gr","ramirez.pcg.tad.gr","ramirez.hp1.tad.gr")

pdf(paste0(OUTFD.d,"Distribution_plot_RELATIVE_DIST_mDOM_DOWN_AND_BEAF_DTK27",DT_K27,".pdf"))
for(TAD.c in TAD.l){ # Iterate through TADs list
  DF <- NULL
  TAD.gr <- get(TAD.c)
  TAD.gr <- addSeqinfo(TAD.gr,g_name = "dm3",celltype="kc")
  BORDL.gr <- bordersGR(TAD.gr)$leftBORD
  BORDR.gr <- bordersGR(TAD.gr)$rightBORD
  
  BEAF_res.gr <- resize(beaf.dm3.gr,1,"center")
  UPK27_res.gr <- resize(upk27.dm3.gr,1,"center")
  mDOM_res.gr <- resize(mDOM.subD.dm3.gr,1,"center")
  NM <- gsub("(.*).tad.gr","\\1",TAD.c)
  
  fol_BEAF <- findOverlaps(BEAF_res.gr,TAD.gr)
  fol_UPK27 <- findOverlaps(UPK27_res.gr,TAD.gr)
  fol_mDOM <- findOverlaps(mDOM_res.gr,TAD.gr)
  
  
  sub_BEAF_res.gr <- BEAF_res.gr[fol_BEAF@from]
  sub_BEAF_res.gr$TAD <- fol_BEAF@to  
  
  sub_UPK27_res.gr <- UPK27_res.gr[fol_UPK27@from]
  sub_UPK27_res.gr$TAD <- fol_UPK27@to  
  
  sub_mDOM_res.gr <- mDOM_res.gr[fol_mDOM@from]             
  sub_mDOM_res.gr$TAD <- fol_mDOM@to              
  
  Beaf_DST <- abs(start(sub_BEAF_res.gr) - start(BORDL.gr[sub_BEAF_res.gr$TAD]))*100/width(TAD.gr[sub_BEAF_res.gr$TAD])
  UPK27_DST <- abs(start(sub_UPK27_res.gr) - start(BORDL.gr[sub_UPK27_res.gr$TAD]))*100/width(TAD.gr[sub_UPK27_res.gr$TAD])
  mDOM_DST <- abs(start(sub_mDOM_res.gr) - start(BORDL.gr[sub_mDOM_res.gr$TAD]))*100/width(TAD.gr[sub_mDOM_res.gr$TAD])
  
  DF <- rbind.data.frame(DF,
                         cbind.data.frame(NM=paste0(NM,"_DOWN"),DST=mDOM_DST),
                         cbind.data.frame(NM=paste0(NM,"_UPK227"),DST=UPK27_DST),
                         cbind.data.frame(NM=paste0(NM,"_BEAF"),DST=Beaf_DST)
  )
  
  p <- gghistogram(DF,y="..density..",x="DST",color=NA,fill="NM",bins = 100,alpha = 0.4,palette="jco") + xlim(0, 100)
  print(p)
  p <- gghistogram(DF,y="..density..",x="DST",color=NA,fill="NM",bins = 50,alpha = 0.4,palette="jco") + xlim(0, 100)
  print(p)
  p <- gghistogram(DF,y="..density..",x="DST",color=NA,fill="NM",bins = 25,alpha = 0.4,palette="jco") + xlim(0, 100)
  print(p)
  p <- ggdensity(DF,y="..density..",x="DST",color=NA,fill="NM",alpha = 0.4,palette="jco") + xlim(0, 100)
  print(p)
  
}

dev.off()


#   -----------------------------------------------------------------------
