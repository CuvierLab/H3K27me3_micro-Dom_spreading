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
source(paste0(OUT.d,"/LIBS/mrc_lib.R"))
TMP.d <- create(paste0(SPATH,"/TMP/"))


#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-
# F7B - PROFILE ---------------------------------------------------------------
SUB <- "FIGURES/Figure_7/Figure_7B/"
OUTF.d <- create(OUT.d,SUB)


BS_CP <- 200
FDR_CP <- 0.01
BS_K27 <- 250
FDR_K27 <- 0.01
TYPE <- "RPKM"
DECTYPE <- "ZS"
NTILE <- 10
MAXGAP <- 2000
asymetric <- T
SD <- 12

# BWP
trt.bwpl <- list(RPKM=paste0(OUT.d,"DATA/PROCESSED/ChIPseq/H3K27me3/MUTANT/h3k27me3_MBImerge_Q10_sorted_dedup_SHFT140_RPKM.bw"))
ctl.bwpl <- list(RPKM=paste0(OUT.d,"DATA/PROCESSED/ChIPseq/H3K27me3/WT/h3k27me3_BImerge_Q10_sorted_dedup_SHFT135_RPKM.bw"))
YLIM <- list(RPKM=c(0,120))
Ylim <- YLIM[[TYPE]]

trt.bwp <- trt.bwpl[[TYPE]]
ctl.bwp <- ctl.bwpl[[TYPE]]

# DOMAIN LIST
microMBI_BI.gr <- GRanges(import(paste0("DATA/PROCESSED/MDOM/MBIvsBI_BIN=",BS_K27,"_FDR=",FDR_K27,"_diffR.bed")))
microMBI_BI.gr <- addSeqinfo(microMBI_BI.gr)
microMBI_BI.gr$dtk27 <- mcols(distanceToNearest(microMBI_BI.gr,Het.gr))$distance
microMBI_BI.gr <- subset(microMBI_BI.gr,dtk27 > 5000)
seqlevelsStyle(microMBI_BI.gr) <- "UCSC"

# COMPUTE SIGNAL FOR K27 CHIP

microMBI_BI.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microMBI_BI.gr)
microMBI_BI.gr <- lfc_zs_dif_decile(anc.gr = microMBI_BI.gr,NTILE=NTILE,ASYM=T)
seqlevelsStyle(microMBI_BI.gr) <- "Ensembl"
microMBI_BI_DOWNMBI.gr <- subset(microMBI_BI.gr,decLFC > NTILE)
microMBI_BI_UPMBI.gr <- subset(microMBI_BI.gr,decLFC <= NTILE)


# CP
cp190_WT.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041448_SRR1041449_Q10_sorted_dedup_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))
cp190_MUT.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/MUTANT_BEAF/SRR1041451_1_Q10_sorted_dedup_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))
cp190_WT_MBI.gr <- reduce(c(cp190_WT.gr))
cp190_WT_MBI.gr <- addSeqinfo(cp190_WT_MBI.gr)
seqlevelsStyle(cp190_WT_MBI.gr) <- "Ensembl"



SD <- 12
OUT_PROF.d <- create(OUT.d,"PROFILE/PROFILE_microDOMAINS_BI_MBI_CENTERED_BOUND_CP190/BS_CP",BS_CP,
                     "/FDRCP_",FDR_CP,"/BSK27_",BS_K27,"/FDRK27_",FDR_K27,
                     "/BW_TYPE_",TYPE,"_DECTYPE_",DECTYPE,
                     "/MAXGAP_",MAXGAP,"/")

microMBI_BI_CP.gr <- subsetByOverlaps(microMBI_BI_DOWNMBI.gr,cp190_WT_MBI.gr,maxgap = MAXGAP)

pdf(paste0(OUTF.d,"PROFILE_K27_BI_MBI_DOWNvsUP_MBI_K27_BOUND_CP.pdf"))
bw.l <- c(trt.bwp,ctl.bwp)
seqPlotSDoutliers(bw.l,TMP.d,c("microMBI_BI_CP.gr"),ylim=Ylim,xlim=c(1000,1000),bin=10L,sd=SD,err = F,smooth = F,type="pf",gnme = "dm3",ignore.strand = F,main="MICRODOM CENTERED DOWN in MBI")
dev.off()

#   -----------------------------------------------------------------------


# FIGURE 7C ---------------------------------------------------------------
# /media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Manuscript/INTEGRATION_CP190_MBI_BI_WT_KD/BOXPLOT/
#BOXPLOT_K27_MBI_BI_RANKED_CP190_DIFF_MBI_BI/NTILE_200/BS_CP200/FDRCP_0.01/BSK27_250/FDRK27_0.01/MAXGAP_2000/
# BOXPLOT_K27_MBI_BI_RANKED_CP190_DIFF_MBI_BI_BSCP200_FDRCP_0.01BSK27_250_FDRK27_0.01_MAXGAP_2000.pdf

SUB <- "FIGURES/Figure_7/Figure_7C/"
OUTF.d <- create(OUT.d,SUB)

BS_CP <- 200
FDR_CP <- 0.01
BS_K27 <- 250
FDR_K27 <- 0.01
TYPE <- "RPKM"
NTILE <- 200
MAXGAP <- 2000
asymetric <- F
type <- "INPUTNORM"

cp190_WT1.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041448_1_Q10_sorted_dedup.bw")
cp190_WT2.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041449_1_Q10_sorted_dedup.bw")
cp190_WT.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041448_SRR1041449_Q10_sorted_dedup_RPKM_log2_inputReadCountNorm.bw")
cp190_MUT.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/MUTANT_BEAF/SRR1041451_1_Q10_sorted_dedup_RPKM_log2_inputReadCountNorm.bw")

h3k27mbi.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/H3K27me3/MUTANT/h3k27me3_MBImerge_Q10_sorted_dedup_RPKM.bw")
h3k27bi.bwp <- paste0(OUT.d,"DATA/PROCESSED/ChIPseq/H3K27me3/WT/h3k27me3_BImerge_Q10_sorted_dedup_RPKM.bw")


# CP
cp190_WT1.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041448_1_Q10_sorted_dedup.bai_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))
cp190_WT2.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041449_1_Q10_sorted_dedup.bai_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))
cp190_WT.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/WT/SRR1041448_SRR1041449_Q10_sorted_dedup_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))
cp190_MUT.gr <- import(paste0(OUT.d,"DATA/PROCESSED/ChIPseq/CP190/MUTANT_BEAF/SRR1041451_1_Q10_sorted_dedup_BIN=",BS_CP,"_FDR=",FDR_CP,"_normR.bed"))

# MICRO DOMAIN MBI/BI
microMBI_BI.gr <- GRanges(import(paste0(OUT.d,"DATA/PROCESSED/MDOM/MBIvsBI_EUC_BIN=",BS_K27,"_FDR=",FDR_K27,"_diffR.bed")))
microMBI_BI.gr <- addSeqinfo(microMBI_BI.gr)
microMBI_BI.gr$dtk27 <- mcols(distanceToNearest(microMBI_BI.gr,Het.gr))$distance
seqlevelsStyle(microMBI_BI.gr) <- "UCSC"

mbi.cov <- importData(h3k27mbi.bwp,"BigWig",microMBI_BI.gr)
mbi.cov <- mbi.cov[[1]]
microMBI_BI.gr$MBIrpkm <- rowMeans(as.matrix(mbi.cov[microMBI_BI.gr]),na.rm = T)
bi.cov <- importData(h3k27bi.bwp,"BigWig",microMBI_BI.gr)
bi.cov <- bi.cov[[1]]
microMBI_BI.gr$BIrpkm <- rowMeans(as.matrix(bi.cov[microMBI_BI.gr]),na.rm = T)
microMBI_BI.gr$LFCrpkm <- log2(microMBI_BI.gr$MBIrpkm / microMBI_BI.gr$BIrpkm)
microMBI_BI.gr <- GRanges(as.data.frame(microMBI_BI.gr) %>% mutate(decile=ntile(LFCrpkm,NTILE)))
if(asymetric==T){
  microMBI_BI_sup.gr <- subset(microMBI_BI.gr,LFCrpkm >= 0)
  microMBI_BI_min.gr <- subset(microMBI_BI.gr,LFCrpkm < 0)
  microMBI_BI_sup.gr <- GRanges(as.data.frame(microMBI_BI_sup.gr) %>% mutate(decile=(NTILE+1)-ntile(LFCrpkm,NTILE)))
  microMBI_BI_min.gr <- GRanges(as.data.frame(microMBI_BI_min.gr) %>% mutate(decile=(NTILE*2+1)-ntile(LFCrpkm,NTILE)))
  microMBI_BI.gr <- c(microMBI_BI_sup.gr,microMBI_BI_min.gr)
  NTILE <- NTILE*2
}
seqlevelsStyle(microMBI_BI.gr) <- "Ensembl"
# LET'S FIRST TRY ON MERGE DATASET FOR CP190
cp190_WT_MBI.gr <- reduce(c(cp190_WT.gr,cp190_MUT.gr))
# cp190_WT_MBI.gr <- subsetByOverlaps(cp190_WT_MBI.gr,myMRC.tss.gr,maxgap=200)
cp_mbi.cov <- importData(cp190_MUT.bwp,"BigWig",cp190_WT_MBI.gr)
cp_mbi.cov <- cp_mbi.cov[[1]]
cp190_WT_MBI.gr$MBIrpkm <- rowMeans(as.matrix(cp_mbi.cov[cp190_WT_MBI.gr]),na.rm = T)
cp_wt.cov <- importData(cp190_WT.bwp,"BigWig",cp190_WT_MBI.gr)
cp_wt.cov <- cp_wt.cov[[1]]
cp190_WT_MBI.gr$BIrpkm <- rowMeans(as.matrix(cp_wt.cov[cp190_WT_MBI.gr]),na.rm = T)
cp190_WT_MBI.gr$LFCrpkm <- cp190_WT_MBI.gr$MBIrpkm - cp190_WT_MBI.gr$BIrpkm
cp190_WT_MBI.gr <- GRanges(as.data.frame(cp190_WT_MBI.gr) %>% mutate(decile=ntile(LFCrpkm,NTILE)))
cp190_WT_MBI.gr <- addSeqinfo(cp190_WT_MBI.gr)

cp190_WT_MBI.gr$decile <- ifelse(cp190_WT_MBI.gr$decile<15,1,
                                 ifelse(cp190_WT_MBI.gr$decile %in% c(16,17,18),2,3))

GRP <- sort(unique(cp190_WT_MBI.gr$decile))
######################################### -
# PLOT
######################################### -

pdf(paste0(OUTF.d,"BOXPLOT_K27_MBI_BI_RANKED_CP190_DIFF_MBI_BI_BSCP",BS_CP,"_FDRCP_",FDR_CP,
           "BSK27_",BS_K27,"_FDRK27_",FDR_K27,"_MAXGAP_",MAXGAP,".pdf"),width=20)
# DST <- 15e3
DST <- 15e3 
DF <- NULL
for(nt in GRP){
  cp190_sub.gr <- subset(cp190_WT_MBI.gr,decile==nt)
  microMBI_BI_sub.gr <- subset(microMBI_BI.gr,dtk27>=DST)
  microMBI_BI_sub.gr <- subsetByOverlaps(microMBI_BI_sub.gr,cp190_sub.gr,maxgap=MAXGAP)
  DF <- rbind.data.frame(DF,
                         cbind.data.frame(WT=microMBI_BI_sub.gr$BIrpkm,KD=microMBI_BI_sub.gr$MBIrpkm,DEC=nt))
  
}
DFm <- melt(DF,id.vars = 3)
DFm$DEC <- as.factor(DFm$DEC)    
p <- ggplot(DFm, aes(x = DEC, y = value,color = variable),
            size=2,width=0.2) +
  coord_cartesian(ylim = c(-20,150)) +
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=1.2,width = 0.1) + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "linerange", size=1.2,linetype="dashed") + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=1.2,width = 0.1) +
  geom_boxplot(width=0.3,outlier.size = 0, coef = 0,size=1.2,outlier.shape=NA,position=position_dodge(0.6)) +
  stat_compare_means(aes(group = variable),paired=T,label.y =-4,size=4,label.x = 1.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     label = "..p.signif..") + 
  scale_color_grey() + theme_classic() + theme(plot.title = element_text(size = 10, face = "bold"))+
  ggtitle(paste0("Boxplot - Wilcoxon TEST \n MBI vs BI levels DST K27=",DST," ranked by CP190 differential binding \n (1 = WT -> 6=MBI)"))
print(p)
dev.off()


pdf(paste0(OUTF.d,"BOXPLOT_CP190_MBI_BI_RANKED_CP190_DIFF_MBI_BI_BSCP",BS_CP,"_FDRCP_",FDR_CP,
           "BSK27_",BS_K27,"_FDRK27_",FDR_K27,"_MAXGAP_",MAXGAP,".pdf"),width=20)
for(nt in GRP){
  cp190_sub.gr <- subset(cp190_WT_MBI.gr,decile==nt)
  DF <- rbind.data.frame(DF,
                         cbind.data.frame(WT=cp190_sub.gr$BIrpkm,KD=cp190_sub.gr$MBIrpkm,DEC=nt))
  
}
DFm <- melt(DF,id.vars = 3)
DFm$DEC <- as.factor(DFm$DEC)    
p <- ggplot(DFm, aes(x = DEC, y = value,color = variable),
            size=2,width=0.2) +
  coord_cartesian(ylim = c(-1.0,3)) +
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=1.2,width = 0.1) + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "linerange", size=1.2,linetype="dashed") + 
  stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.25,na.rm=T)-0.2*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
               geom = "crossbar", size=1.2,width = 0.1) +
  geom_boxplot(width=0.3,outlier.size = 0, coef = 0,size=1.2,outlier.shape=NA,position=position_dodge(0.6)) +
  stat_compare_means(aes(group = variable),paired=T,label.y =-4,size=4,label.x = 1.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     label = "..p.signif..") + 
  scale_color_grey() + theme_classic() + theme(plot.title = element_text(size = 10, face = "bold"))+
  ggtitle(paste0("Boxplot - Wilcoxon TEST \n MBI vs BI levels DST K27=",DST," ranked by CP190 differential binding \n (1 = WT -> 6=MBI)"))
print(p)
dev.off()



#   -----------------------------------------------------------------------


# Figure 7 D----------------------------------------------------------------

SUB <- "FIGURES/Figure_7/Figure_7D/"
OUTF.d <- create(OUT.d,SUB)

BS_K27 <- 40
FDR_K27 <- 0.1
BS_K27_MUT <- 250
FDR_K27_MUT <- 0.1
TYPE <- "RPKM"
SQLVSTYLE <- "UCSC" 
DECTYPE_CHIP <- "ZS"
SPREAD <- "BOTH"
NTILE <- 10
baseSet <- "baseSet_TSS"
asymetric <- T
DST_BORDER <- 10000
bInf <- 40
bSup <- 1500
XMIN <- 2000
XMAX <- 2000
microDOMTYPE <- "WT_KD"



# RNASEQ
rna.lst <- list(WT=matrecap$X.WTPaul+1,KD=matrecap$X.KDBPaul+1)
rnaCOUNTS.lst <- list(UR=ur_paul,DR=dr_paul)

mrc_tss_tmp.gr <- mrc_tss.gr


# K27 CHIPSEQ
# BWP
YLIM <- bw.ll[[TYPE]][["ylim"]]
ctl.bwp <- bw.ll[[TYPE]][["sig"]][1]
trt.bwp <- bw.ll[[TYPE]][["sig"]][2]

# BWL
BW.l <- c(trt.bwp,ctl.bwp)

# MICRO DOMAINS WT/KD
microKD_WT_ALL.gr <- mDOM.dm3.gr
microKD_WT.gr <- subset(microKD_WT_ALL.gr,dtk27 > DST_BORDER & width < bSup & width >= bInf)
microKD_WT_BORDER.gr <- subset(microKD_WT_ALL.gr,dtk27 < 1000 & width < bSup & width >= bInf)
seqlevelsStyle(microKD_WT.gr) <- SQLVSTYLE
seqlevelsStyle(microKD_WT_BORDER.gr) <- SQLVSTYLE

# COMPUTE SIGNAL FOR K27 CHIP
microKD_WT.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microKD_WT.gr)
microKD_WT.gr <- lfc_zs_dif_decile(anc.gr = microKD_WT.gr,NTILE=NTILE,ASYM=T)
seqlevelsStyle(microKD_WT.gr) <- "Ensembl"
microKD_WT_DOWNKD.gr <- GRanges(as.data.frame(microKD_WT.gr) %>% filter(UQ(as.name(paste0("dec",DECTYPE_CHIP))) > NTILE + 3))
microKD_WT_UPKD.gr <- GRanges(as.data.frame(microKD_WT.gr) %>% filter(UQ(as.name(paste0("dec",DECTYPE_CHIP))) <= NTILE-6))
microKD_WT_noUPDOWN.gr <- GRanges(as.data.frame(microKD_WT.gr) %>% filter(UQ(as.name(paste0("dec",DECTYPE_CHIP))) %in% ((NTILE-5) : (NTILE +5))))
microKD_WT_ALL.gr <- microKD_WT.gr

mrc_tss_tmp.gr$dtk27 <- 0
dtn <- distanceToNearest(mrc_tss_tmp.gr,Het.gr)
mrc_tss_tmp.gr[dtn@from]$dtk27 <- mcols(dtn)$distance
mrc_tss_tmp.gr$ID <- 1:length(mrc_tss_tmp.gr)


DF <- NULL
MAXGAPTSS <- 1e3  
tss.IDX <- mrc_tss_tmp.gr$ID
rnaLV <- log2((rna.lst[["KD"]][tss.IDX])/(rna.lst[["WT"]][tss.IDX]))
DF <- rbind.data.frame(DF,cbind.data.frame(type="ALL_TSS",val=rnaLV))

# 
tssNO.IDX <- mrc_tss_tmp.gr[nobg & nocp & !bs & nogaf_kc]$ID
rnaLV <- log2((rna.lst[["KD"]][tssNO.IDX])/(rna.lst[["WT"]][tssNO.IDX]))
DF <- rbind.data.frame(DF,cbind.data.frame(type="TSS_NOINSUL",val=rnaLV))


tssBG_CP.IDX <- mrc_tss_tmp.gr[bg & matrecap$dist2k27Bg < 1000 & cp & upK27bPM1kb]$ID
rnaLV <- log2((rna.lst[["KD"]][tssBG_CP.IDX])/(rna.lst[["WT"]][tssBG_CP.IDX]))
DF <- rbind.data.frame(DF,cbind.data.frame(type="TSS_BG_CP",val=rnaLV))

myFOL <- findOverlaps(microKD_WT_DOWNKD.gr,mrc_tss_tmp.gr,maxgap=MAXGAPTSS)
tss_micro.IDX <- subset(mrc_tss_tmp.gr,ID %in% subjectHits(myFOL))$ID
rnaLV <- log2((rna.lst[["KD"]][tss_micro.IDX])/(rna.lst[["WT"]][tss_micro.IDX]))
DF <- rbind.data.frame(DF,cbind.data.frame(type="TSS_MICRO_DOWN",val=rnaLV))


tss_gaf_micro.IDX <- subset(mrc_tss_tmp.gr,ID %in% subjectHits(myFOL) & gaf_kc)$ID
rnaLV <- log2((rna.lst[["KD"]][tss_gaf_micro.IDX])/(rna.lst[["WT"]][tss_gaf_micro.IDX]))
DF <- rbind.data.frame(DF,cbind.data.frame(type="TSS_MICRO_DOWN_GAF",val=rnaLV))




pdf(paste0(OUTF.d,"BOXPLOT_RNASEQ_WT_KD_ALL_TSS_AND_MICRODOM_TSS_WTKD_MAXGAPTSS_",MAXGAPTSS,".pdf"))
print(ggboxplot(DF,x="type",y="val",col="type") + theme(axis.text.x = element_text(angle=90)) +
        stat_compare_means(paired=F,method.args = list(alternative="greater"),
                           ref.group = "ALL_TSS",
                           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                           label = "..p.signif.."))
dev.off()




#   -----------------------------------------------------------------------



# FIGURE 7E ---------------------------------------------------------------
SUB <- "FIGURES/Figure_7/Figure_7E/"
OUTF.d <- create(OUT.d,SUB)


feat.l <- c("bg|bs")
of.dir <- create(paste0(OUTF.d,"Features_binding_influences_whole_domains_in_K27_LV/"))

all_nocp_nobg <- all & !cp & !bg
`bg|cp|bs` <- bg | cp | bs
`bg|bs` <- bg | bs
nobgorcporbs <- !(bg | cp | bs)
`cp&noBg&Bs` <- cp & !bg & bs
`cp&noBg` <- cp & !bg

YLIM <- c(-8,8)

for(mxgap in c(500,1200,1500)){
  for(xtend in c(1000)){
    all_dom.gr <- k27domains.grl$h3SqdNorm
    all_dom.gr$wdth <- width(all_dom.gr)
    eucdom <- subset(all_dom.gr,type=="Euc" & wdth < 2000) 
    mcols(eucdom)[,5:(ncol(mcols(eucdom))-1)] <- NA 
    all_dom.gr[which(all_dom.gr$type == "Euc" & all_dom.gr$wdth < 2000)] <- eucdom
    all_dom.gr$len <- 1:length(all_dom.gr)
    all_dom.gr$nb <- 0
    all_dom.gr$cnt <- 0
    all_dom.gr$Lcnt <- 0
    all_dom.gr$Rcnt <- 0
    all_dom.grl <- split(all_dom.gr,seqnames(all_dom.gr))
    # Avoid first/last domain of each chromosome, issue xhen doing +/- 1
    all_dom_noend.grl <- lapply(all_dom.grl,function(domChrX.gr){ domChrX.gr[-c(1,length(domChrX.gr))] })
    all_dom_noend.gr <- Reduce("c",all_dom_noend.grl)
    het_noend.gr <- all_dom_noend.gr[all_dom_noend.gr$type=="Het"]
    euc_noend.gr <- all_dom_noend.gr[all_dom_noend.gr$type=="Euc"]
    # Dataframe containing the results that will be plotted
    
    for(feat in feat.l){
      myDf <- NULL
      hetDom.gr <- all_dom_noend.gr[all_dom_noend.gr$type=="Het"]
      eucDom.gr <- all_dom_noend.gr[all_dom_noend.gr$type=="Euc"]
      # Constraint domains to have a TSS at each border
      hetDomL.gr <- resize(hetDom.gr,1,"start"); hetDomR.gr <- resize(hetDom.gr,1,"end")
      eucDomL.gr <- resize(eucDom.gr,1,"start"); eucDomR.gr <- resize(eucDom.gr,1,"end")
      toKeepHet <- intersect(unique(queryHits(findOverlaps(hetDomL.gr,mrc_tss.gr,maxgap=mxgap))),
                             unique(queryHits(findOverlaps(hetDomR.gr,mrc_tss.gr,maxgap=mxgap))))
      toKeepEuc <- intersect(unique(queryHits(findOverlaps(eucDomL.gr,mrc_tss.gr,maxgap=mxgap))),
                             unique(queryHits(findOverlaps(eucDomR.gr,mrc_tss.gr,maxgap=mxgap))))
      hetDomL.gr <- hetDomL.gr[toKeepHet];hetDomR.gr <- hetDomR.gr[toKeepHet]
      eucDomL.gr <- eucDomL.gr[toKeepEuc];eucDomR.gr <- eucDomR.gr[toKeepEuc]
      hetDom.gr <- hetDom.gr[toKeepHet]
      eucDom.gr <- eucDom.gr[toKeepEuc]
      
      qhHetL <- unique(queryHits(findOverlaps(hetDomL.gr,mrc_tss.gr[get(feat)],maxgap=mxgap)))
      qhHetR <- unique(queryHits(findOverlaps(hetDomR.gr,mrc_tss.gr[get(feat)],maxgap=mxgap)))
      qhEucL <- unique(queryHits(findOverlaps(eucDomL.gr,mrc_tss.gr[get(feat)],maxgap=mxgap)))
      qhEucR <- unique(queryHits(findOverlaps(eucDomR.gr,mrc_tss.gr[get(feat)],maxgap=mxgap)))
      
      hetDom.gr[qhHetL]$Lcnt <- 1 
      hetDom.gr[qhHetR]$Rcnt <- 1 
      eucDom.gr[qhEucL]$Lcnt <- 1 
      eucDom.gr[qhEucR]$Rcnt <- 1 
      
      eucIDX <- eucDom.gr$len
      hetIDX <- hetDom.gr$len
      
      
      myDf <- rbind.data.frame(myDf,
                               cbind.data.frame(WTk27dens=all_dom.gr[eucIDX]$WT_domain_density,
                                                KDk27dens=all_dom.gr[eucIDX]$KD_domain_density,
                                                group_name="EUC",
                                                dens=paste0(eucDom.gr$Lcnt,eucDom.gr$Rcnt),
                                                feat=feat),
                               cbind.data.frame(WTk27dens=all_dom.gr[hetIDX]$WT_domain_density,
                                                KDk27dens=all_dom.gr[hetIDX]$KD_domain_density,
                                                group_name="HET",
                                                dens=paste0(hetDom.gr$Lcnt,hetDom.gr$Rcnt),
                                                feat=feat))
      
      nof_nof <- paste0("∅",toupper(feat),"_∅",toupper(feat))
      f_nof <- paste0(toupper(feat),"_∅",toupper(feat))
      nof_f <- paste0("∅",toupper(feat),"_",toupper(feat))  
      f_f <- paste0(toupper(feat),"_",toupper(feat))
      myDf$dens <- as.character(myDf$dens)
      myDf <- myDf %>% mutate_cond(dens == "00", dens=nof_nof)
      myDf <- myDf %>% mutate_cond(dens == "01", dens=nof_f)
      myDf <- myDf %>% mutate_cond(dens == "10", dens=f_nof)
      myDf <- myDf %>% mutate_cond(dens == "11", dens=f_f)
      myDf$dens <- as.factor(myDf$dens)
      myDf$dens <- factor(myDf$dens, levels = c(nof_nof, f_nof,nof_f, f_f))
      myDfmelt <- melt(myDf,id=c("group_name","dens","feat"))
      
      
      for(mygp in unique(myDfmelt$group_name)){
        cairo_pdf(paste0(OUTF.d,"Boxplot_K27WTLV_vs_K27kdLV_Center_",mygp,"_Gap=",mxgap,"_",toupper(feat),".pdf"),onefile = T,family="Arial Unicode MS",height=15,width=20)
        myDfrmeltCGp <- subset(myDfmelt,group_name==mygp)
        myDfrmeltCGp <- as.data.frame(myDfrmeltCGp %>%
                                        group_by(dens,variable) %>%
                                        mutate(med=median(value,na.rm=T))) 
        tmp.df <- as.data.frame(myDfrmeltCGp %>% group_by(dens) %>% summarise(rat=unique((med[variable=="KDk27dens"]+10)-(med[variable=="WTk27dens"]+10))))
        p <- ggplot(myDfrmeltCGp, aes(x = dens, y = value,color = variable),
                    size=2,width=0.2) +
          coord_cartesian(ylim = YLIM) + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "crossbar", size=2,width = 0.2) + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "linerange", size=2,linetype="dashed") + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "crossbar", size=2,width = 0.2) +
          geom_boxplot(width=0.2,outlier.size = 0, coef = 0,size=2,outlier.shape=NA,position=position_dodge(0.6)) +
          stat_compare_means(aes(group = variable),paired=T,label.y = 6.5,size=10,label.x = 1.5,label = "..p.signif..") + 
          geom_text(data=tmp.df, aes(label=paste0("Median Difference \n ",round(rat,2)), x=dens,y=4.8),size=8, inherit.aes = F,position=position_dodge(0.9)) +
          scale_color_grey() + theme_classic() + ggtitle(paste0("CENTER=",mygp)) + theme(plot.title = element_text(size = 40, face = "bold"))
        print(p)
        
        p <- ggplot(myDfrmeltCGp, aes(x = dens, y = value,color = variable),
                    size=2,width=0.2) +
          coord_cartesian(ylim = YLIM) + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "crossbar", size=2,width = 0.2) + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.75,na.rm=T)+1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "linerange", size=2,linetype="dashed") + 
          stat_summary(fun.y = median, fun.ymin = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), fun.ymax = function(x)quantile(x,0.25,na.rm=T)-1.5*(quantile(x,0.75,na.rm=T)-quantile(x,0.25,na.rm=T)), position = position_dodge(0.6),
                       geom = "crossbar", size=2,width = 0.2) +
          geom_boxplot(width=0.2,outlier.size = 0, coef = 0,size=2,outlier.shape=NA,position=position_dodge(0.6)) +
          stat_compare_means(aes(group = variable),paired=T,label.y =5.5,size=4,label.x = 1.5,label = "..p.format..") + 
          scale_color_grey() + theme_classic() + ggtitle(paste0("CENTER=",mygp)) + theme(plot.title = element_text(size = 40, face = "bold"))
        print(p)
        
        dev.off()
      }
      
    }
    
  }
}



#   -----------------------------------------------------------------------


