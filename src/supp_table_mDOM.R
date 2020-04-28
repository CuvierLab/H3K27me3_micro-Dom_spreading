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


# Table S1/S2 ----------------------------------------------------------------
SUB <- "DATA/PROCESSED/SUP/"
OUTF.d <- create(OUT.d,SUB)
MAX <- 1500
MIN <- 80

microWT_KD.gr <- subset(mDOM.dm3.gr,width < 1500 & width >= 80)
microWT_KD.gr$DST_HET <- mcols(distanceToNearest(microWT_KD.gr,Het.gr))$distance
microWT_KD.gr$DST_EUC <- mcols(distanceToNearest(microWT_KD.gr,Euc.gr))$distance
seqlevelsStyle(microWT_KD.gr) <- "Ensembl"
ctl.bwp <- bw.ll[["RAW"]][["sig"]][1]
trt.bwp <- bw.ll[["RAW"]][["sig"]][2]
microWT_KD.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microWT_KD.gr)
microWT_KD.gr <- lfc_zs_dif_decile(anc.gr = microWT_KD.gr,NTILE=NTILE,ASYM=T,nm = "_WT_KD")
microWT_KD.gr$DOWN <- 0
microWT_KD.gr$ID <- 1:length(microWT_KD.gr)
microWT_KD.gr$USED <- 0
microWT_KD.gr <- GRanges(as.data.frame(microWT_KD.gr) %>% mutate_cond(ZS_WT_KD< 0.002,DOWN=1))
MYID <- subset(microWT_KD.gr[-unique(findOverlaps(microWT_KD.gr,mrc_tss.gr[bg],maxgap=2e3,ignore.strand=T)@from)],dtk27 > 1e3 + MIN)$ID
microWT_KD.gr[MYID]$USED <- 1



# MATRECAP GENE CENTERED
mrc_tss_tmp.gr <- mrc_tss.gr
# mcols(mrc_tss_tmp.gr) <- NULL
mrc_tss_tmp.gr$mDOM <- 0
mrc_tss_tmp.gr$BEAF <- 0
mrc_tss_tmp.gr$CTCF <- 0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
mrc_tss_tmp.gr$CP190 <- 0
mrc_tss_tmp.gr$GAF <- 0
mrc_tss_tmp.gr$COH <- 0
mrc_tss_tmp.gr$DREF <- 0
mrc_tss_tmp.gr$UR <- 0


# TSS CENTERED
seqlevelsStyle(microWT_KD.gr) <- "Ensembl"
fol <- findOverlaps(mrc_tss_tmp.gr,microWT_KD.gr,maxgap=1500)
mrc_tss_tmp.gr[unique(fol@from)]$mDOM <- 1
mrc_tss_tmp.gr[bg]$BEAF <- 1
mrc_tss_tmp.gr[ctcf]$CTCF <- 1
mrc_tss_tmp.gr[cp]$CP190 <- 1
mrc_tss_tmp.gr[gaf_kc]$GAF <- 1
mrc_tss_tmp.gr[coh]$COH <- 1
mrc_tss_tmp.gr[dref]$DREF <- 1
mrc_tss_tmp.gr[ur_paul]$UR <- 1

mrc_tss_tmp.df <- as.data.frame(mrc_tss_tmp.gr)
mrc_tss_tmp.df$FBgnID <- rownames(mrc_tss_tmp.df)
# mrc_tss_tmp.df <- mrc_tss_tmp.df[,c(12,6:11)]
nrow(subset(mrc_tss_tmp.df,mDOM==1 & UR==1))
write.csv(mrc_tss_tmp.df,paste0(OUTF.d,"tss_features_s2_dm3_tableS2.csv"),row.names = F)

# MATRECAP normR CENTERED
microWT_KD.gr$BEAF <- 0
microWT_KD.gr$CTCF <- 0
microWT_KD.gr$CP190 <- 0
microWT_KD.gr$GAF <- 0
microWT_KD.gr$DREF <- 0

seqlevelsStyle(microWT_KD.gr) <- "Ensembl"
microWT_KD.gr[unique(findOverlaps(beaf.dm3.gr,microWT_KD.gr,maxgap=1e3)@to)]$BEAF <- 1
microWT_KD.gr[unique(findOverlaps(ctcf.dm3.gr,microWT_KD.gr,maxgap=1e3)@to)]$CTCF <- 1
microWT_KD.gr[unique(findOverlaps(cp190.dm3.gr,microWT_KD.gr,maxgap=1e3)@to)]$CP190 <- 1
microWT_KD.gr[unique(findOverlaps(gaf.dm3.gr,microWT_KD.gr,maxgap=1e3)@to)]$GAF <- 1
microWT_KD.gr[unique(findOverlaps(dref.dm3.gr,microWT_KD.gr,maxgap=1e3)@to)]$DREF <- 1
microWT_KD.df <- as.data.frame(microWT_KD.gr)
write.csv(microWT_KD.df,paste0(OUTF.d,"microdomains_features_s2_dm3_tableS1.csv"))


#   -----------------------------------------------------------------------
