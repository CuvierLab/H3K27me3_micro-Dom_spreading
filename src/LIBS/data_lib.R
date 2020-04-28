
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
TMP.d <- create(paste0(SPATH,"/TMP/"))


#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-


# - K27 domains
# See h3k27_domains_processing.R script
k27domains.grl <- readRDS(paste0(OUT.d,"/DATA/PROCESSED/h3k27me3_hmmxnormr_domains.rds"))
normRxHMM.gr <- k27domains.grl$raw
Het.gr <- subset(normRxHMM.gr,type=="Het")
Euc.gr <- subset(normRxHMM.gr,type=="Euc")
mcols(Het.gr) <- NULL
mcols(Euc.gr) <- NULL
seqlevelsStyle(Het.gr) <- "Ensembl"
seqlevelsStyle(Euc.gr) <- "Ensembl"
HetL_BORD.gr <- resize(Het.gr,1,"start") + 1500
HetR_BORD.gr <- resize(Het.gr,1,"end") + 1500


# K27 MULTI NORMALIZATION
# Published data
h3k27kd_RPKM.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/K27_KDbeaf_2013_Q10_sorted_SHFT130_RPKM.bw")
h3k27wt_RPKM.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/K27_WT_2013_Q10_sorted_SHFT130_RPKM.bw")
h3k27kd_raw.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_rawSignal_KD_cuvier_s2_dm3.bw")
h3k27wt_raw.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_rawSignal_WT_cuvier_s2_dm3.bw")
h3k27kd_log.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/K27_KDbeaf_2013_Q10_sorted_None.bw")
h3k27wt_log.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/K27_WT_2013_Q10_sorted_None.bw")
h3k27kd_norm.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_normSeqDepth_KD_cuvier_s2_dm3.bw")
h3k27wt_norm.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_normSeqDepth_WT_cuvier_s2_dm3.bw")
h3k27wt_normH3.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_LFC_normH3_normSeqDepth_wtxh3_cuvier_s2_dm3.bw")
h3k27kd_normH3.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/chipK27_LFC_normH3_normSeqDepth_kdxh3_cuvier_s2_dm3.bw")
nuc_lfc.bwp <- paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/H3K27me3/mnase_LFC_kdxwt_cuvier_s2_dm3.bw")

bw.ll <- list(RPKM=list(sig=c(h3k27wt_RPKM.bwp,h3k27kd_RPKM.bwp),ylim=c(-10,50)),
              RAW=list(sig=c(h3k27wt_raw.bwp,h3k27kd_raw.bwp),ylim=c(0,0.5)),
              SQDNORM=list(sig=c(h3k27wt_norm.bwp,h3k27kd_norm.bwp),ylim=c(0,0.2)),
              NORMH3=list(sig=c(h3k27wt_normH3.bwp,h3k27kd_normH3.bwp,nuc_lfc.bwp),ylim=c(-5,5)),
              NONE=list(sig=c(h3k27wt_log.bwp,h3k27kd_log.bwp),ylim=c(-0,5)))




# MICRO DOMAINS
# See h3k27_domains_processing.R script
mDOM.dm3.gr <- addSeqinfo(readRDS(paste0(OUT.d,"/DATA/PROCESSED/microdomains_normr.dm3.rds")))
mDOM.dm3.gr$dtk27 <- mcols(distanceToNearest(mDOM.dm3.gr,Het.gr))$distance
mDOM.dm3.gr$dteuc <- mcols(distanceToNearest(mDOM.dm3.gr,Euc.gr))$distance
seqlevelsStyle(mDOM.dm3.gr) <- "Ensembl"
export.bed(mDOM.dm3.gr,paste0(OUT.d,"/DATA/PROCESSED/MDOM/WT_INPUT_microdomains_normr.dm3.bed"))



# TADS
corces.tad.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/corces_TADS.dm3.rds"))
corces.act.tad.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/corces_active_TADS.dm3.rds"))
eagen.tad.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/eagen_TAD.dm3.rds"))
eagen.act.tad.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/eagen_active_TADS.dm3.rds"))
wang.tad.gr <- addSeqinfo(readRDS(paste0(OUT.d,"/DATA/PUBLISHED/wang_TADS.dm3.rds")))
ramirez.all.tad.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/ramirez_TADS.dm3.rds"))
ramirez.act.tad.gr <- subset(ramirez.all.tad.gr,name=="active")
ramirez.ina.tad.gr <- subset(ramirez.all.tad.gr,name=="inactive")
ramirez.pcg.tad.gr <- subset(ramirez.all.tad.gr,name=="PcG")
ramirez.hp1.tad.gr <- subset(ramirez.all.tad.gr,name=="HP1")

# HIC
HIC.l <-  list(WT=paste0(OUT.d,"/DATA/PUBLISHED/Ramirez_inter_1000_allCHR.rds"),
               KD=paste0(OUT.d,"/DATA/PUBLISHED/Ramirez_BeafKD_inter_1000_allCHR.rds"))
mydothiclist <- list(ramirezS2bKD_1000=paste0(OUT.d,"/DATA/PUBLISHED/Ramirez_bKD_inter_1000.hic"),
                  ramirezS2WT_1000=paste0(OUT.d,"/DATA/PUBLISHED/Ramirez_WT_inter_1000.hic"))


#chipseq
beaf.dm3.gr <- addSeqinfo(GenomicRanges::reduce(import(paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/beaf32.s2.dm3.bed"))))
k27ac.dm3.gr <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/h3k27ac_Kc167_dm3.gr.RData"))
pc.dm3.gr <-  addSeqinfo(import(paste0(OUT.d,"/DATA/PUBLISHED/CHIPSEQ/GSE55303_PC.bed")))

# Compartments
ab_compRamirezWT.bwp <- paste0(OUT.d,"/PROCESSING_SCRIPTS/COMPARTMENT_COMPUTATION/Compartments_REANNOTATED/beafWT/Ramirez_compartment_res10000_S2_dm3.wig_ab_comp_res10000_S2_dm3.gr.bw")
