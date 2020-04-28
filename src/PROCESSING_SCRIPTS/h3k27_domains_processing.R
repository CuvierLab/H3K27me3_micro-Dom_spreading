
#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          PATH  ----
#####################################################################################-
SPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
OUT.d <- create(paste0(SPATH,"/"))
TMP.d <- create(paste0(SPATH,"/TMP/"))
#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source(paste0(OUT.d,"/LIBS/lrc_lib.R"))
source(paste0(OUT.d,"LIBS/mrc_lib.R"))
source(paste0(OUT.d,"/LIBS/data_lib.R")) # load data


#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-


trt <- paste0(OUT.d,"DATA/RAW/GA0750_H3K27me3_WT_filtered_sorted_nodup.bam")
ctl <- paste0(OUT.d,"DATA/RAW/input_FCC43F9ACXX_L6_WHAPPI001534-129_1_filtered_sorted_nodup.bam")

# PROCESS HMM DOMAINS OF H3K27me3 -----------------------------------------
domK27pm.gr <- readRDS(paste0(OUT.d,"./DATA/PUBLISHED/cuvier_HMM_K27me3_domains.RData"))
K27_dom.gr <- resize(domK27pm.gr,width(domK27pm.gr)+1250,"start") # We have seen that there is a shift at 3' end
K27_dom.gr <- fill_uncovered_domain(K27_dom.gr,paste0(OUT.d,"/DATA/OTHERS/iHMM.M1K16.fly_L3.bed"))
# Control with the number of Beaf fallin in right border
domK27L.gr <- resize(K27_dom.gr,1,"start")
domK27R.gr <- resize(K27_dom.gr,1,"end")
# Use the dm3 granges to make a set diff to obtain noK27 domains (k27.gr + nok27.gr = allGenome.gr)
noK27_dom.gr <- trim(GenomicRanges::setdiff(dm3_mychr.gr,K27_dom.gr))

# Save all
K27_dom.gr$type <- "Het"
noK27_dom.gr$type <- "Euc"
allDOM.gr <- c(K27_dom.gr,noK27_dom.gr)
allDOM.gr <- sortSeqlevels(allDOM.gr);allDOM.gr <- sort(allDOM.gr)

h3k27wt.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/WT/bw/chipK27_rawSignal_WT_cuvier_s2_dm3.cov.RData")
h3k27kd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27//H3K27_WT_bKD_CUVIER/Beaf32_KD/bw/chipK27_rawSignal_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/WT/bw/chipK27_normSeqDepth_WT_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/Beaf32_KD/bw/chipK27_normSeqDepth_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/WT/bw/chipK27_LFC_normH3_normSeqDepth_wtxh3_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/Beaf32_KD/bw/chipK27_LFC_normH3_normSeqDepth_kdxh3_cuvier_s2_dm3.cov.RData")
# Load K27 domains (find by HMM and shiftd from 1300bp on right side (bias from HMM see parametrization))
# allDOM.gr <- readRDS(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Manuscript/Parameterization/k27domainrecap_generator_HMM/current/",full.names = T,pattern=".RData"))
seqlevelsStyle(mrc_tss.gr) <- "Ensembl"
allDOM.gr$k27KDdens <- 0
allDOM.gr$k27WTdens <- 0
# Give each domain a density of K27 level
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, list("allDOM.gr","h3k27kd.cov","h3k27wt.cov","h3k27wt_sqd.cov","h3k27kd_sqd.cov","h3k27wt_sqd_h3.cov","h3k27kd_sqd_h3.cov"),envir = environment())
clusterEvalQ(cl, library(GenomicRanges))
kdDens <- parSapply(cl,1:length(allDOM.gr),
                    function(idx){
                      sum(as.matrix(h3k27kd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                      # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                    })
wtDens <- parSapply(cl,1:length(allDOM.gr),
                    function(idx){
                      sum(as.matrix(h3k27wt.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                      # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                    })
kdDens_sqd <- parSapply(cl,1:length(allDOM.gr),
                        function(idx){
                          sum(as.matrix(h3k27kd_sqd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                          # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                        })
wtDens_sqd <- parSapply(cl,1:length(allDOM.gr),
                        function(idx){
                          sum(as.matrix(h3k27wt_sqd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                          # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                        })
kdDens_sqd_h3 <- parSapply(cl,1:length(allDOM.gr),
                           function(idx){
                             sum(as.matrix(h3k27kd_sqd_h3.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                             # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                           })

wtDens_sqd_h3 <- parSapply(cl,1:length(allDOM.gr),
                           function(idx){
                             sum(as.matrix(h3k27wt_sqd_h3.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                             # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                           })

stopCluster(cl)
k27domains.grl <- list(raw=NULL,sqdNorm=NULL,h3SqdNorm=NULL)
# Then assign
k27domains.gr <- allDOM.gr
mcols(k27domains.gr)[c(2,3)] <- NULL
k27domains.gr$ID <- 1:length(k27domains.gr)
k27domains.gr$KD_domain_density <- 0
k27domains.gr$WT_domain_density <- 0
k27domains.gr$KD_domain_density <- kdDens
k27domains.gr$WT_domain_density <- wtDens
k27domains.grl[["raw"]] <- k27domains.gr
k27domains.gr$KD_domain_density <- kdDens_sqd
k27domains.gr$WT_domain_density <- wtDens_sqd
k27domains.grl[["sqdNorm"]] <- k27domains.gr
k27domains.gr$KD_domain_density <- kdDens_sqd_h3
k27domains.gr$WT_domain_density <- wtDens_sqd_h3
k27domains.grl[["h3SqdNorm"]] <- k27domains.gr



# MODIF3 -  ASSOCIATE K27 CHIP DENSITY AT BORDER OF EACH DOMAINS FROM 1KB STARTIN AT THE BORDER
h3k27wt.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/WT/bw/chipK27_rawSignal_WT_cuvier_s2_dm3.cov.RData")
h3k27kd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Beaf32_KD/bw/chipK27_rawSignal_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/WT/bw/chipK27_normSeqDepth_WT_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Beaf32_KD/bw/chipK27_normSeqDepth_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/WT/bw/chipK27_LFC_normH3_normSeqDepth_wtxh3_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Beaf32_KD/bw/chipK27_LFC_normH3_normSeqDepth_kdxh3_cuvier_s2_dm3.cov.RData")

myList <- NULL
myList[["raw"]] <- list(wt_cov=h3k27wt.cov,kd_cov=h3k27kd.cov)
myList[["sqdNorm"]] <- list(wt_cov=h3k27wt_sqd.cov,kd_cov=h3k27kd_sqd.cov)
myList[["h3SqdNorm"]] <- list(wt_cov=h3k27wt_sqd_h3.cov,kd_cov=h3k27kd_sqd_h3.cov)
for(range in c(1000,2000)){
  for(type in names(k27domains.grl)){
    kd.cov <- myList[[type]][["kd_cov"]]
    wt.cov <- myList[[type]][["wt_cov"]]
    k27domains.gr <- k27domains.grl[[type]]
    k27KDdensL <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,range,"start"))]),na.rm=T)/range
    k27WTdensL <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,range,"start"))]),na.rm=T)/range
    k27KDdensR <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,range,"end"))]),na.rm=T)/range
    k27WTdensR <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,range,"end"))]),na.rm=T)/range
    k27KDdensM <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,1,"center")+range/2)]),na.rm=T)/range
    k27WTdensM <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,1,"center")+range/2)]),na.rm=T)/range
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensL)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_L_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensL)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_L_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensR)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_R_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensR)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_R_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensM)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_M_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensM)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_M_border_density_",range))
    
    k27domains.grl[[type]] <- k27domains.gr
  }
}

saveRDS(k27domains.grl,paste0(OUT.d,"DATA/PROCESSED/HMM_K27me3_domains_rproc.RData"))

#   -----------------------------------------------------------------------


#  PROCESS NormR H3K27me3 DOMAINS -----------------------------------------------

countConfigSE <- countConfigSingleEnd(binsize=700,shift=73,mapq=30,filteredFlag = 1024)
regime <- regimeR(treatment = trt, control = ctl, genome=gr,models = 3,countConfig = countConfigSE,
                  iterations = 50, procs = 8, verbose = TRUE)

normR.gr <- addSeqinfo(getRanges(regime,fdr=0.1),"dm3","kc",T)
saveRDS(normR.gr,paste0(OUT.d,"/DATA/PROCESSED/h3k27me3_normr_domains.dm3.rds"))
normR.gr <- readRDS(paste0(OUT.d,"/DATA/PROCESSED/h3k27me3_normr_domains.dm3.rds"))
normR.gr <- trim(GenomicRanges::reduce(normR.gr+100))
normR.gr$wdth <- width(normR.gr)



# - HMM
hmmK27_alldom.gr <- readRDS(paste0(OUT.d,"DATA/PROCESSED/HMM_K27me3_domains_rproc.RData"))
hmmK27_alldom.gr <- hmmK27_alldom.gr$raw
hmmK27dom.gr <- subset(hmmK27_alldom.gr,type=="Het")
# - OVERLAP
normR.ovlpHMM.gr <- normR.gr[unique(queryHits(findOverlaps(normR.gr,hmmK27dom.gr)))]
summary(width(hmmK27dom.gr))
# So let's take the same values as the first quartiles of hmm domains (1st Q = 13kb -> we take NORMR domains > 10 kb)
# To avoid cofounding with microDOMAINS and letting enough margin to deal with smaller domains (mDOM)
normR.ovlpHMM.10kbsup.gr <-  subset(normR.ovlpHMM.gr,wdth > 10000)
# OK so let's take normR.gr > 10000 as default K27 domaions
K27_dom.gr <- normR.ovlpHMM.10kbsup.gr
K27_dom.gr <- addSeqinfo(K27_dom.gr)
mcols(K27_dom.gr) <- NULL
# And make noK27 domains
# Use the dm3 granges to make a set diff to obtain noK27 domains (k27.gr + nok27.gr = allGenome.gr)
noK27_dom.gr <- trim(GenomicRanges::setdiff(dm3_mychr.gr,K27_dom.gr))

# Save all
K27_dom.gr$type <- "Het"
noK27_dom.gr$type <- "Euc"
allDOM.gr <- c(K27_dom.gr,noK27_dom.gr)
allDOM.gr <- sortSeqlevels(allDOM.gr);allDOM.gr <- sort(allDOM.gr)



# MODIF2 - ASSOCIATE K27 CHIP DENSITY TO EACH EUC/HET DOMAINS
# LOAD ChipSeq levels (.coverage file)
h3k27wt.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/WT/coverage/k27pm_s2dm3_WT.cov.RData")
h3k27kd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/Beaf32_KD/coverage/k27pm_s2dm3_KD.cov.RData")
h3k27wt_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER//WT/bw/chipK27_normSeqDepth_WT_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/Beaf32_KD/bw/chipK27_normSeqDepth_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/WT/bw/chipK27_LFC_normH3_normSeqDepth_wtxh3_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/H3K27_WT_bKD_CUVIER/Beaf32_KD/bw/chipK27_LFC_normH3_normSeqDepth_kdxh3_cuvier_s2_dm3.cov.RData")
# Load K27 domains (find by normRxHMM and shiftd from 1300bp on right side (bias from HMM see parametrization))
seqlevelsStyle(mrc_tss.gr) <- "Ensembl"
allDOM.gr$k27KDdens <- 0
allDOM.gr$k27WTdens <- 0
# Give each domain a density of K27 level
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, list("allDOM.gr","h3k27kd.cov","h3k27wt.cov","h3k27wt_sqd.cov","h3k27kd_sqd.cov","h3k27wt_sqd_h3.cov","h3k27kd_sqd_h3.cov"),envir = environment())
clusterEvalQ(cl, library(GenomicRanges))
kdDens <- parSapply(cl,1:length(allDOM.gr),
                    function(idx){
                      sum(as.matrix(h3k27kd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                      # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                    })
wtDens <- parSapply(cl,1:length(allDOM.gr),
                    function(idx){
                      sum(as.matrix(h3k27wt.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                      # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                    })
kdDens_sqd <- parSapply(cl,1:length(allDOM.gr),
                        function(idx){
                          sum(as.matrix(h3k27kd_sqd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                          # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                        })
wtDens_sqd <- parSapply(cl,1:length(allDOM.gr),
                        function(idx){
                          sum(as.matrix(h3k27wt_sqd.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                          # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                        })
kdDens_sqd_h3 <- parSapply(cl,1:length(allDOM.gr),
                           function(idx){
                             sum(as.matrix(h3k27kd_sqd_h3.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                             # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                           })

wtDens_sqd_h3 <- parSapply(cl,1:length(allDOM.gr),
                           function(idx){
                             sum(as.matrix(h3k27wt_sqd_h3.cov[allDOM.gr[idx]]),na.rm=T)*10/width(allDOM.gr[idx])
                             # sum(as.matrix(h3k27kd.cov[all_dom.gr[idx]]),na.rm=T)*10/width(all_dom.gr[idx])
                           })



stopCluster(cl)
k27domains.grl <- list(raw=NULL,sqdNorm=NULL,h3SqdNorm=NULL)
# Then assign
k27domains.gr <- allDOM.gr
mcols(k27domains.gr)[c(-1)] <- NULL
k27domains.gr$ID <- 1:length(k27domains.gr)
k27domains.gr$KD_domain_density <- 0
k27domains.gr$WT_domain_density <- 0
k27domains.gr$KD_domain_density <- kdDens
k27domains.gr$WT_domain_density <- wtDens
k27domains.grl[["raw"]] <- k27domains.gr
k27domains.gr$KD_domain_density <- kdDens_sqd
k27domains.gr$WT_domain_density <- wtDens_sqd
k27domains.grl[["sqdNorm"]] <- k27domains.gr
k27domains.gr$KD_domain_density <- kdDens_sqd_h3
k27domains.gr$WT_domain_density <- wtDens_sqd_h3
k27domains.grl[["h3SqdNorm"]] <- k27domains.gr

# MODIF3 -  ASSOCIATE K27 CHIP DENSITY AT BORDER OF EACH DOMAINS FROM 1KB STARTIN AT THE BORDER
myList <- NULL
myList[["raw"]] <- list(wt_cov=h3k27wt.cov,kd_cov=h3k27kd.cov)
myList[["sqdNorm"]] <- list(wt_cov=h3k27wt_sqd.cov,kd_cov=h3k27kd_sqd.cov)
myList[["h3SqdNorm"]] <- list(wt_cov=h3k27wt_sqd_h3.cov,kd_cov=h3k27kd_sqd_h3.cov)
for(range in c(1000,2000)){
  for(type in names(k27domains.grl)){
    kd.cov <- myList[[type]][["kd_cov"]]
    wt.cov <- myList[[type]][["wt_cov"]]
    k27domains.gr <- k27domains.grl[[type]]
    k27KDdensL <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,range,"start"))]),na.rm=T)/range
    k27WTdensL <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,range,"start"))]),na.rm=T)/range
    k27KDdensR <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,range,"end"))]),na.rm=T)/range
    k27WTdensR <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,range,"end"))]),na.rm=T)/range
    k27KDdensM <- rowSums(as.matrix(kd.cov[trim(resize(k27domains.gr,1,"center")+range/2)]),na.rm=T)/range
    k27WTdensM <- rowSums(as.matrix(wt.cov[trim(resize(k27domains.gr,1,"center")+range/2)]),na.rm=T)/range
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensL)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_L_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensL)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_L_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensR)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_R_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensR)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_R_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27KDdensM)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("KD_M_border_density_",range))
    
    mcols(k27domains.gr) <- cbind(mcols(k27domains.gr),k27WTdensM)
    colnames(mcols(k27domains.gr)) <- c(colnames(mcols(k27domains.gr))[-length(colnames(mcols(k27domains.gr)))],paste0("WT_M_border_density_",range))
    
    k27domains.grl[[type]] <- k27domains.gr
  }
}
saveRDS(k27domains.grl,paste0(OUT.d,"/DATA/PROCESSED/h3k27me3_hmmxnormr_domains.rds"))
#   -----------------------------------------------------------------------


# # MDOM ------------------------------------------------------------------


gr <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("dm3")
idx <- which(!gr$circular & gr$SequenceRole=="assembled-molecule")
gr <- gr[idx,1:2]

countConfigSE <- countConfigSingleEnd(binsize=40,mapq=30,filteredFlag = 1024)
k27wtfit <- enrichR(treatment = trt,
                    control   = ctl,
                    genome    = gr,procs=8,countConfig=countConfigSE)

normR.gr <- getRanges(k27wtfit,fdr=0.1)
normR.gr <- addSeqinfo(normR.gr)
normR.gr <- GenomicRanges::reduce(normR.gr+700)-700
saveRDS(normR.gr,paste0(OUT.d,"/DATA/PROCESSED/microdomains_normr.dm3.rds"))

#   -----------------------------------------------------------------------


