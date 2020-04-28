
#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          PATH  ----
#####################################################################################-
SPATH <- dirname(rstudioapi::getSourceEditorContext())
OUT.d <- create(paste0(SPATH,"/"))

#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source("../LIBS/lrc_lib.R")


#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-
# USEFULL FEATURES (TSS, MATRECAP)
tss.gr <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),1,"start")
Dmt=transcriptsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene,by="gene")
matrecap <- read.csv(paste0(OUT.d,"/DATA/PUBLISHED/gene_mat_recap.csv"))

# Genes class
# 1 INSIDE K27me3 domains
# 2 OUTSIDE K27me3 domains
# 3 > 2.5 Kb Pointing Upstream K27me3 domains
# 4 < 2.5 Kb Pointing Downstram K27me3 domains
k27_class <- readRDS(paste0(OUT.d,"/DATA/PUBLISHED/K27_genes_class.rds"))



#####################################################################################-
#
#          RUN
#
#####################################################################################-
# First, filter genes from matrecap that are not in the list of genes from TxDb ----
tss.gr <- sortSeqlevels(tss.gr)
sum(matrecap$Fbgnid %ni% names(tss.gr))
# we remove 1279 genes that are not in UCSC fbgnid
# It seems that they are only peaks (motif for instance)
tss.gr <- sort(tss.gr,ignore.strand=T)
g_fbgn <- as.character(matrecap[,"Fbgnid"])
mtch <- match(g_fbgn,names(tss.gr))
mtch <- mtch[!is.na(mtch)]
tssMTC.gr <- tss.gr[mtch]
mtch2 <- match(names(tssMTC.gr),g_fbgn)
mtch2 <- mtch2[!is.na(mtch2)]
matrecap <- matrecap[mtch2,]
tssMTC.gr$class <- matrecap[mtch2,]$k27Class

# ADD K27me3 WT vs KD Zscore 
zscoreM500 <- readRDS(paste0(OUT.d,"DATA/PROCESSED/ZSCORE/K27_wt_kd_zscoreM500_matrecap.RData"))
zscorePM1kb <- readRDS(paste0(OUT.d,"DATA/PROCESSED/ZSCORE/K27_wt_kd_zscorePM1KB_matrecap.RData"))
zscorePM2kb <- readRDS(paste0(OUT.d,"DATA/PROCESSED/ZSCORE/K27_wt_kd_zscorePM2KB_matrecap.RData"))
sum(matrecap$gene %ni% zscoreM500$genes) # same goes for PM 1,2Kb
# We remove 475 genes here
# So in the end we have 14753 - 1279-453 = 12999 genes
myzs <- c("zscoreM500","zscorePM1kb","zscorePM2kb")
for(zs in myzs){
  
  zs.new <- get(zs)
  sze <- gsub("zscore(.*)","\\1",zs)
  colnames(zs.new) <- c("gene",zs,paste0("pv",sze),paste0("pv_adj",sze))
  matrecap <- merge(matrecap,zs.new,by.x="gene",by.y="gene")
  matrecap[,paste0("upK27b",sze)] <- 0
  matrecap[,paste0("dwnK27b",sze)] <- 0
  
  # Function to compute Ntile
  tomuteDec <- paste0('ntile(',paste0(zs),',10)')
  tomuteQuint <- paste0('ntile(',paste0(zs),',5)')
  
  matrecap[which((matrecap[,zs] > 0) & (matrecap[,paste0("pv",sze)] < 0.05)),paste0("upK27b",sze)] <- 1 
  matrecap[which((matrecap[,zs] < 0) & (matrecap[,paste0("pv",sze)] < 0.05)),paste0("dwnK27b",sze)] <- 1 
  matrecap <- matrecap %>% mutate_(.dots=setNames(tomuteDec,paste0(zs,"_decile"))) %>% mutate_(.dots=setNames(tomuteQuint,paste0(zs,"_quintile"))) 
}


# Add K27 Class (In, out, etc...) ----
matrecap$k27Class <- 0
sapply(names(k27_class),
       function(nms){
         class <- gsub(".*class(.)","\\1",nms)
         print(class)
         print(matrecap[as.character(matrecap$gene) %in% k27_class[[nms]],]$k27Class)
         matrecap[which(as.character(matrecap$gene) %in% k27_class[[nms]]),]$k27Class <<- class
       })
saveRDS(matrecap,paste0(OUT.d,"DATA/PROCESSED/matrecap.dm3.rds"))

#   -----------------------------------------------------------------------


















# 5 - Add BEAF-32  modencode Data (Karpen Gary, ID 274 922 3745) ----
matrecap <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_211217_zscore.rds")
toto <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/4_Analysis/Fisher/raw/l_BS_GRange.rds")
g_Beaf_bs <- toto$`BEAF-32`
g_bbs_cgata <- g_Beaf_bs[which(g_Beaf_bs %in% matrecap[which(matrecap$cgata==1),"Fbgnid"])]
matrecap$beafg_s2_modencode <- 0
matrecap[matrecap$Fbgnid %in% g_bbs_cgata,]$beafg_s2_modencode <- 1
sum(matrecap$beafg_s2_modencode) # 536 tss with cgata & beaf from karpen
sum(matrecap$BEAFg) # 3810 tss with cgata & beaf from homemade detection (hand annotated > 10 reads, not with MACS2)
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/raw/Matrecap/mrc_120118_beaf_modencode.rds")

#   -----------------------------------------------------------------------


# 6 Add cohesin KC167 ----
matrecap <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/raw/Matrecap/mrc_120118_beaf_modencode.rds")
tss.gr <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),1,"start")
mtch <- match(matrecap$Fbgnid,tss.gr$gene_id)
# Don't sort so we can use binary way to retrieve GRange of the corresponding TSS
mrc_tss.gr <- tss.gr[mtch]
coh.gr <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/raw/Data_resource/RDS/All_rds/Rad21_macs2_Kc167_dm3.gr.RData")
seqlevelsStyle(coh.gr) <- "UCSC"
coh_idx <- unique(queryHits(findOverlaps(mrc_tss.gr,coh.gr,maxgap = 250)))
matrecap$coh_kc167 <- 0 
matrecap[coh_idx,"coh_kc167"] <- 1 
sum(matrecap$coh_kc167) # 5197
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/raw/Matrecap/mrc_190118_rad21_kc167.rds")

#   -----------------------------------------------------------------------


# 7 Adding KC167 TAD border -----------------------------------------------
matrecap <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_190118_rad21_kc167.rds")
tss.gr <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),1,"start")
mtch <- match(matrecap$Fbgnid,tss.gr$gene_id)
mrc_tss.gr <- tss.gr[mtch]


# Filters ----
sexton_tad <- "/media/alexandre/Data/Recherche/LBME/Database/Drosophila/Kc167/dm3/tads/Corces/allTAD_corces_kc167_dm3.GR.RData"
ramirez_tad <- "/media/alexandre/Data/Recherche/LBME/Database/Drosophila/Kc167/dm3/tads/Ramirez/HiCexplorer/Ramirez_TADs_kc167_dm3.GR.RData"
wang_tad <- "/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2R/TAD/wang_dm3_s2r.GR.RData" 
xtend <- c(250,500,1000,2000,4000)
sxtad.gr <- readRDS(sexton_tad);seqlevelsStyle(sxtad.gr) <- "UCSC"
ramtad.gr <- readRDS(ramirez_tad);seqlevelsStyle(ramtad.gr) <- "UCSC"
wangtad.gr <- readRDS(wang_tad);seqlevelsStyle(wangtad.gr) <- "UCSC"

for(i in xtend){
  sxtadbor.gr <- trim(c(resize(sxtad.gr,1,"start")+i,resize(sxtad.gr,1,"end")+i))
  ramtadbor.gr <- trim(c(resize(ramtad.gr,1,"start")+i,resize(ramtad.gr,1,"end")+i))
  wangtadbor.gr <- trim(c(resize(wangtad.gr,1,"start")+i,resize(wangtad.gr,1,"end")+i))
  sxt_vec <- data.frame(sext=rep(0,length(mrc_tss.gr)))
  ramt_vec <- data.frame(ram=rep(0,length(mrc_tss.gr)))
  wang_vec <- data.frame(wang=rep(0,length(mrc_tss.gr)))
  
  idx_sx_1 <- unique(queryHits(findOverlaps(mrc_tss.gr,sxtadbor.gr)))
  idx_ram_1 <- unique(queryHits(findOverlaps(mrc_tss.gr,ramtadbor.gr)))
  idx_wang_1 <- unique(queryHits(findOverlaps(mrc_tss.gr,wangtadbor.gr)))
  
  sxt_vec[idx_sx_1,1] <- 1
  ramt_vec[idx_ram_1,1] <- 1
  wang_vec[idx_wang_1,1] <- 1
  
  colnames(sxt_vec) <- paste0("sexton_TADbord_Kc167_pm",i)
  colnames(ramt_vec) <- paste0("ramirez_TADbord_Kc167_pm",i)
  colnames(wang_vec) <- paste0("wang_TADbord_S2_pm",i)
  
  matrecap <- cbind(matrecap,sxt_vec,ramt_vec,wang_vec)
}

file.remove(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_290118_TAd_kc167.rds")
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/mrc_290118_TAd_kc167.rds")
#   -----------------------------------------------------------------------



# 8 - adding Newly treated RNAseq from node2 with Deseq2 and STAR
# matrecap <- readRDS("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_290118_TAd_kc167.rds")
# dif.tbl <- read.csv("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/RNAseq/Team/Beaf32/Lhoumaud_paper/STAR_aln/DEseq/Gmax_DESeq2-results-with-normalized.csv")
# 
# newdif.tbl <- dif.tbl[which(dif.tbl$gene %in% matrecap$Fbgnid),]
# FbgnP <- newdif.tbl[which(newdif.tbl$log2FoldChange <= 0),"gene"]
# FbgnM <- newdif.tbl[which(newdif.tbl$log2FoldChange > 0),"gene"]
# 
# matrecap$beafMalex <- 0
# matrecap$beafPalex <- 0
# 
# matrecap[match(FbgnM,matrecap$Fbgnid),"beafMalex"] <- 1
# matrecap[match(FbgnP,matrecap$Fbgnid),"beafPalex"] <- 1
# file.remove(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
# saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_190418_TAd_kc167.rds")
# saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/mrc_190418_TAd_kc167.rds")
# 
# 
# sum(matrecap$beafPalex)



# 9 - Add new K27 domain border and adjust depending parameters -----------
matrecap <- readRDS(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
allDOM.gr <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Domains/K27_PM/k27_noK27_domains_s2dm3_Chipdensity.gr.RData")
het_dom.gr <- allDOM.gr[allDOM.gr$type=="Het"]
euc_dom.gr <- allDOM.gr[allDOM.gr$type=="Euc"]

# Now recreate the 4 Classes of genes depending on their overlap with border, inside the domains, outside.
# Class1 : Tss have to be inside the domain -1.5 kb from L/R borders
tss_c1.gr <- resize(het_dom.gr,1,"center") + (width(het_dom.gr)-1500*2)/2
# Class2 : The opposite of class1, this time it's inside euchromatin domains
tss_c2.gr <- resize(euc_dom.gr,1,"center") + (width(euc_dom.gr)-1500*2)/2

which(width(euc_dom.gr) < 3000)

file.remove(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_290118_TAd_kc167.rds")
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/mrc_290118_TAd_kc167.rds")

#   -----------------------------------------------------------------------



# 9 - Add new K27 domain border and adjust depending parameters -----------
matrecap <- readRDS(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))

toKeep <- which(as.vector(seqnames(mrc_tss.gr)) %in% s2.chr)
matrecap <- matrecap[toKeep,]
file.remove(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/mrc_010618_s2chr.rds")
saveRDS(matrecap,"/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/mrc_010618_s2chr.rds")

#   -----------------------------------------------------------------------



# ADD ZSCORE COMPUTED WITH H3normSQD & normSQD K27 ChIPseq ----------------
matrecap <- readRDS(dir("/media/alexandre/Data/Recherche/LBME/Projects/Beaf32_LRC/exp/Proper_pipeline/00_Data_processing/matrecap_cuvier_s2/results/Matrecap/current/",full.names = T))
tss.gr <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),1,"start")
seqlevelsStyle(tss.gr) <- "Ensembl"
mtch <- match(matrecap$Fbgnid,tss.gr$gene_id)
mrc_tss.gr <- tss.gr[mtch]

h3k27wt_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/WT/bw/chipK27_normSeqDepth_WT_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Beaf32_KD/bw/chipK27_normSeqDepth_KD_cuvier_s2_dm3.cov.RData")
h3k27wt_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/WT/bw/chipK27_LFC_normH3_normSeqDepth_wtxh3_cuvier_s2_dm3.cov.RData")
h3k27kd_sqd_h3.cov <- readRDS("/media/alexandre/Data/Recherche/LBME/Database/Drosophila/S2/dm3/histone_marks/h3k27/h3k27_wt_beaf32-kd_mutant_cuvier/Beaf32_KD/bw/chipK27_LFC_normH3_normSeqDepth_kdxh3_cuvier_s2_dm3.cov.RData")

myList <- NULL
myList[["sqdNorm"]] <- list(wt_cov=h3k27wt_sqd.cov,kd_cov=h3k27kd_sqd.cov)
myList[["h3SqdNorm"]] <- list(wt_cov=h3k27wt_sqd_h3.cov,kd_cov=h3k27kd_sqd_h3.cov)
for(type in c("sqdNorm","h3SqdNorm")){
  kd.cov <- myList[[type]][["kd_cov"]]
  wt.cov <- myList[[type]][["wt_cov"]]
  mrc_tss_m500.gr <- extend(mrc_tss.gr,500,0)
  mrc_tss_pm1kb.gr <- extend(mrc_tss.gr,1000,1000)
  mrc_tss_pm2kb.gr <- extend(mrc_tss.gr,2000,2000)
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

#   -----------------------------------------------------------------------


