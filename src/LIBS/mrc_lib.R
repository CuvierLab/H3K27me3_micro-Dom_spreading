
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



#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-


# EXCTRACT DATA FROM MATRECAP -------------------------------------------

matrecap <- readRDS(paste0(OUT.d,"DATA/PROCESSED/matrecap.dm3.rds"))
tss.gr <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene),1,"start")
seqlevelsStyle(tss.gr) <- "Ensembl"
mtch <- match(matrecap$Fbgnid,tss.gr$gene_id)
mrc_tss.gr <- tss.gr[mtch]

# Up/Down K27 Chipseq signal
for(var in c("M500","PM1kb","PM2kb")){
  assign(paste0("upK27b",var),  matrecap[,paste0("upK27b",var)]==1)
  assign(paste0("dwnK27b",var),  matrecap[,paste0("dwnK27b",var)]==1)
  assign(paste0("equK27b",var) ,!get(paste0("dwnK27b",var)) & !get(paste0("dwnK27b",var)))
  assign(paste0("zscore",var),  matrecap[,paste0("zscore",var)])
}

# Top X%
for(var in c("M500","PM1kb","PM2kb")){
  mrc_zs_bord <- cbind(zs=matrecap[,paste0("zscore",var)] ,idx=1:nrow(matrecap[,]))
  mrc_abs_zs_bord <- cbind(zs=abs(matrecap[,paste0("zscore",var)]) ,idx=1:nrow(matrecap[,]))
  
  dwn2up <- mrc_zs_bord[order(mrc_zs_bord[,"zs"],decreasing = F),]
  mid2max <- mrc_abs_zs_bord[order(mrc_abs_zs_bord[,"zs"],decreasing = F),]
  for(p in c(0.01,.02,.05,.1)){
    # TOP DOWNK27
    top <- rep(0,nrow(matrecap))
    topXp_dwn_idx <- dwn2up[ dwn2up[,"zs"] < quantile(dwn2up[,"zs"],p,type=5), "idx"]
    top[topXp_dwn_idx] <- 1
    assign(paste0("top_",p*100,"p_downK27",var),top)
    # TOP UPK27
    top <- rep(0,nrow(matrecap))
    topXp_up_idx <- dwn2up[ dwn2up[,"zs"] > quantile(dwn2up[,"zs"],1-p,type=5), "idx"]
    top[topXp_up_idx] <- 1
    assign(paste0("top_",p*100,"p_upK27",var),top)
    # TOP MIDK27
    top <- rep(0,nrow(matrecap))
    topXp_mid_idx <- mid2max[ mid2max[,"zs"] > quantile(mid2max[,"zs"],1-p,type=5), "idx"]
    top[topXp_mid_idx] <- 1
    assign(paste0("top_",p*100,"p_midK27",var),top)
  }
}

# Create basic feature of border K27 (4 classes: in, out, border with gene downstream the domain, upstream the domain)
ink <- matrecap$k27Class==1
outk <- matrecap$k27Class==2
upsk <- matrecap$k27Class==3
dnsk <- matrecap$k27Class==4
bordK27 <- matrecap$k27Class==4 | matrecap$k27Class==3
noborderK27 <- !bordK27


upK27bM500 <- matrecap$upK27bM500
dwnK27bM500 <- matrecap$dwnK27bM500

# Generate variable for all available RNAseq in beafKD
dr_pric <- matrecap$BeafmPric==1
ur_pric <- matrecap$BeafpPric==1
dr_mag <- matrecap$BeafmMag==1
ur_mag <- matrecap$BeafpMag==1
dr_paul <- matrecap$BeafmPaul==1
ur_paul <- matrecap$BeafpPaul==1
dr_rna <- matrecap$RNAm==1
ur_rna <- matrecap$RNAp==1
urdr_l <- list(pric=cbind(dr=dr_pric,ur=ur_pric),
               mag=cbind(dr=dr_mag,ur=ur_mag),
               paul=cbind(dr=dr_paul,ur=ur_paul),
               rna=cbind(dr=dr_rna,ur=ur_rna))


# Create protein binding related variables
# h for homemade beaf, m for modencode beaf
bg <- matrecap$BEAFg==1
bg_m <- matrecap$beafg_s2_modencode==1
bs<- matrecap$BEAFs==1
coh<- matrecap$cohesin==1
ctcf<- matrecap$ctcfg==1
dref<- matrecap$drefg==1
cp<- matrecap$CP190==1
gaf<- matrecap$gafg==1

matrecap$polycomb <- 0
idx <- unique(queryHits(findOverlaps(mrc_tss.gr,pc.dm3.gr,maxgap=500)))
matrecap[idx,"polycomb"] <- 1 
polycomb <- matrecap[,"polycomb"] == 1

all<- rep(T,nrow(matrecap))

nobg <- !bg
noctcf <- !ctcf
nocp <- !cp
nogaf <- !gaf
nodref <- !dref
nocoh <- !coh
noall <- !all


gaf_kc.gr <- addSeqinfo(readRDS(paste0(OUT.d,"DATA/PUBLISHED/CHIPSEQ/Gaf_kc167_dm3.GR.RData")),g_name = "dm3","kc")
matrecap$gaf_kc <- 0
idx <- unique(queryHits(findOverlaps(mrc_tss.gr,gaf_kc.gr,maxgap=500)))
matrecap[idx,"gaf_kc"] <- 1 
gaf_kc <- matrecap[,"gaf_kc"] == 1
nogaf_kc <- matrecap[,"gaf_kc"] == 0


# SuHW
suhw.dm3.gr <- addSeqinfo(import(paste0(OUT.d,"DATA/PUBLISHED/CHIPSEQ/GSE41950_Suhw_S2.bed.gz")),g_name = "dm3","kc")
matrecap$suhw <- 0
idx <- unique(queryHits(findOverlaps(mrc_tss.gr,suhw.dm3.gr,maxgap=500)))
matrecap[idx,"suhw"] <- 1 
suhw <- matrecap[,"suhw"] == 1
nosuhw <- matrecap[,"suhw"] == 0


# M1BP
m1bp.dm3.gr <- addSeqinfo(import(paste0(OUT.d,"DATA/PUBLISHED/CHIPSEQ/GSE101554_S2_M1BP_ChIP_peaks.bed.gz")),g_name = "dm3","kc")
matrecap$m1bp <- 0
idx <- unique(queryHits(findOverlaps(mrc_tss.gr,m1bp.dm3.gr,maxgap=500)))
matrecap[idx,"m1bp"] <- 1 
m1bp <- matrecap[,"m1bp"] == 1
nom1bp <- matrecap[,"m1bp"] == 0




gafnob <- gaf & !bg
bgbord <- bg & bordK27
nobgbord <- !bg & bordK27
nobgbs <- !bg & bs
nobggafctcfcpbs <- !bg & (gaf | ctcf | cp | bs)



