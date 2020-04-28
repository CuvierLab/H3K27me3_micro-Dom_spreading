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



# NORMALIZE CONTROL HEATMAP WT by KD TST --------------------------------------
# Compute normFAC
# Preprocessing was erformed as in figure 5A except that we dumped hic interaction between TSSs 

nms <- "int=pcnoBG_anc=gafnoBG"
source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000_CTL_noBG/apa.cfg"))
WT.mtx <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000_CTL_noBG/myMat.rds"))
WT.norm.mtx <- WT.mtx[[nms]]
source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000_CTL_noBG/apa.cfg"))
KD.mtx <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000_CTL_noBG/myMat.rds"))
KD.norm.mtx <- KD.mtx[[nms]]

Y <- do.call(cbind, WT.norm.mtx)
Y_WT <- array(Y, dim=c(dim(WT.norm.mtx[[1]]), length(WT.norm.mtx)))
Y <- do.call(cbind, KD.norm.mtx)
Y_KD <- array(Y, dim=c(dim(KD.norm.mtx[[1]]), length(KD.norm.mtx)))

m_size <- (bin-1)/2 # Middle pixel value - 1
x_msp <- (bin-1)/2+1;y_msp <- (bin-1)/2+1 # Center bin interaction of size (1bin;1bin)
x_mid <- m_size:(m_size+2);y_mid <- m_size:(m_size+2)
y_ll <- 1:3;x_ll <- (bin-2):bin
y_lr <- (bin-2):bin;x_lr <- (bin-2):bin
y_ul <- 1:3;x_ul <- 1:3
y_ur <- (bin-2):bin;x_ur <- 1:3
x_ctrl1 <- ((m_size):(m_size+2));y_ctrl1 <- 1:3
x_ctrl2 <- ((m_size):(m_size+2));y_ctrl2 <- (bin-2):(bin)


UL_WT.mtx <-Y_WT[x_ul,y_ul,]
UL_KD.mtx <- Y_KD[x_ul,y_ul,]
LR_WT.mtx <- Y_WT[x_lr,y_lr,]
LR_KD.mtx <- Y_KD[x_lr,y_lr,]
UL_WT.vec <- apply(UL_WT.mtx,c(3),mean,na.rm=T)
UL_KD.vec <- apply(UL_KD.mtx,c(3),mean,na.rm=T)
LR_WT.vec <- apply(LR_WT.mtx,c(3),mean,na.rm=T)
LR_KD.vec <- apply(LR_KD.mtx,c(3),mean,na.rm=T)

normFAC <- (mean(c(UL_KD.vec,LR_KD.vec),na.rm=T)) / (mean(c(UL_WT.vec,LR_WT.vec),na.rm=T))

#   -----------------------------------------------------------------------


# F6A - PC vs GAF ---------------------------------------------------------------
SUB <- "FIGURES/Figure_6/Figure_6A/"
OUTF.d <- create(OUT.d,SUB)

source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000_CTL_noBG/apa.cfg"))
myMat <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000_CTL_noBG/myMat.rds"))
apa_mat_WT.mtx <- NULL
for(nm in "int=pcnoBG_anc=gafnoBG"){
  MAT <- myMat[[nm]]
  Y <- do.call(cbind, MAT)
  Y <- array(Y, dim=c(dim(MAT[[1]]), length(MAT)))
  # Do the trimmed mean on all the matrix to obtain the APA
  apa_mat_WT.mtx[["sum"]][[nm]] <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  apa_mat_WT.mtx[["mean"]][[nm]] <- apply(Y, c(1, 2), mean, trim=0.025, na.rm = TRUE)
}

source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000_CTL_noBG/apa.cfg"))
myMat <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000_CTL_noBG/myMat.rds"))
apa_mat_KD.mtx <- NULL
for(nm in "int=pcnoBG_anc=gafnoBG"){
  MAT <- myMat[[nm]]
  Y <- do.call(cbind, MAT)
  Y <- array(Y, dim=c(dim(MAT[[1]]), length(MAT)))
  # Do the trimmed mean on all the matrix to obtain the APA
  apa_mat_KD.mtx[["sum"]][[nm]] <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  apa_mat_KD.mtx[["mean"]][[nm]] <- apply(Y, c(1, 2), mean, trim=0.025, na.rm = TRUE)
}


compXpname <- gsub("(.*)_ramirez.*","\\1",xpname)
out_apa.d <- create(paste0(OUTF.d,"/HEATMAP/",compXpname,"_WT_normalized_KD_CTL/"))

my_palette <- colorRampPalette(c("white","black"))(n = 256)
for(type in names(apa_mat_WT.mtx)){
  mat_WT.mtx <- apa_mat_WT.mtx[[type]]
  mat_KD.mtx <- apa_mat_KD.mtx[[type]]
  
  for(nms in names(mat_WT.mtx)){
    od <- create(paste0(out_apa.d,"2D/",nms,"_",type,"/"))
    normMAT.mtx <- mat_WT.mtx[[nms]] * normFAC
    # normMAT.mtx <- apa_mat_KD.mtx[[nms]] * (sum(apa_mat_WT.mtx[[nms]])/sum(apa_mat_KD.mtx[[nms]]))
    
    if(type=="sum"){
      pdf(paste0(od,"heatmap_tadsummits_OE_bin",get("bin"),".pdf"))
      heatmap(as.matrix(normMAT.mtx),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2000,10000))
      heatmap(as.matrix(mat_KD.mtx[[nms]]),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2000,10000))
      dev.off()
    }else{
      pdf(paste0(od,"heatmap_tadsummits_OE_bin",get("bin"),".pdf"))
      heatmap(as.matrix(normMAT.mtx),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2.5,7))
      heatmap(as.matrix(mat_KD.mtx[[nms]]),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2.5,7))
      dev.off()
    }
    
    od <- create(paste0(out_apa.d,"3D/",nms,"/"))
    apa3dPlot(as.matrix(normMAT.mtx),od)
  }
}              

#   -----------------------------------------------------------------------



# F6B - BEAF vs GAF   -----------------------------------------------------------------------
SUB <- "FIGURES/Figure_6/Figure_6B/"
OUTF.d <- create(OUT.d,SUB)

source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000/apa.cfg"))
myMat <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2WT_1000/myMat.rds"))
apa_mat_WT.mtx <- NULL
for(nm in "int=bg_anc=gaf"){
  MAT <- myMat[[nm]]
  Y <- do.call(cbind, MAT)
  Y <- array(Y, dim=c(dim(MAT[[1]]), length(MAT)))
  # Do the trimmed mean on all the matrix to obtain the APA
  apa_mat_WT.mtx[["sum"]][[nm]] <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  apa_mat_WT.mtx[["mean"]][[nm]] <- apply(Y, c(1, 2), mean, trim=0.025, na.rm = TRUE)
}

source(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000/apa.cfg"))
myMat <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/gaf_cp_bg_bs_41_1000_ramirezS2bKD_1000/myMat.rds"))
apa_mat_KD.mtx <- NULL
for(nm in "int=bg_anc=gaf"){
  MAT <- myMat[[nm]]
  Y <- do.call(cbind, MAT)
  Y <- array(Y, dim=c(dim(MAT[[1]]), length(MAT)))
  # Do the trimmed mean on all the matrix to obtain the APA
  apa_mat_KD.mtx[["sum"]][[nm]] <- apply(Y, c(1, 2), sum, na.rm = TRUE)
  apa_mat_KD.mtx[["mean"]][[nm]] <- apply(Y, c(1, 2), mean, trim=0.025, na.rm = TRUE)
}


compXpname <- gsub("(.*)_ramirez.*","\\1",xpname)
out_apa.d <- create(paste0(OUTF.d,"/HEATMAP/",compXpname,"_WT_normalized_KD_CTL/"))

my_palette <- colorRampPalette(c("white","black"))(n = 256)
for(type in names(apa_mat_WT.mtx)){
  mat_WT.mtx <- apa_mat_WT.mtx[[type]]
  mat_KD.mtx <- apa_mat_KD.mtx[[type]]
  
  for(nms in names(mat_WT.mtx)){
    od <- create(paste0(out_apa.d,"2D/",nms,"_",type,"/"))
    normMAT.mtx <- mat_WT.mtx[[nms]] * normFAC
    # normMAT.mtx <- apa_mat_KD.mtx[[nms]] * (sum(apa_mat_WT.mtx[[nms]])/sum(apa_mat_KD.mtx[[nms]]))
    
    if(type=="sum"){
      pdf(paste0(od,"heatmap_tadsummits_OE_bin",get("bin"),".pdf"))
      heatmap(as.matrix(normMAT.mtx),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2000,10000))
      heatmap(as.matrix(mat_KD.mtx[[nms]]),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2000,10000))
      dev.off()
    }else{
      pdf(paste0(od,"heatmap_tadsummits_OE_bin",get("bin"),".pdf"))
      heatmap(as.matrix(normMAT.mtx),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2.5,7))
      heatmap(as.matrix(mat_KD.mtx[[nms]]),Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette,zlim=c(2.5,7))
      dev.off()
    }
    
    od <- create(paste0(out_apa.d,"3D/",nms,"/"))
    apa3dPlot(as.matrix(normMAT.mtx),od)
  }
}              
#   -----------------------------------------------------------------------


# F6C ---------------------------------------------------------------------
SUB <- "FIGURES/Figure_6/Figure_6C/"
OUTF.d <- create(OUT.d,SUB)

source(paste0(OUT.d,"DATA/PROCESSED/APA/bgBordk27_noBg_41_1000_ramirezS2WT_1000/apa.cfg"))
WT.dfp <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/bgBordk27_noBg_41_1000_ramirezS2WT_1000/matrecap_LRC_exctracted.rds"))
source(paste0(OUT.d,"DATA/PROCESSED/APA/bgBordk27_noBg_41_1000_ramirezS2bKD_1000/apa.cfg"))
KD.dfp <- readRDS(paste0(OUT.d,"DATA/PROCESSED/APA/bgBordk27_noBg_41_1000_ramirezS2bKD_1000/matrecap_LRC_exctracted.rds"))


# ASSOCIATE GENES FROM MATRECAP WITH THE CURRENT DATA FRAME
WT.dfp$ancIDmrcTss <- match(WT.dfp$ancID,mrc_tss.gr$gene_id)
KD.dfp$ancIDmrcTss <- match(KD.dfp$ancID,mrc_tss.gr$gene_id)
WT.dfp$intIDmrcTss <- match(WT.dfp$intID,mrc_tss.gr$gene_id)
KD.dfp$intIDmrcTss <- match(KD.dfp$intID,mrc_tss.gr$gene_id)


WT.dfp$msNorm <- WT.dfp$ms / WT.dfp$norm
WT.dfp$mspNorm <- WT.dfp$msp / WT.dfp$norm
WT.dfp$llNorm <- WT.dfp$ll / WT.dfp$norm
WT.dfp$urNorm <- WT.dfp$ur / WT.dfp$norm

KD.dfp$msNorm <- KD.dfp$ms / KD.dfp$norm
KD.dfp$mspNorm <- KD.dfp$msp / KD.dfp$norm
KD.dfp$llNorm <- KD.dfp$ll / KD.dfp$norm
KD.dfp$urNorm <- KD.dfp$ur / KD.dfp$norm

# TEST WT & KD DIST
KD_save.dfp <- KD.dfp
KD_save.dfp$type <- "KD"
WT_save.dfp <- WT.dfp
WT_save.dfp$type <- "WT"
WT_KD.df <- rbind(WT_save.dfp,KD_save.dfp)


idCOLN <- c("intID","ancID", "nme","ancID_net","ancIDmrcTss","intIDmrcTss")
valCOLN <- colnames(WT.dfp)[colnames(WT.dfp) %ni% idCOLN]

for(met in valCOLN){
  WT_KD.df$metrics <- WT_KD.df[,met]
  # print(ggboxplot(WT_KD.df,x="type",y="metrics")+ ylim(c(0,10)) + ggtitle(paste0(toupper(met))))
}



# CREATE DIFFERENTIAL DATA FRAME
# check if wt and kd df are sorted the same way
all.equal(KD.dfp$intID,WT.dfp$intID)
# T
idCOLN <- c("intID","ancID", "nme","ancID_net","ancIDmrcTss","intIDmrcTss")
valCOLN <- colnames(WT.dfp)[colnames(WT.dfp) %ni% idCOLN]

DIVZS.dfp <- WT.dfp[,c("intID","ancID", "nme","ancID_net","ancIDmrcTss","intIDmrcTss")]
for(valC in valCOLN){
  trtSIGNAL <- KD.dfp[,valC]+1e-5
  ctlSIGNAL <- WT.dfp[,valC]+1e-5
  DIVZS.dfp[,valC] <- -(trtSIGNAL - ctlSIGNAL)/sqrt(rowMeans(cbind(trtSIGNAL,ctlSIGNAL)))
}


gafbs <- gaf_kc & bs
gafcp <- gaf_kc & cp
gafcpbs <- gaf_kc & cp & bs

bw.l <- ctl.bwp
pdf(paste0(OUTF.d,"boxplot.pdf"))
for(metrics in paste0(c("mspNorm","urNorm","llNorm"))){
  df_PLOT <- NULL
  
  mydf <- NULL
  for(feat in c("all","gaf_kc","gafcp","gafbs","gafcpbs","cp","coh","ctcf","bs")){
    
    
    # norm_gid <- mrc_tss.gr[unique(subjectHits(myFOL))]$gene_id
    df_pTMP <- DIVZS.dfp %>% filter(ancID_net %in% unique(mrc_tss.gr[get(feat) & nobg]$gene_id)) #Base set is only TSS microDOMAINS
    
    df_pTMP$metrics <- df_pTMP[,metrics]
    df_pTMP$type <- feat
    
    df_PLOT <- rbind(df_PLOT,df_pTMP)
    
    
  }
  mygrp <- list(c("gaf_kc","cp"),c("gaf_kc","coh"),c("gaf_kc","ctcf"),c("gaf_kc","bs"))
  print(ggboxplot(df_PLOT,x="type",y="metrics")+ ggtitle(paste0(toupper(paste0(metrics)))) +
          stat_compare_means(aes(group = type),paired=F,ref.group = "all",method.args = list(alternative="greater"),
                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             label = "..p.signif..") + ylim(c(-10,10))+
          geom_hline(yintercept = 0,col="red")
  )
}
dev.off()



#   -----------------------------------------------------------------------


