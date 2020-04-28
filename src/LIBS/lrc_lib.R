library(GenomicRanges)
library(seqplots)
library(ggplot2)
library(ggpubr)
library("plot3D")
library(plyr)
library(dplyr)
library(reshape2)
library(Matrix)
library("stringr")
library("gridExtra")
require("eulerr")
library(BSgenome.Dmelanogaster.UCSC.dm3)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(scales)
library(trackViewer)
library(RColorBrewer)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(pheatmap)
library(fgsea)

library(foreach)
library(doParallel)
numCores <- detectCores()




`%ni%` <-  Negate('%in%')
# DM 3
kc.chr <- c("2L","2R","3L","3R","X")
s2.chr <- c("2L","2R","3L","3R","X")
sqfDM3 <- SeqinfoForBSGenome(genome="dm3")
sqfDM6 <- SeqinfoForBSGenome(genome="dm6")
dm3.gr <- GRanges(sqfDM3)
dm6.gr <- GRanges(sqfDM6)

# CREATE RANDOM GENOMIC RANGES BASED ON SEQINFO DATA
randomGR <- function(sqINFO=NULL,n=1000) {
  Reduce("c",sample(tileGenome(sqINFO,tilewidth = 200),n,F))
}


MYGENOME <- list(dm3=c(txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene,kc.chr=list(kc.chr),s2.chr=list(s2.chr),sqf=sqfDM3),
                 dm6=c(txdb=TxDb.Dmelanogaster.UCSC.dm6.ensGene,kc.chr=list(kc.chr),s2.chr=list(s2.chr),sqf=sqfDM6)
)


create <- function (...){
  tmp <- list(...)
  dir.create(paste0(...),showWarnings = F,recursive = T)
  
  return(paste0(...))
}



hic_pipeline <- function(l_IN=l_IN,apaCPLS.df=NULL,mBin=20) {
  
  
  HIC.lmtx <- l_IN$HIC
  res <- l_IN$RES
  NM <- l_IN$NM
  
  
  # EXCTRACT HIC MATRICES BINxBIN
  myMAT <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  myMAT <- foreach::foreach(X=1:nrow(apaCPLS.df)) %dopar% {
    binANC <- apaCPLS.df[X,"bin_a_IDX"]
    binINT <- apaCPLS.df[X,"bin_i_IDX"]
    CHR <- apaCPLS.df[X,"chr"]
    if(binINT < binANC){
      TMP <- binINT
      binINT <- binANC
      binANC <- TMP
    }
    # Compute the APA matrix and return it in myMAT as a list of matrix
    feat_apa.mtx <- Matrix(HIC.lmtx[[CHR]][((binANC-mBin):(binANC+mBin)),((binINT-mBin):(binINT+mBin))],sparse=F)
    feat_apa.mtx[!is.finite(feat_apa.mtx) | feat_apa.mtx==0] <- NA
    list(a_IDX=apaCPLS.df[X,"a_IDX"],i_IDX=apaCPLS.df[X,"i_IDX"],MAT=as.matrix(feat_apa.mtx))
  }
  stopImplicitCluster()
  
  
  
  
  # CREATE DATA FRAME OF METRICS
  myMAT.l <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  myMAT.l <- foreach::foreach(X=1:length(myMAT)) %dopar% {
    a_IDX <- myMAT[[X]][["a_IDX"]]
    i_IDX <- myMAT[[X]][["i_IDX"]]
    MAT <- myMAT[[X]][["MAT"]]
    l_ <- list()
    l_[["CP_V"]] <- MAT[i_CP,j_CP]
    l_[["CS_V"]] <- MAT[i_CS,j_CS]
    l_[["ULS_V"]] <- MAT[i_ULS,j_ULS]
    l_[["URS_V"]] <- MAT[i_URS,j_URS]
    l_[["BLS_V"]] <- MAT[i_BLS,j_BLS]
    l_[["BRS_V"]] <- MAT[i_BRS,j_BRS]
    TMP.df <- data.frame(a_IDX=a_IDX,i_IDX=i_IDX)
    for(PAR in names(l_)){
      TMP.df[,PAR] <- mean(as.array(l_[[PAR]]),na.rm=T)
    }
    TMP.df
  }
  stopImplicitCluster()
  METRICS.df <- do.call("rbind",myMAT.l)
  
  
  myMAT_QUANT.l <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  
  
  myMAT_QUANT.l <- foreach::foreach(X=1:length(myMAT)) %dopar% {
    a_IDX <- myMAT[[X]][["a_IDX"]]
    i_IDX <- myMAT[[X]][["i_IDX"]]
    MAT <- myMAT[[X]][["MAT"]]
    MAT_QUANT <- matrix(ntile(MAT,500),APA_SIZE,APA_SIZE)
    l_ <- list()
    l_[["CP_V"]] <- MAT_QUANT[i_CP,j_CP]
    l_[["CS_V"]] <- MAT_QUANT[i_CS,j_CS]
    l_[["ULS_V"]] <- MAT_QUANT[i_ULS,j_ULS]
    l_[["URS_V"]] <- MAT_QUANT[i_URS,j_URS]
    l_[["BLS_V"]] <- MAT_QUANT[i_BLS,j_BLS]
    l_[["BRS_V"]] <- MAT_QUANT[i_BRS,j_BRS]
    TMP.df <- data.frame(a_IDX=a_IDX,i_IDX=i_IDX)
    for(PAR in names(l_)){
      TMP.df[,PAR] <- mean(as.array(l_[[PAR]]),na.rm=T)
    }
    TMP.df
  }
  stopImplicitCluster()
  METRICS_QUANT.df <- do.call("rbind",myMAT_QUANT.l)
  
  
  return(list(METRICS.df=METRICS.df,METRICS_QUANT.df=METRICS_QUANT.df,myMAT=myMAT))
  
  
}



fill_uncovered_domain <- function(mydom.gr,map.grp){
  mappable.gr <- import(map.grp)
  mappable.gr <- addSeqinfo(mappable.gr)
  plotlist <- NULL
  for(nme in unique(mappable.gr$name)){
    assign(paste0(nme,".gr"),mappable.gr[mappable.gr$name==nme])
    plotlist <- c(plotlist,paste0(nme,".gr"))
  }
  xtendUnmp.gr <- GenomicRanges::reduce(`17_Unmap.gr`+ 1000)
  gap.gr <- setdiff(dm3_mychr.gr,mydom.gr)
  mygr_filled.gr <- GenomicRanges::reduce(c(mydom.gr,subsetByOverlaps(gap.gr,xtendUnmp.gr,type="within")))
  return(mygr_filled.gr)
}




fisher_ntiles_plot <- function(F.df,X,Y,ratio=.75) {
  colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))
  
  p <- ggplot(F.df, aes(x = X, y = Y,fill=lfc)) +      geom_tile(col="black") +
    theme_minimal()+
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle = 90,size=25, vjust=1, hjust=.5,face = "bold",
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=20, margin=margin(0,-3,0,0),face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          axis.title.y = element_blank(),
          legend.title.align=0.1,
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    
    geom_text(aes(label =sign),size=4) + 
    scale_fill_gradientn(colours = colfunc(3) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 10,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(-2,0.2,2),
                         labels=c("-1.5","0","1.5"),
                         limits=c(-2,2)) 
  return(p)
}





addSeqinfo <- function(my.gr,g_name="dm3",celltype="kc",keepSeq=T,chrTYPE="Ensembl"){
  sqf <- MYGENOME[[g_name]]$sqf
  chr <- as.character(unlist(MYGENOME[[g_name]][paste0(celltype,".chr")]))
  seqlevelsStyle(my.gr) <- "UCSC"
  seqinfo(my.gr) <- sqf[seqlevels(my.gr),]
  seqlevelsStyle(my.gr) <- "Ensembl"
  if(keepSeq) my.gr <- keepSeqlevels(my.gr,chr,pruning.mode="coarse")
  my.gr <- sortSeqlevels(my.gr)
  my.gr <- sort(my.gr,ignore.strand=T)
  seqlevelsStyle(my.gr) <- chrTYPE
  return(my.gr)
}

surf.colors <- function(x, col = terrain.colors(20)) {
  
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  
  return(colors)
}

apa3dPlot <- function(feat_apa.mtx,o_plot){
  # z <- feat_apa.mtx
  z <- feat_apa.mtx
  jpeg(paste0(o_plot,"_3d_plot.jpeg"),width=5, height=5, units="in", res=500)
  persp(z,col = surf.colors(z,col = rev(heat.colors(80))),theta = 145,phi = 30, border = NA, shade = 2.5)
  dev.off()
  
  jpeg(paste0(o_plot,"_plot3D_plot.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = FALSE)
  
  # arrows3D(11,11,myMat[11,11]+max(myMat)/3,
  #          z1=myMat[11,11]+max(myMat)/8,
  #          col="black",
  #          border="white",
  #          lwd = 5,
  #          add=T,
  #          plot=T,bty = "10",type = "cone")
  dev.off()
  jpeg(paste0(o_plot,"_plot3D_legend_plot.jpeg"),width=5, height=5, units="in", res=500)
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = list(side = 2, length = 0.7))
  
  arrows3D(3,3,myMat[3,3]+max(myMat)/3,
           z1=myMat[3,3]+max(myMat)/8,
           col="black",
           border="white",
           lwd = 5,
           add=T,
           plot=T,bty = "10",type = "cone")
  dev.off()
}


myFisher <- function(a_mat1,b_mat1,a_mat2,b_mat2,meth){
  col1_perc <- (a_mat1*100)/(a_mat1+b_mat1) # 58%
  line1_perc <- (a_mat1*100)/(a_mat1+a_mat2) # 58%
  col2_perc <- (a_mat2*100)/(a_mat2+b_mat2) # 58%
  line2_perc <- (b_mat1*100)/(b_mat1+b_mat2) # 58%
  
  
  
  f_mat <- matrix(c(line_col=a_mat1,line_nocol=a_mat2,noline_col=b_mat1,noline_nocol=b_mat2),nrow=2,ncol=2,byrow=T)
  ft <- fisher.test(f_mat,alternative = meth)
  ft.pv <- ft$p.value # 0.93
  ft.fc <- ft$estimate # 0.97
  return(list(ft.pv=ft.pv,ft.fc=ft.fc,
              line1_perc=line1_perc,
              col1_perc=col1_perc,
              line2_perc=line2_perc,
              col2_perc=col2_perc,
              a_num=a_mat1,b_num=b_mat1,matF=f_mat))
}


mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
}



bordersGR <- function(GR) {
  s.gr <- resize(GR,1,"start")
  e.gr <- resize(GR,1,"end")
  a.gr <- c(e.gr,s.gr)
  a.gr <- sort(a.gr)
  return(list(allBORD=a.gr,leftBORD=s.gr,rightBORD=e.gr))
  
}





fisherGBplot <- function(df,title,out_f,ratio=.15){
  l_mat <- matrix(c(1,2),nrow=2,byrow=T)
  mylay <- layout(mat=l_mat,
                  heights = c(2,10))
  
  # MAT PVAL PLOT ---
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.text.x=element_text(size=10, vjust=1, hjust=.5,
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    geom_text(aes(label = paste0(pv,"\n",odds)),col="black",size=2) +
    
    scale_fill_gradientn(colours=c("firebrick","white","dodgerblue"),
                         breaks=c(-2,0,2),labels=c("-2","0","+2"),
                         limits=c(-2,2))
  fp <- ggarrange(p,nrow=1,heights = c(0.3,0.7))
  ggexport(fp,filename=paste0(out_f,"_MAT_PVAL_BOXPL.pdf"))
  # TIFF PART
  tiff(paste0(out_f,"_MAT_PVAL_PAPER.tiff"),units="in", width=5, height=5, res=300)
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank()) +
    scale_fill_gradient2(low = "firebrick",mid="white", high = "dodgerblue") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  print(p)
  dev.off()
}






seqPlotSDoutliers <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,
                              sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                              ignore.strand=F,leg_pos="topright",
                              COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                              main=""){
  mySMPLseq <- seq(0,1,1e-3)
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  print("totp")
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      # To avoid column entirely covered with 0 (issue with scale and smooth)
      gpsa.mtx <- data.frame(apply(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]]
                                   ,c(1,2),function(X) X + sample(mySMPLseq,1)))
      
      gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      gpsa.scl.mtx[apply(gpsa.scl.mtx,c(1,2),is.nan)] <- NA
      means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC)
  
}







bw_signal <- function(trt.bwp,ctl.bwp,anc.gr) {
  trt.cov <- importData(trt.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)["trtSIGNAL"] <- rowSums(as.matrix(trt.cov[anc.gr]),na.rm = T)/width(anc.gr)
  ctl.cov <- importData(ctl.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)["ctlSIGNAL"] <- rowSums(as.matrix(ctl.cov[anc.gr]),na.rm = T)/width(anc.gr)
  return(anc.gr)
}

bw_signal_simple <- function(ctl.bwp,anc.gr,name="") {
  ctl.cov <- importData(ctl.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)[paste0(name)] <- rowSums(as.matrix(ctl.cov[anc.gr]),na.rm = T)/width(anc.gr)
  return(anc.gr)
}

lfc_zs_dif_decile <- function(anc.gr,NTILE,ASYM=T,nm="") {
  anc.gr$LFC <- log2((anc.gr$trtSIGNAL + 1e-5 ) / (anc.gr$ctlSIGNAL + 1e-5))
  anc.gr$ZS <- (anc.gr$trtSIGNAL -  anc.gr$ctlSIGNAL)/ sqrt(mean(c(anc.gr$trtSIGNAL,anc.gr$ctlSIGNAL)))
  anc.gr$DIF <- (anc.gr$trtSIGNAL -  anc.gr$ctlSIGNAL)
  if(ASYM){
    anc_sup.gr <- subset(anc.gr,LFC >= 0)
    anc_min.gr <- subset(anc.gr,LFC < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decLFC=(NTILE+1)-ntile(LFC,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decLFC=(NTILE*2+1)-ntile(LFC,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
    anc_sup.gr <- subset(anc.gr,ZS >= 0)
    anc_min.gr <- subset(anc.gr,ZS < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decZS=(NTILE+1)-ntile(ZS,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decZS=(NTILE*2+1)-ntile(ZS,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
    anc_sup.gr <- subset(anc.gr,DIF >= 0)
    anc_min.gr <- subset(anc.gr,DIF < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decDIF=(NTILE+1)-ntile(DIF,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decDIF=(NTILE*2+1)-ntile(DIF,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
  }else{
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decLFC=(NTILE+1)-ntile(LFC,NTILE)))
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decZS=(NTILE+1)-ntile(ZS,NTILE)))
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decDIF=(NTILE+1)-ntile(DIF,NTILE)))
  }
  clnm <- colnames(mcols(anc.gr))
  lfc_c <- which(clnm == "LFC");clnm[lfc_c] <- paste0("LFC",nm)
  lfc_c <- which(clnm == "decLFC");clnm[lfc_c] <- paste0("decLFC",nm)
  lfc_c <- which(clnm == "ZS");clnm[lfc_c] <- paste0("ZS",nm)
  lfc_c <- which(clnm == "decZS");clnm[lfc_c] <- paste0("decZS",nm)
  lfc_c <- which(clnm == "DIF");clnm[lfc_c] <- paste0("DIF",nm)
  lfc_c <- which(clnm == "decDIF");clnm[lfc_c] <- paste0("decDIF",nm)
  
  colnames(mcols(anc.gr)) <- clnm
  return(anc.gr)
}



seqPlotFastPERC_LST <- function(struct.ll,tmp,ylim,xlim=xlim,bin=bin,pscale="linear",
                                err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                                ignore.strand=F,leg_pos="topright",
                                COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                                main="",sum_type="mean"){
  
  bw.n <- NULL
  o.tmp <- NULL
  BW.l <- NULL
  trash <- sapply(names(struct.ll),function(X){
    subS.l <- struct.ll[[X]]
    BW <- subS.l[["BW"]]
    GR <- subS.l[["GR"]]
    bw.n[[X]] <<- gsub("(.*).bw","\\1",basename(BW))
    
    sze <- length(get(GR))
    o.tmp <<- unique(c(o.tmp,toString(export.bed(get(GR),paste0(tmp,"/",GR,"_#",sze,"peaks.bed")))))
    BW.l <<- unique(c(BW.l,BW))
  })
  
  gpsa <- getPlotSetArray(BW.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  
  trash <- sapply(names(struct.ll),function(X){
    subS.l <- struct.ll[[X]]
    BW <- subS.l[["BW"]]
    GR <- subS.l[["GR"]]
    bw.n[[X]] <<- gsub("(.*).bw","\\1",basename(BW))
    normF <- subS.l[["NORM"]]
    CENT <- subS.l[["CENTILE"]]
    
    sze <- length(get(GR))
    gpsa.mtx <- data.frame(gpsa.data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["heatmap"]])
    RS <- data.frame(RS=rowSums(gpsa.mtx,na.rm=T))
    RS <- RS %>% mutate(CENT=ntile(RS,100))
    ROW_KEEP <- RS[,"CENT"] %in% CENT
    
    gpsa.mtx <- gpsa.mtx[ROW_KEEP,] 
    if(sum_type=="mean"){
      if(normF[["BOOL"]]==T){
        if(normF[["TYPE"]]=="exp"){
          means <- exp(colMeans(gpsa.mtx,na.rm=T))*normF[["LINFAC"]] # Now you can do the mean on original data without 3 SD away outliers
        }else if(normF[["TYPE"]]=="pow"){
          means <- `^`(colMeans(gpsa.mtx,na.rm=T),normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- log(colMeans(gpsa.mtx,na.rm=T),base = normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="MIN"){
          means <- colMeans(gpsa.mtx,na.rm=T)-normF[["LINFAC"]]
        }else{
          means <- colMeans(gpsa.mtx,na.rm=T)*normF[["LINFAC"]]
        }
      }else{
        means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      }
    }else{
      mat_c <- gpsa.mtx 
      
      if(normF[["BOOL"]]==T){
        if(normF[["TYPE"]]=="exp"){
          means <- exp(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }))*normF[["LINFAC"]] # Now you can do the mean on original data without 3 SD away outliers
        }else if(normF[["TYPE"]]=="pow"){
          means <- `^`(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }),normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- log(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }),base = normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- apply(mat_c,2,function(X){
            median(X,na.rm=T)
          })-normF[["LINFAC"]]
        }else{
          means <- apply(mat_c,2,function(X){
            median(X,na.rm=T)
          })*normF[["LINFAC"]]
        }
      }else{
        means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      }
    }
    print(means)
    if(smooth){
      means = smooth.spline(1:(length(means)), means, spar=spar)$y
    }
    stderror <- apply(gpsa.mtx,2,function(n){
      sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
    })
    conint <- apply(gpsa.mtx , 2, function(n) {
      qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
    })
    stderror[is.na(stderror)] <- 0
    conint[is.na(conint)] <- 0
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["means"]] <<- means # change the means vector from getPlotSetArray object
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["stderror"]] <<- stderror # change the means vector from getPlotSetArray object
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["conint"]] <<- conint  # change the means vector from getPlotSetArray object
    
  })
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              plotScale = pscale,
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC) 
}


createGRcouplesAPA <- function(apa_size=41,hic_res=1000,anc.gr=NULL,int.gr=NULL,constraint.gr=NULL,SQF=sqfDM3,MYCHR=dm3CHR){
  
  if(is.null(constraint.gr)){
    print("NO DOMAIN TO RESTRICT INTERACTIONS, USE ANOTHER FUNCTION OR PROVIDE DOMAIN GRANGES")
    return(F)
  }
  
  # DELETE MCOLS TO AVOID FIRST.X.BLAHBLAHBLAH
  mcols(anc.gr) <- NULL
  mcols(int.gr) <- NULL
  mcols(constraint.gr) <- NULL
  names(anc.gr) <- NULL
  names(int.gr) <- NULL
  names(constraint.gr) <- NULL
  seqlevelsStyle(anc.gr) <- "Ensembl"
  seqlevelsStyle(int.gr) <- "Ensembl"
  seqlevelsStyle(constraint.gr) <- "Ensembl"
  
  # Bin size 
  mBin <- (apa_size-1)/2
  #Interacting feature
  #Interacting feature
  feati.gr <- resize(int.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feati.gr$bin_i_IDX <- ceiling((start(feati.gr))/hic_res)
  feati.gr$i_IDX <- 1:length(feati.gr)
  dfi <- data.frame(findOverlapPairs(constraint.gr,feati.gr))
  dfi$filt <- paste0(dfi$first.seqnames,"_",dfi$first.start)
  dfi <- dfi[,c("second.bin_i_IDX","second.i_IDX","filt")]
  colnames(dfi) <- c("bin_i_IDX","i_IDX","cons_info")
  
  #Anchoring feature
  feata.gr <- resize(anc.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feata.gr$bin_a_IDX <- ceiling((start(feata.gr))/hic_res)
  feata.gr$a_IDX <- 1:length(feata.gr)
  dfa <- data.frame(findOverlapPairs(constraint.gr,feata.gr))
  dfa$filt <- paste0(dfa$first.seqnames,"_",dfa$first.start)
  dfa <- dfa[,c("second.bin_a_IDX","second.a_IDX","filt")]
  colnames(dfa) <- c("bin_a_IDX","a_IDX","cons_info")
  
  # Merge dataframes
  mia <- merge(dfi,dfa,by = "cons_info")
  # Remove couples that are < 15kb/binS fragment away
  if(length(which(abs(mia$bin_i_IDX-mia$bin_a_IDX)<(15000/hic_res))))
    mia <- mia[-which(abs(mia$bin_i_IDX-mia$bin_a_IDX)<(15000/hic_res)),]
  
  # Add distance between int and anc in bin 
  mia$dst_bin <- abs(mia$bin_i_IDX-mia$bin_a_IDX)
  # Add chromosome info
  mia$chr <- gsub("(.*)_.*","\\1",mia$cons_info)
  
  
  
  # remove peaks that fall less than bin+1 away from start and end on each chromosome
  MYDF <- NULL
  for(CHR in MYCHR){
    print(CHR)
    miaTMP <- subset(mia,chr==CHR)
    DIM <- as.numeric(ceiling(seqlengths(SQF[paste0("chr",CHR)])/hic_res))
    toDEL <- which(miaTMP$bin_i_IDX >= DIM-mBin | miaTMP$bin_i_IDX <= mBin | miaTMP$bin_a_IDX >= DIM-mBin | miaTMP$bin_a_IDX <= mBin  )
    if(length(toDEL)) miaTMP <- miaTMP[-toDEL,]
    MYDF <- rbind(MYDF,miaTMP)
  }
  return(MYDF)
}





fisherDistPlotBothSided <- function(anchor,bait,inter,maxrange,step,window){
  # Output setup
  main <- paste0("Anchor=",anchor,"_Bait=",bait,"_Intersect=",inter,"_maxrange=",maxrange,"KB_step=",step,"_Win=",window)
  wdth <- 25000/(maxrange/step)
  #Get the data : from character to vector
  anchor <- get(anchor);bait <- get(bait);inter <- get(inter)
  # Turn TSS anchor point to Granges vector
  anchor.gr <- mrc_tss.gr[anchor]
  # Create a reference set of matrecap GRanges index that fall into the range of analysis
  refset_idx <- unique(subjectHits(findOverlaps(trim(anchor.gr+(maxrange+window)),mrc_tss.gr)))
  # get filter set
  bait_idx <- which(bait==1)
  # Apply refset and filterset
  tokeep <- intersect(refset_idx,bait_idx)
  mrc_filt_refeset.df <- matrecap[tokeep,]
  mrc_filt_refset.gr <- mrc_tss.gr[tokeep]
  refset_size <- length(mrc_filt_refset.gr)
  inter_filt_refset.gid <- mrc_filt_refset.gr[which(inter[tokeep]==1)]$gene_id
  inter_size <- length(inter_filt_refset.gid)
  shft <- seq(0,maxrange,step)
  window.gr <- resize(anchor.gr,window,"start",ignore.strand=T)
  myvec <- sapply(shft,
                  function(sshft){
                    sbstp.gr <- subsetByOverlaps(mrc_filt_refset.gr,trim(shift(window.gr,sshft)))
                    sbstn.gr <- subsetByOverlaps(mrc_filt_refset.gr,trim(shift(window.gr,-(sshft+window)))) # Added
                    sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
                    if(length(sbst.gr)){
                      win_set <- sbst.gr$gene_id
                      a1 <- length(which(win_set %in% inter_filt_refset.gid))
                      a2 <- length(win_set) - a1
                      b1 <- inter_size - a1
                      b2 <- refset_size - (inter_size + length(win_set) - a1)
                      # print(matrix(c(a1,a2,b1,b2),2,2))
                      pv <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")$p.value
                      perc <- 100*(a1/(a1+a2))
                      return(cbind(pv,perc))
                    }
                    else
                      return(cbind(1,NA))
                  })
  perc.df <- NULL
  myvec <- t(myvec)
  
  trsh <- sapply(1:(nrow(myvec)-3),
                 function(idx){
                   # print(idx)
                   perc.df <<- rbind.data.frame(perc.df,cbind.data.frame(x=mean(myvec[c(idx:(idx+3)),2]),sd=sd(myvec[c(idx:(idx+3)),2])))
                 })
  
  perc.df <- cbind(perc.df,shft[-((length(shft)-2):length(shft))]+window/2)
  colnames(perc.df) <- c("mean","sd","x")
  return(perc.df)
}

