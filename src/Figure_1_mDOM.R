#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          PATH  ----
#####################################################################################-
SPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
OUT.d <- (paste0(SPATH,"/"))

#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

source(paste0(OUT.d,"/LIBS/lrc_lib.R")) # Global lib
source(paste0(OUT.d,"/LIBS/data_lib.R")) # load data
source(paste0(OUT.d,"/LIBS/mrc_lib.R")) # matrecap lib

#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-  

#####################################################################################-
#  FIGURE 1B 1C -----------------------
# Plot profile K27me3 at beaf border
# OPTS
SUB <- "FIGURES/Figure_1/Figure_1B_1C"
SD <- 2
mxgap <- 1000
norm <- "NORMH3"

OUTF.d <- create(OUT.d,SUB,"/PROFILE_K27_WT_KD_SIG_AT_TSS_BEAF_NORMR_K27DOMAINS/NORM_TYPE_",norm,"/MAXGAP_",mxgap,"/")
HET.L.gr <- resize(Het.gr,1,"start")
HET.R.gr <- resize(Het.gr,1,"end")

# Test
bg.borderL.gr <- subsetByOverlaps(mrc_tss.gr[bg],HET.L.gr,maxgap = mxgap)
bg.borderR.gr <- subsetByOverlaps(mrc_tss.gr[bg],HET.R.gr,maxgap = mxgap)
strand(bg.borderL.gr) <- "-"
strand(bg.borderR.gr) <- "+"
bg.border.gr <- unique(c(bg.borderL.gr,bg.borderR.gr))

# Ctrl
nobg.borderL.gr <- subsetByOverlaps(mrc_tss.gr[!bg & !cp & !bs & !gaf_kc],HET.L.gr,maxgap = mxgap)
nobg.borderR.gr <- subsetByOverlaps(mrc_tss.gr[!bg & !cp & !bs & !gaf_kc],HET.R.gr,maxgap = mxgap)
strand(nobg.borderL.gr) <- "-"
strand(nobg.borderR.gr) <- "+"
nobg.border.gr <- unique(c(nobg.borderL.gr,nobg.borderR.gr))

pdf(paste0(OUTF.d,"Profile_K27_WT_vs_KDBeaf_at_Beaf_BORDERS.pdf"))
bw.l <- bw.ll[[norm]][["sig"]]
YLIM <- bw.ll[[norm]][["ylim"]]
XLIM <- c(10000,10000)
seqPlotSDoutliers(bw.l,TMP.d,c("bg.border.gr"),ylim=YLIM,xlim=XLIM,bin=50L,sd=SD,err = T,smooth = T,spar = 0.2,type="pf",gnme = "dm3",ignore.strand = F,COLVEC = c("dodgerblue","firebrick","black"))
seqPlotSDoutliers(bw.l,TMP.d,c("nobg.border.gr"),ylim=YLIM,xlim=XLIM,bin=50L,sd=SD,err = T,smooth = T,spar=0.2,type="pf",gnme = "dm3",ignore.strand = F,COLVEC = c("dodgerblue","firebrick","black"))
dev.off()

#   -----------------------------------------------------------------------


#  FIGURE 1D -----------------------
#   UP K27 FOR ALL TSS DISTANCE PLOT - BOTHSIDED
SUB <- "FIGURES/Figure_1/Figure_1D"
MAX <- 60e3
SHFT <- 2e3
WINSIZE <- 1e3
DIF_K27_BIN <- "M500"
OUTF.d <- create(OUT.d,SUB,"/DISTRIBUTION_DISTANCE_BASED/BOTH_SIDED/MAXRANGE_",MAX,"_SHIFTWINDOW_",SHFT,"_WINSIZE_",WINSIZE,
                 "/DIFFERENTIAL_K27_ON_",DIF_K27_BIN,"/")

# COMPUTE PERCENTAGE OF OVERLAPS ON A SLIDING WINDOW
myL <- NULL
mybait <- c("all")
myfilt <- c(paste0("dwnK27b",DIF_K27_BIN),paste0("upK27b",DIF_K27_BIN))#,"upK27bM500")
myanchor <- c("bgbord","nobgbord")
for(anchor in myanchor){
  for(bait in mybait){
    for(filt in myfilt){
      myL[[bait]][[anchor]][[filt]] <- fisherDistPlotBothSided(anchor,bait,filt,4000,50,250)
    }
  }
}

# CREATE DATA FRAME TO PLOT WITH GPLOT
for(FILT in c("upK27bM500")){
  
  bg_gafnobg.df <- myL$all$bgbord[[FILT]]
  bg_gafnobg.df <- cbind(bg_gafnobg.df,NM="BG_BORD")
  nobg_gafnobg.df <- myL$all$nobgbord[[FILT]]
  nobg_gafnobg.df <- cbind(nobg_gafnobg.df,NM="noBG_BORD")
  perc.df <- rbind(bg_gafnobg.df,nobg_gafnobg.df) 
  perc.df <- perc.df %>% replace_na(list(mean=1,sd=0.1))
  # GPLOT OPTS
  wdth <- 4000/(3*(4000/50))
  # PLOTS
  pdf(paste0(OUTF.d,"Comparison_BGbord_ALL_VS_noBGbord_ALL_",FILT,".pdf"),width=20)
  print(ggplot(perc.df, aes(x=x , y=mean,fill=NM)) + geom_bar(position=position_dodge(wdth+wdth/2),stat="identity",width = wdth)   +
          scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
          geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(wdth+wdth/2),size=0.5,width=.15) + 
          scale_fill_manual(values = c("dodgerblue","firebrick")) +
          xlim(0,1500))
    
  dev.off()
  tiff(paste0(OUTF.d,"Comparison_BGbord_ALL_VS_noBGbord_ALL_",FILT,".tiff"),units="in", width=10, height=5, res=300)
  print(ggplot(perc.df, aes(x=x , y=mean,fill=NM)) + geom_bar(position=position_dodge(1000),stat="identity",width = wdth)   +
          scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
          geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),position=position_dodge(1200),size=0.5,width=.15) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.text = element_blank(),
                legend.title = element_blank()) + theme_tiff+ scale_fill_manual(values = c("dodgerblue","firebrick")))
  dev.off()
  
}


#   -----------------------------------------------------------------------


# Figure 1E ---------------------------------------------------------------
SUB <- "FIGURES/Figure_1/Figure_1E"
MAX <- 60e3
SHFT <- 2e3
WINSIZE <- 1e3
DIF_K27_BIN <- "M500"
OUTF.d <- create(OUT.d,SUB,"/FISHER_TEST_CHIPSEQ_DIFFERENTIAL_K27_LEVELS/")

dfb <- NULL
dfg<- NULL
gaf_cp <- gaf_kc & cp & nobg
gaf_nocp <- gaf_kc & !cp & nobg
nogaf_cp <- !gaf_kc & cp & nobg
nogaf_nocp <- !gaf_kc & !cp & nobg
dfb <- NULL
dfg <- NULL
ctrl <- c("bg", "coh" , "gaf", "ctcf" ,"dref" ,"suhw")
for(line in ctrl){
  for(col in sort(unique(matrecap$zscorePM1kb_quintile))){
    typ <- "h"
    
    val <- matrecap$zscorePM1kb_quintile == col
    row <- line
    prot <- paste0(line,"_h")
    print(row)
    a_mat1 <- nrow(matrecap[val & get(line) ,])
    b_mat1 <- nrow(matrecap[val & !get(line),])
    a_mat2 <- nrow(matrecap[!val & get(line),])
    b_mat2 <- nrow(matrecap[!val & !get(line),])
    myFish.df <- myFisher(a_mat1,b_mat1,a_mat2,b_mat2,"two.sided")
    print(myFish.df)
    l_perc <- myFish.df$a_num
    lc_pv <- formatC(myFish.df$ft.pv, format = "e", digits = 2)
    odds <- round(myFish.df$ft.fc,2)
    lgt <- a_mat1
    sign <- sign(odds-1)
    dfb <- rbind.data.frame(dfb,cbind.data.frame(row=row,col=col,pv=lc_pv,odds=odds,perc=l_perc,type=typ,lgt=lgt,logpv=sign*(-log10(as.numeric(lc_pv)))))
    dfg <- rbind.data.frame(dfg,cbind.data.frame(row=prot,col=col,pv=lc_pv,odds=odds,perc=l_perc,type=typ,lgt=lgt,logpv=sign*(-log10(as.numeric(lc_pv)))))
    dfg$logpv <- ifelse(dfg$logpv > 2, 2,dfg$logpv)
    dfg$logpv <- ifelse(dfg$logpv < -2, -2,dfg$logpv)
  }
}

# Pvalue colored scale
this <- "INTER_MAT_Col=quint_k27_Line=ChIPseq"
ointer_this <- create(paste0(OUTF.d,"Quintile_k27/",this,"/All_genes_",this,"/"))
out_f <- paste0(ointer_this,"All_genes_",this)

dfg$row <- factor(dfg$row,levels = rev(levels(dfg$row)))

fisherGBplot(dfg,"Intersection plot : ChIPseq :",out_f,ratio=1)

#   -----------------------------------------------------------------------






