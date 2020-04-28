
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

# Figure 3A   -----------------------------------------------------------------------
SUB <- "FIGURES/Figure_3/Figure_3B/"
OUTF.d <- create(OUT.d,SUB)



# GAF x CP190

dfb <- NULL
dfg<- NULL
gaf_cp <- gaf_kc & cp & nobg
gaf_nocp <- gaf_kc & !cp & nobg
nogaf_cp <- !gaf_kc & cp & nobg
nogaf_nocp <- !gaf_kc & !cp & nobg
ctrl <- c("gaf_cp","gaf_nocp","nogaf_cp","nogaf_nocp")
for(line in ctrl){
  for(col in sort(unique(matrecap$zscorePM1kb_quintile))){
    typ <- "h"
    val <- matrecap$zscorePM1kb_quintile == col
    row <- line
    prot <- paste0(line,"_h")
    a_mat1 <- nrow(matrecap[val & get(line) & nobg,])
    b_mat1 <- nrow(matrecap[val & !get(line)& nobg,])
    a_mat2 <- nrow(matrecap[!val & get(line)& nobg,])
    b_mat2 <- nrow(matrecap[!val & !get(line)& nobg,])
    myFish.df <- myFisher(a_mat1,b_mat1,a_mat2,b_mat2,"two.sided")
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
NM <- "INTER_MAT_Col=quint_k27_Line=Gaf_CoFac=CP190_NOBGG"
OUTFN.p <- paste0(OUTF.d,"/",NM)
fisherGBplot(dfg,"Intersection plot : Beaf & Gaf :",OUTFN.p,ratio=.5)



# GAF x RAD21

dfb <- NULL
dfg<- NULL
gaf_coh <- gaf_kc & coh & nobg
gaf_nocoh <- gaf_kc & !coh & nobg
nogaf_coh <- !gaf_kc & coh & nobg
nogaf_nocoh <- !gaf_kc & !coh & nobg
ctrl <- c("gaf_coh","gaf_nocoh","nogaf_coh","nogaf_nocoh")
for(line in ctrl){
  for(col in sort(unique(matrecap$zscorePM1kb_quintile))){
    typ <- "h"
    val <- matrecap$zscorePM1kb_quintile == col
    row <- line
    prot <- paste0(line,"_h")
    a_mat1 <- nrow(matrecap[val & get(line) & nobg,])
    b_mat1 <- nrow(matrecap[val & !get(line)& nobg,])
    a_mat2 <- nrow(matrecap[!val & get(line)& nobg,])
    b_mat2 <- nrow(matrecap[!val & !get(line)& nobg,])
    myFish.df <- myFisher(a_mat1,b_mat1,a_mat2,b_mat2,"two.sided")
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
NM <- "INTER_MAT_Col=quint_k27_Line=Gaf_CoFac=RAD21_NOBG"
OUTFN.p <- paste0(OUTF.d,"/",NM)
fisherGBplot(dfg,"Intersection plot : Beaf & Gaf :",OUTFN.p,ratio=.5)
#   -----------------------------------------------------------------------



