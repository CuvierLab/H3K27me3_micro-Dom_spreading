
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

# Figure 4C   -----------------------------------------------------------------------
SUB <- "FIGURES/Figure_4/Figure_4C/"
OUTF.d <- create(OUT.d,SUB)
DST_BORDER <- 5000
dataTYPE <- "WT_KD"
SPREAD <- "ALL"
MAXGAP <- 1500
MAXGAPTSS <- 2000
TYPE <- "RAW"
DIFF_SCORE <- "LFC"
NTILE <- 10
FDR <- 0.1
FDR_WT_KD <- 0.1
BS_WT_KD <- 40

ctl.bwp <- bw.ll[[TYPE]][["sig"]][1]
trt.bwp <- bw.ll[[TYPE]][["sig"]][2]


microWT_KD.gr <- subset(mDOM.dm3.gr,width < 1500 & width >= 80)
microWT_KD.gr$dtk27 <- mcols(distanceToNearest(microWT_KD.gr,Het.gr))$distance
microWT_KD.gr$dteuc <- mcols(distanceToNearest(microWT_KD.gr,Euc.gr))$distance

seqlevelsStyle(microWT_KD.gr) <- "Ensembl"

microWT_KD.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microWT_KD.gr)
microWT_KD.gr <- lfc_zs_dif_decile(anc.gr = microWT_KD.gr,NTILE=NTILE,ASYM=T)


microWT_KD_EUC.gr <- subset(microWT_KD.gr,dtk27 > DST_BORDER)
microWT_KD_HET.gr <- subset(microWT_KD.gr,dteuc > DST_BORDER)

microWT_KD_DOWN_EUC.gr <- GRanges(data.frame(microWT_KD_EUC.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) > NTILE )) 
microWT_KD_UP_EUC.gr <- GRanges(data.frame(microWT_KD_EUC.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) <= NTILE )) 
microWT_KD_DOWN_HET.gr <- GRanges(data.frame(microWT_KD_HET.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) > NTILE + 5)) 
microWT_KD_UP_HET.gr <- GRanges(data.frame(microWT_KD_HET.gr) %>% filter(UQ(as.name(paste0("dec",DIFF_SCORE))) < NTILE - 5)) 

# - OTHERS
microDOM_EUC.grl <- list(WT_KD=list(UP=microWT_KD_UP_EUC.gr,DOWN=microWT_KD_DOWN_EUC.gr,ALL=microWT_KD_EUC.gr))
microDOM_HET.grl <- list(WT_KD=list(UP=microWT_KD_UP_HET.gr,DOWN=microWT_KD_DOWN_HET.gr))


BS.l <- list(WT_KD=BS_WT_KD)
FDR.l <- list(WT_KD=FDR_WT_KD)

# - Associate noBg to all chipseq peaks
mylist <- c("cp","gaf","coh","ctcf")
for(nm in mylist){
  assign(paste0(nm,"_nobg"),get(nm) & !bg)
}


# BOOLEAN
`bg|cp|bs` <- bg | cp | bs
theme_update(plot.title = element_text(hjust = 0.5))
K27.gr <- normRxHMM.gr
K27.gr$len <- 1:length(K27.gr)
K27.gr$Lcnt <- 0
K27.gr$Rcnt <- 0

K27.grl <- split(K27.gr,seqnames(K27.gr))

# Avoid first/last domain of each chromosome, issue xhen doing +/- 1
K27_noEND.grl <- lapply(K27.grl,function(domChrX.gr){domChrX.gr[-c(1,length(domChrX.gr))]})
K27_noEND.gr <- Reduce("c",K27_noEND.grl)
Het_noEND.gr <- K27_noEND.gr[K27_noEND.gr$type=="Het"]
Euc_noEND.gr <- K27_noEND.gr[K27_noEND.gr$type=="Euc"]


# Features that will bracket
myFF <- "bg"


K27_TMP.gr <- K27.gr
allL.gr <- resize(K27_TMP.gr,1,"start")
allR.gr <- resize(K27_TMP.gr,1,"end")

allL.qH <- unique(queryHits(findOverlaps(allL.gr,mrc_tss.gr[get(myFF)],maxgap = MAXGAP)))
allR.qH <- unique(queryHits(findOverlaps(allR.gr,mrc_tss.gr[get(myFF)],maxgap = MAXGAP)))

K27_TMP.gr[allL.qH]$Lcnt <- 1
K27_TMP.gr[allR.qH]$Rcnt <- 1

het_idx <- Het_noEND.gr$len
euc_idx <- Euc_noEND.gr$len

Euc_noEND.gr$m_bord <- paste0(K27_TMP.gr[euc_idx-1]$Lcnt,
                              K27_TMP.gr[euc_idx]$Lcnt,
                              K27_TMP.gr[euc_idx]$Rcnt,
                              K27_TMP.gr[euc_idx+1]$Rcnt)

Euc_noEND.gr$LeftANDRight <- 0
Euc_noEND.gr$Left <- 0
Euc_noEND.gr$Right <- 0
Euc_noEND.gr$None <- 0

torep <- which(Euc_noEND.gr$m_bord=="0111" |
                 Euc_noEND.gr$m_bord=="1110" |
                 Euc_noEND.gr$m_bord=="1100" |
                 Euc_noEND.gr$m_bord=="0011" |
                 Euc_noEND.gr$m_bord=="1101" |
                 Euc_noEND.gr$m_bord=="1011" |
                 Euc_noEND.gr$m_bord=="1111" 
)
if(length(torep))Euc_noEND.gr[torep]$LeftANDRight <- 1



torep <- which(Euc_noEND.gr$m_bord=="0001"|
                 Euc_noEND.gr$m_bord=="1000"| 
                 Euc_noEND.gr$m_bord=="1010"| 
                 Euc_noEND.gr$m_bord=="1001")
if(length(torep))Euc_noEND.gr[torep]$Left <- 1

torep <- which(Euc_noEND.gr$m_bord=="0100"|
                 Euc_noEND.gr$m_bord=="0101"| 
                 Euc_noEND.gr$m_bord=="0010"| 
                 Euc_noEND.gr$m_bord=="0110")
if(length(torep))Euc_noEND.gr[torep]$Right <- 1

torep <- which(Euc_noEND.gr$m_bord=="0000")
if(length(torep))Euc_noEND.gr[torep]$None <- 1

mydf <- NULL
bracket.l <- c("Left","Right","LeftANDRight","None")
for(feat in c("all","cp","bs","bg|cp|bs")){
  mydf <- rbind.data.frame(mydf,cbind.data.frame(feat=feat,
                                                 pv=1,lfcb=0,
                                                 sign=" ",
                                                 x=bracket.l))
  
}

for(OVLPTYPE in c("all","cp","bs","bg|cp|bs")){
  
  microDOM_.gr <- microDOM_EUC.grl[[dataTYPE]][[SPREAD]]
  BS <- BS.l[[dataTYPE]]
  FDR <- FDR.l[[dataTYPE]]
  
  microDOM.gr <- subsetByOverlaps(microDOM_.gr,mrc_tss.gr[get(OVLPTYPE)],maxgap=MAXGAPTSS)
  mydf$feat <- as.character(mydf$feat)
  mydf$sign <- as.character(mydf$sign)
  mydf$x <- as.character(mydf$x)
  
  
  Euc_noEND.df <- as.data.frame(Euc_noEND.gr)
  Euc_noEND.df$micro <- 0
  row.IDX <- unique(queryHits(findOverlaps(Euc_noEND.gr,microDOM.gr)))
  Euc_noEND.df[row.IDX,"micro"] <- 1
  
  Euc_noEND.df$m_bord <- as.character(Euc_noEND.df$m_bord)
  for(comb in c("Left","Right","LeftANDRight")){
    a_mat1 <- as.integer(Euc_noEND.df %>% filter(UQ(as.name(paste0(comb)))==1,micro==1) %>% summarise(n()))
    b_mat1 <- as.integer(Euc_noEND.df %>% filter(UQ(as.name(paste0("None")))==1,micro==1) %>% summarise(n()))
    a_mat2 <- as.integer(Euc_noEND.df %>% filter(UQ(as.name(paste0(comb)))==1,micro==0) %>% summarise(n()))
    b_mat2 <- as.integer(Euc_noEND.df %>% filter(UQ(as.name(paste0("None")))==1,micro==0) %>% summarise(n()))
    
    myFish.df <- myFisher(a_mat1,a_mat2,b_mat1,b_mat2,"two.sided")
    
    mypv <- formatC(as.numeric(myFish.df$ft.pv, format = "e", digits = 2))
    lfc <- round(log2(myFish.df$ft.fc),2)
    mysign <- ifelse(myFish.df$ft.pv <= 0.05,
                     ifelse(myFish.df$ft.pv < 0.01,
                            ifelse(myFish.df$ft.pv < 0.001,"***","**"),"*"),
                     ifelse(myFish.df$ft.pv > .05 & is.infinite(lfc)," ", "NS"))
    
    
    mydf <- mydf %>% mutate_cond(feat==OVLPTYPE & x==comb ,pv=mypv ,lfcb=lfc,sign=mysign)
  }
}


mydf$intSign <- sign(mydf$lfcb)
mydf <- mydf %>% mutate(lfc=ifelse(abs(lfcb) > 2,intSign*2,intSign*abs(lfcb)))
mydf <- mydf %>% mutate(lfc=ifelse(is.infinite(lfcb), NA,lfc))
colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))



pdf(paste0(OUTF.d,"Heatmap_Ftest_FeaturesCOMBINATION_MICROTYPE_",myFF,".pdf"))
p <- ggplot(mydf, aes(x = x, y = feat,fill=lfc)) +      geom_tile(col="black") +
  theme( panel.background=element_rect(fill = "white"),
         panel.grid.major=element_blank(),
         panel.border = element_rect(fill = NA,
                                     colour = "grey20"))+
  # panel.grid.minor=element_blank()) +
  coord_equal(ratio=.5) +
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),position = "top") +
  geom_text(aes(label =sign)) +
  scale_fill_gradientn(colours = colfunc(3) ,
                       guide = guide_colorbar(barwidth = 0.8,
                                              title = 'Log2 OddsRatio',
                                              label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                              barheight = 10,
                                              nbin = 10,
                                              draw.ulim = FALSE,
                                              draw.llim = FALSE,
                                              ticks = FALSE),
                       breaks=c(-2,0.2,2),
                       labels=c("-1.5","0","1.5"),
                       limits=c(-2,2),na.value = "white")
print(p)
dev.off()




#   -----------------------------------------------------------------------



# Figure 4D  -----------------------------------------------------------------------
SUB <- "FIGURES/Figure_4/Figure_4D/"
OUTA.d <- create(OUT.d,SUB)
OUTSCA.d <- create(OUT.d,SUB,"Scatterplot/")

constraint.gr <- corces.tad.gr 
bin <- 41
hic_res <- 1000
APA_SIZE <- 41
mBin <- (APA_SIZE-1)/2

# Coordinate of central pixel
i_CP <- (bin-1)/2+1;j_CP <- (bin-1)/2+1 # Center bin interaction of size (1bin;1bin)
# Coordinate of 3x3 square centered on central pixel
m_size <- (bin-1)/2 # Middle pixel value - 1
i_CS <- m_size:(m_size+2);j_CS <- m_size:(m_size+2)
# Coordinate of Upper left 3x3 square from 3x3 square at center
i_ULS <- i_CS-3;j_ULS <- j_CS-3
# Coordinate of  Upper Right 3x3 square from 3x3 square at center (Compartment)
i_URS <- i_CS-3;j_URS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BRS <- i_CS+3;j_BRS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BLS <- i_CS+3;j_BLS <- j_CS-3






# Create a Granges of features as anchors
anchor.gr <- mrc_tss.gr
anchor.gr$a_IDX <- 1:length(anchor.gr)

# BAIT / INTERACTION POINT 
bait.gr <- mrc_tss.gr
bait.gr$i_IDX <- 1:length(bait.gr)

# CREATE INTERACTION DATA FRAME
apaCPLS.df <- createGRcouplesAPA(anc.gr = anchor.gr,int.gr = bait.gr,hic_res = hic_res,SQF = sqfDM3,constraint.gr = constraint.gr,MYCHR = kc.chr,apa_size = 41)
# AS THE BLS NEED COUPLES TO BE >= 14 SQUARES AWAY FILTER THEM
apaCPLS.df <- subset(apaCPLS.df,dst_bin >= 14)

# Sample some of the interaction 
set.seed(123)
apaCPLS_SAMP.df <- apaCPLS.df[sample(1:nrow(apaCPLS.df),100e3,F),]
nrow(apaCPLS_SAMP.df)

# EXCTRACT INTERACTION BETWEEN ALL DF COUPLES
l_IN <- readRDS(HIC.l$WT)
NM <- l_IN$NM
hic_RES.lst <- hic_pipeline(l_IN=l_IN,apaCPLS.df = apaCPLS_SAMP.df ,mBin = mBin)
OUT_METRIC.d <- create(paste0(OUTA.d,"/METRICS_HIC/WT/"))
saveRDS(hic_RES.lst,paste0(OUT_METRIC.d,"/TSS_vs_TSS_interaction.rds"))
gc()

l_IN <- readRDS(HIC.l$KD)
NM <- l_IN$NM
hic_RES.lst <- hic_pipeline(l_IN=l_IN,apaCPLS.df = apaCPLS_SAMP.df ,mBin = mBin)
OUT_METRIC.d <- create(paste0(OUTA.d,"/METRICS_HIC/KD/"))
saveRDS(hic_RES.lst,paste0(OUT_METRIC.d,"/TSS_vs_TSS_interaction.rds"))
gc()

GR.d <- create(paste0(OUTA.d,"/METRICS_HIC/ANCHOR_INTERACTION_GR/"))
saveRDS(anchor.gr,paste0(GR.d,"anchor.rds"))
saveRDS(bait.gr,paste0(GR.d,"bait.rds"))



# LOAD
OUT_METRIC.d <- paste0(OUTA.d,"/METRICS_HIC/WT/")
GR.d <- paste0(OUTA.d,"/METRICS_HIC/ANCHOR_INTERACTION_GR/")
hic_RES.lst <- readRDS(paste0(OUT_METRIC.d,"/TSS_vs_TSS_interaction.rds"))
anchor.gr <- readRDS(paste0(GR.d,"anchor.rds"))
bait.gr <- readRDS(paste0(GR.d,"bait.rds"))



myLST <- c("bg","cp","gaf","ctcf","coh")
l_DF <- hic_RES.lst$METRICS_QUANT.df


MYDF <- NULL
tst <- combn(myLST,2,
      function(X){
        COMB <- paste0(X[1],"_",X[2])
        
        A_IDX <- anchor.gr[get(X[1])]$a_IDX
        I_IDX <- bait.gr[get(X[2])]$i_IDX
        AR_IDX <- anchor.gr[get(X[2])]$a_IDX
        IR_IDX <- bait.gr[get(X[1])]$i_IDX
        
        CP_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX))|((a_IDX %in% AR_IDX) & (i_IDX %in% IR_IDX)))$CP_V,na.rm=T)
        CS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX))|((a_IDX %in% AR_IDX) & (i_IDX %in% IR_IDX)))$CS_V,na.rm=T)
        
        BLS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX))|((a_IDX %in% AR_IDX) & (i_IDX %in% IR_IDX)))$BLS_V,na.rm=T)
        MYDF <<- rbind.data.frame(MYDF,
                                  cbind.data.frame(CPL="HET",
                                                   NM=COMB,
                                                   CP_V_M=CP_V_M,
                                                   BLS_V_M=BLS_V_M,
                                                   CS_V_M=CS_V_M))
        
      })

tst <- combn(myLST,1,
             function(X){
               COMB <- X[1]
               print(COMB)
               A_IDX <- anchor.gr[get(X[1])]$a_IDX
               I_IDX <- bait.gr[get(X[1])]$i_IDX
             
               CP_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$CP_V,na.rm=T)
               CS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$CS_V,na.rm=T)
               BLS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$BLS_V,na.rm=T)
               print(CP_V_M)
               print(BLS_V_M)
               print(CS_V_M)
               MYDF <<- rbind.data.frame(MYDF,
                                         cbind.data.frame(CPL="HOM",
                                                          NM=COMB,
                                                          CP_V_M=CP_V_M,
                                                          BLS_V_M=BLS_V_M,
                                                          CS_V_M=CS_V_M))
             },simplify = F)

for(CTL in 1:1000){
  A_IDX <- sample(1:nrow(l_DF),10000,T)
  I_IDX <- sample(1:nrow(l_DF),10000,T)
  CP_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$CP_V,na.rm=T)
  CS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$CS_V,na.rm=T)
  BLS_V_M <- mean(subset(l_DF,((a_IDX %in% A_IDX) & (i_IDX %in% I_IDX)))$BLS_V,na.rm=T)
  MYDF <<- rbind.data.frame(MYDF,
                            cbind.data.frame(CPL="CTL",
                                             NM=" ",
                                             CP_V_M=CP_V_M,
                                             BLS_V_M=BLS_V_M,
                                             CS_V_M=CS_V_M))
}


head(MYDF)
CS_V_Q <- quantile(subset(MYDF,CPL=="CTL")$CS_V,seq(0,1,0.01),na.rm=T)[95]
CP_V_Q <- quantile(subset(MYDF,CPL=="CTL")$CP_V,seq(0,1,0.01),na.rm=T)[95]
BLS_V_Q <- quantile(subset(MYDF,CPL=="CTL")$BLS_V,seq(0,1,0.01),na.rm=T)[95]

MYDF_T <- subset(MYDF,CPL=="CTL")
MYDF_T$CPL <- as.factor("CTL")
MYDF_A <- subset(MYDF,CPL!="CTL")
MYDF_T$CPL <- as.factor("CTL")

pdf(paste0(OUTSCA.d,"All_scatterplot_all_comb_insulators.pdf"))
pGGS <- ggscatter(MYDF_T,x="BLS_V_M",y="CS_V_M",color = "CPL",size=1, font.label = c(6, "plain"),repel = F,palette = c("grey")) +
  ggtitle("MEAN RAW DATA") + coord_cartesian(xlim=c(230,290),ylim=c(230,290))+
  geom_vline(aes(xintercept=BLS_V_Q),linetype="dotted") + geom_hline(aes(yintercept=CS_V_Q),linetype="dotted")
ggscatter(MYDF_A,x="BLS_V_M",y="CS_V_M",ggp=pGGS,color = "CPL",label="NM",size=1,font.label = c(6, "plain"),repel = T,palette = c("grey","dodgerblue","firebrick"))

pGGS <- ggscatter(MYDF_T,x="BLS_V_M",y="CS_V_M",color = "CPL",size=1, font.label = c(6, "plain"),repel = F,palette = c("grey")) +
  ggtitle("MEAN RAW DATA") + coord_cartesian(xlim=c(230,290),ylim=c(230,290))+
  geom_vline(aes(xintercept=BLS_V_Q),linetype="dotted") + geom_hline(aes(yintercept=CS_V_Q),linetype="dotted")
ggscatter(MYDF_A,x="BLS_V_M",y="CS_V_M",ggp=pGGS,color = "CPL",size=1,font.label = c(6, "plain"),repel = T,palette = c("grey","dodgerblue","firebrick"))




pGGS <- ggscatter(MYDF_T,x="BLS_V_M",y="CP_V_M",color = "CPL",size=1, font.label = c(6, "plain"),repel = F,palette = c("grey")) +
  ggtitle("MEAN RAW DATA") + coord_cartesian(xlim=c(230,290),ylim=c(230,330))+
  geom_vline(aes(xintercept=BLS_V_Q),linetype="dotted") + geom_hline(aes(yintercept=CS_V_Q),linetype="dotted")
ggscatter(MYDF_A,x="BLS_V_M",y="CP_V_M",ggp=pGGS,color = "CPL",label="NM",size=1,font.label = c(6, "plain"),repel = T,palette = c("grey","dodgerblue","firebrick"))

pGGS <- ggscatter(MYDF_T,x="BLS_V_M",y="CP_V_M",color = "CPL",size=1, font.label = c(6, "plain"),repel = F,palette = c("grey")) +
  ggtitle("MEAN RAW DATA") + coord_cartesian(xlim=c(230,290),ylim=c(230,330))+
  geom_vline(aes(xintercept=BLS_V_Q),linetype="dotted") + geom_hline(aes(yintercept=CS_V_Q),linetype="dotted")
ggscatter(MYDF_A,x="BLS_V_M",y="CP_V_M",ggp=pGGS,color = "CPL",size=1,font.label = c(6, "plain"),repel = T,palette = c("grey","dodgerblue","firebrick"))

dev.off()

#   -----------------------------------------------------------------------



# Figure 4E ---------------------------------------------------------------
SUB <- "FIGURES/Figure_4/Figure_4E/"
OUTF.d <- create(OUT.d,SUB)
WDTH_MAX <- 1500
WDTH_MIN <- 40
APA_SIZE <- 41
mBin <- (APA_SIZE-1)/2
filt_DIST <- T
MIN_K27 <- 60
MAXGAP_TSS_BORD <- 1000
TYPE <- "RPKM"

# COORDINATE
i_CP <- (APA_SIZE-1)/2+1;j_CP <- (APA_SIZE-1)/2+1 # Center APA_SIZE interaction of size (1APA_SIZE;1APA_SIZE)
# Coordinate of 3x3 square centered on central pixel
m_size <- (APA_SIZE-1)/2 # Middle pixel value - 1
i_CS <- m_size:(m_size+2);j_CS <- m_size:(m_size+2)
# Coordinate of Upper left 3x3 square from 3x3 square at center
i_ULS <- i_CS-3;j_ULS <- j_CS-3
# Coordinate of  Upper Right 3x3 square from 3x3 square at center (Compartment)
i_URS <- i_CS-3;j_URS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BRS <- i_CS+3;j_BRS <- j_CS+3
# Coordinate of  Bottom Right 3x3 square from 3x3 square at center (Compartment)
i_BLS <- i_CS+3;j_BLS <- j_CS-3



ctl.bwp <- bw.ll[[TYPE]][["sig"]][1]
trt.bwp <- bw.ll[[TYPE]][["sig"]][2]


# MICRO DOMAIN LOADING
microWT_KD.gr <- subset(mDOM.dm3.gr,width < WDTH_MAX & width >= WDTH_MIN)
seqlevelsStyle(microWT_KD.gr) <- "UCSC"
microWT_KD.gr <- bw_signal(trt.bwp = trt.bwp,ctl.bwp = ctl.bwp,anc.gr = microWT_KD.gr)
microWT_KD.gr <- lfc_zs_dif_decile(anc.gr = microWT_KD.gr,NTILE=NTILE,ASYM=T)
seqlevelsStyle(microWT_KD.gr) <- "Ensembl"
microWT_KD_DOWN.gr <- subset(microWT_KD.gr,decZS > 10)
microWT_KD_UP.gr <- subset(microWT_KD.gr,decZS <= 10)

# HIC LOADING
l_IN <- readRDS(HIC.l$WT)
HIC.WT.lmtx <- l_IN$HIC
HIC.WT.res <- l_IN$RES
HIW.WT.NM <- l_IN$NM
l_IN <- readRDS(HIC.l$KD)
HIC.KD.lmtx <- l_IN$HIC
HIC.KD.res <- l_IN$RES
HIW.KD.NM <- l_IN$NM
RES <- HIC.WT.res

# EXCTRACT GSEA FROM TAD SUMMIT INTERACTION DEPENDING ON MDOM PRESENCE
ramirez.tad.dm3.gr <- addSeqinfo(ramirez.all.tad.gr,"dm3","s2")
myTAD.l <- list(RAMI=ramirez.tad.dm3.gr)
for(TADN in names(myTAD.l)){
  NM_TAD_O <- toupper(gsub("(.*).tad.dm3.gr","\\1",TADN))
  tad.gr <- myTAD.l[[TADN]]
  tad.dm3.gr$IDX <- 1:length(tad.dm3.gr)
  mydom.grl <- split(tad.dm3.gr,seqnames(tad.dm3.gr))
  tad.noXtrm.dm3.gr <- Reduce("c",(lapply(mydom.grl,function(X){X[-c(1,length(X))]})))
  
  folp <- findOverlapPairs(tad.noXtrm.dm3.gr,Het.gr)
  tad.K27.gr <- subsetByOverlaps(tad.noXtrm.dm3.gr,Het.gr)
  tad.noK27.gr <- subsetByOverlaps(tad.noXtrm.dm3.gr,Het.gr,invert=T)
  
  pint <- pintersect(folp)
  folp.df <- as.data.frame(folp)
  folp.df$inter <- width(pint)
  folp.df$nme <- paste0(folp.df$first.X.start,"_",folp.df$first.X.end)
  dens.df <- folp.df %>% dplyr::select(nme,inter) %>% group_by(nme) %>% summarise(tot=sum(inter))
  
  tad.K27.gr$nm <- paste0(start(tad.K27.gr),"_",end(tad.K27.gr))
  
  print(NM)
  tad.K27.gr$K27wdth <- 0 
  tad.K27.gr[match(dens.df$nme,tad.K27.gr$nm)]$K27wdth <- dens.df$tot 
  tad.K27.gr$K27dens <- round(tad.K27.gr$K27wdth*100/width(tad.K27.gr))
  
  tad.noK27.gr$nm <- paste0(start(tad.noK27.gr),"_",end(tad.noK27.gr))
  tad.noK27.gr$K27wdth <- 0
  tad.noK27.gr$K27dens <- 0
  tad.gr <- sort(c(tad.noK27.gr,tad.K27.gr))
  
  tad.gr$IDX <- 1:length(tad.gr)
  # CREATE COUPLES OF INTERACTION TO DUMP
  tadL.gr <- bordersGR(tad.gr)$leftBORD
  tadR.gr <- bordersGR(tad.gr)$rightBORD
  TADTMP.df <- data.frame(LBORD=ceiling(start(tadL.gr)/RES),RBORD=ceiling(start(tadR.gr)/RES),IDX=tadL.gr$IDX,CHR=seqnames(tadL.gr))
  TAD.df <- NULL
  
  
  
  for(chr in s2.chr){
    TMP.df <- subset(TADTMP.df,CHR==chr)
    DIM <- as.numeric(ceiling(seqlengths(sqfDM3[paste0("chr",chr)])/RES))
    print(DIM)
    toDEL <- which(TMP.df$LBORD > DIM-(mBin) | TMP.df$LBORD < (mBin) | TMP.df$RBORD > DIM-(mBin) | TMP.df$RBORD < (mBin)  )
    if(length(toDEL)) TMP.df <- TMP.df[-toDEL,]
    toDEL <- which(abs(TMP.df$LBORD - TMP.df$RBORD) < APA_SIZE  )
    if(length(toDEL)) TMP.df <- TMP.df[-toDEL,]
    TAD.df <- rbind(TAD.df,TMP.df)
  }
  
  
  
  # EXCTRACT HIC MATRICES BINxBIN
  myMAT <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  myMAT <- foreach::foreach(X=1:nrow(TAD.df)) %dopar% {
    binANC <- TAD.df[X,1]
    binINT <- TAD.df[X,2]
    CHR <- TAD.df[X,4]
    if(binINT < binANC){
      TMP <- binINT
      binINT <- binANC
      binANC <- TMP
    }
    # Compute the APA matrix and return it in myMAT as a list of matrix
    feat_apa.WT.mtx <- Matrix(HIC.WT.lmtx[[CHR]][((binANC-mBin):(binANC+mBin)),((binINT-mBin):(binINT+mBin))],sparse=F)
    feat_apa.WT.mtx[!is.finite(feat_apa.WT.mtx) | feat_apa.WT.mtx==0] <- NA
    
    feat_apa.KD.mtx <- Matrix(HIC.KD.lmtx[[CHR]][((binANC-mBin):(binANC+mBin)),((binINT-mBin):(binINT+mBin))],sparse=F)
    #rremove NaN and 0
    feat_apa.KD.mtx[!is.finite(feat_apa.KD.mtx) | feat_apa.KD.mtx==0] <- NA
    
    
    list(IDX=TAD.df[X,3],HIC.WT=as.matrix(feat_apa.WT.mtx),HIC.KD=as.matrix(feat_apa.KD.mtx))
  }
  stopImplicitCluster()
  
  # EXCTRACT APA METRICS AS A DATAFRAME
  # OBS/EXP VERSION
  TAD_OE.df <- NULL 
  invisible(lapply(myMAT, 
                   function(X){
                     IDX <- X[[1]]
                     TYPE <- c("WT","KD")
                     for(V in 2:3){
                       MAT <- X[[V]]
                       l_ <- list()
                       l_[["CP_V"]] <- MAT[i_CP,j_CP]
                       l_[["CS_V"]] <- MAT[i_CS,j_CS]
                       l_[["ULS_V"]] <- MAT[i_ULS,j_ULS]
                       l_[["URS_V"]] <- MAT[i_URS,j_URS]
                       l_[["BLS_V"]] <- MAT[i_BLS,j_BLS]
                       l_[["BRS_V"]] <- MAT[i_BRS,j_BRS]
                       TMP.df <- data.frame(IDX=IDX)
                       for(PAR in names(l_)){
                         TMP.df[,PAR] <- mean(as.array(l_[[PAR]]),na.rm=T)
                       }
                       TMP.df[,"TYPE"] <- TYPE[V-1]
                       TAD_OE.df <<- rbind(TAD_OE.df,TMP.df)
                     }
                     
                   }))
  
  
  
  # QUANTILE VERSION
  TAD_OE_QUANT.df <- NULL 
  invisible(lapply(myMAT, 
                   function(X){
                     IDX <- X[[1]]
                     TYPE <- c("WT","KD")
                     for(V in 2:3){
                       MAT_QUANT <- X[[V]]
                       MAT_QUANT <- matrix(ntile(MAT_QUANT,500),APA_SIZE,APA_SIZE)
                       l_ <- list()
                       l_[["CP_V"]] <- MAT_QUANT[i_CP,j_CP]
                       l_[["CS_V"]] <- MAT_QUANT[i_CS,j_CS]
                       l_[["ULS_V"]] <- MAT_QUANT[i_ULS,j_ULS]
                       l_[["URS_V"]] <- MAT_QUANT[i_URS,j_URS]
                       l_[["BLS_V"]] <- MAT_QUANT[i_BLS,j_BLS]
                       l_[["BRS_V"]] <- MAT_QUANT[i_BRS,j_BRS]
                       TMP.df <- data.frame(IDX=IDX)
                       for(PAR in names(l_)){
                         TMP.df[,PAR] <- round(mean(as.array(l_[[PAR]]),na.rm=T))
                       }
                       TMP.df[,"TYPE"] <- TYPE[V-1]
                       TAD_OE_QUANT.df <<- rbind(TAD_OE_QUANT.df,TMP.df)
                     }
                   }))
  
  # CREATE DATA FRAME WITH FEATURES BOUND AT BORDERS
  GLOBAL.df <- data.frame(mcols(tad.gr))
  GLOBAL.df$CHR <- as.vector(seqnames(tad.gr))
  
  
  GLOBAL.df$BEAF_L <- 0 ; GLOBAL.df$BEAF_R <- 0
  GLOBAL.df$GAF_L <- 0 ; GLOBAL.df$GAF_R <- 0
  GLOBAL.df$GAFKC_L <- 0 ; GLOBAL.df$GAFKC_R <- 0
  GLOBAL.df$CP190_L <- 0 ; GLOBAL.df$CP190_R <- 0
  GLOBAL.df$CTCF_L <- 0 ; GLOBAL.df$CTCF_R <- 0
  GLOBAL.df$COH_L <- 0 ; GLOBAL.df$COH_R <- 0
  GLOBAL.df$MDOM <- 0
  GLOBAL.df$ALL <- 1
  
  
  
  matrix(ntile(as.vector(myMAT[[2]]$HIC.WT),500),41,41)
  
  # Associate tss from mrc with TAD border L& R
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[bg],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[bg],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"BEAF_L"] <- 1
  GLOBAL.df[fol_R,"BEAF_R"] <- 1
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[gaf],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[gaf],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"GAF_L"] <- 1
  GLOBAL.df[fol_R,"GAF_R"] <- 1
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[gaf_kc],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[gaf_kc],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"GAFKC_L"] <- 1
  GLOBAL.df[fol_R,"GAFKC_R"] <- 1
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[cp],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[cp],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"CP190_L"] <- 1
  GLOBAL.df[fol_R,"CP190_R"] <- 1
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[ctcf],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[ctcf],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"CTCF_L"] <- 1
  GLOBAL.df[fol_R,"CTCF_R"] <- 1
  fol_L <- unique(queryHits(findOverlaps(tadL.gr,mrc_tss.gr[coh],maxgap=MAXGAP_TSS_BORD)))
  fol_R <- unique(queryHits(findOverlaps(tadR.gr,mrc_tss.gr[coh],maxgap=MAXGAP_TSS_BORD)))
  GLOBAL.df[fol_L,"COH_L"] <- 1
  GLOBAL.df[fol_R,"COH_R"] <- 1
  
  # 
  md_fol <- unique(queryHits(findOverlaps(tad.gr,microWT_KD_DOWN.gr)))
  GLOBAL.df[md_fol,"MDOM"] <- 1
  
  
  GLOBAL.df$MDOM_EUC_RIGHT_TAD <- 0
  
  for(IDX in 1:(nrow(GLOBAL.df)-1)){
    if(GLOBAL.df[IDX,"K27dens"]>50 & GLOBAL.df[IDX+1,"MDOM"]==1 & GLOBAL.df[IDX,"CHR"]==GLOBAL.df[IDX+1,"CHR"]){
      GLOBAL.df[IDX,"MDOM_EUC_RIGHT_TAD"] <- 1
    }
  }
  
  # MERGE WITH HIC METRIC
  GLOBAL_TAD_OE.df <- merge(GLOBAL.df,TAD_OE.df,by="IDX")
  GLOBAL_TAD_OE.df[is.na(GLOBAL_TAD_OE.df)] <- 0  
  GLOBAL_TAD_OE_QUANT.df <- merge(GLOBAL.df,TAD_OE_QUANT.df,by="IDX")
  GLOBAL_TAD_OE_QUANT.df[is.na(GLOBAL_TAD_OE_QUANT.df)] <- 0  
  
  # DO THE GSEA ON RAW HIC METRICS
  METRICS <- "BLS_V"
  
  # DO THE GSEA ON QUANTILIZED HIC METRICS
  for(FILTERS in c("ALL")){
    OUTG.d <- create(paste0(OUTF.d,"/GSEA_DIFFERENTIAL/"))
    for(metrics in METRICS){
      
      FILTER <- "ALL"
      # EXCTRACT WT AND KD METRICS VALUES
      WT.TMP <- subset(GLOBAL_TAD_OE_QUANT.df,TYPE=="WT")
      KD.TMP <- subset(GLOBAL_TAD_OE_QUANT.df,TYPE=="KD")
      
      #COMPUTE RATIO
      RATIO <- WT.TMP[,metrics]-KD.TMP[,metrics]
      WT.TMP$RATIO <- RATIO
      
      # CREATE GSEA RANKED BY THE METRICS AND NAM THE LIST WITH TAD IDX
      DF.TMP <- WT.TMP
      DF.TMP <- DF.TMP[order(RATIO,decreasing = T),]
      GSEA <- DF.TMP[,"RATIO"]
      names(GSEA) <- DF.TMP[,"IDX"]
      GSEA.lst <- NULL
      GSEA.lst[["TOT"]] <- DF.TMP[DF.TMP$MDOM_EUC_RIGHT_TAD==1 & DF.TMP[,FILTERS]==1,"IDX"]
      # GSEA.lst[["TOT"]] <- c(head(names(GSEA),100))#,tail(names(GSEA),100))
      # GSEA.lst[["TOT"]] <- names(GSEA)[seq(1,length(GSEA),by=2)]
      # GSEA.lst[["TOT"]] <- sample(names(GSEA),200,replace = F)
      
      
      fgseaRes <- fgsea(pathways = GSEA.lst, 
                        stats = GSEA,
                        minSize=5,
                        maxSize=max(unlist(lapply(GSEA.lst,length))),
                        nperm=10000,nproc=8)
      
      sink(paste0(OUTG.d,FILTERS,"_",toupper(metrics),"_GSEA_RANKED_TAD_CENTERED_MDOM_DOWN_ON_EUCHROMATIN_TAD_ON_LEFT.txt"))
      print(fgseaRes)
      sink()
      
      pdf(paste0(OUTG.d,FILTERS,"_",toupper(metrics),"_GSEA_RANKED_TAD_CENTERED_MDOM_DOWN_ON_EUCHROMATIN_TAD_ON_LEFT.pdf"))
      print(plotEnrichment(GSEA.lst[["TOT"]],GSEA) + labs(title=paste0("GSEA on MDOM DOWN")))
      dev.off()
      
    }
  }
}



#   -----------------------------------------------------------------------










