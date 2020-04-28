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
source(paste0(OUT.d,"/LIBS/mrc_lib.R"))
source(paste0(OUT.d,"/LIBS/data_lib.R")) # load data
TMP.d <- create(paste0(SPATH,"/TMP/"))

#####################################################################################-
#
#          LOAD DATA
#
#####################################################################################-
# 1. STEP  : TURN JUICER EIGEN VECTOR TO WIG FILE ----
OUT.NM <- "beafWT"
IN.d <- paste0(OUT.d,"/PROCESSING_SCRIPTS/COMPARTMENT_COMPUTATION/Compartments/")
OUTC.d <- create(paste0(OUT.d,"/PROCESSING_SCRIPTS/COMPARTMENT_COMPUTATION/Compartments_WIG/",OUT.NM,"/"))

for(comp in dir(IN.d,full.names = T)){
  if(!grepl("*.txt",comp)){
    next()
  }
  res <- gsub(".*res([0-9]+)bp.*","\\1",basename(comp))
  chr <- paste0(gsub(".*(compartment_|compartment_chr)(.+)_res.*","\\2",basename(comp)))
  NM <- paste0(gsub("(.*)_(compartment_|compartment_chr)(.+)_res.*","\\1",basename(comp)))
  comp_wig <- paste0(NM,"_compartment_res",res,"_S2_dm3.wig")
  header <- ifelse(file.exists(paste0(OUTC.d,comp_wig)), '','track type=wiggle_0 visibility=full \n')
  
  header=paste0(header,'fixedStep chrom=chr',chr,' start=1 step=',res,' span=',res)
  write.table(header,paste0(OUTC.d,comp_wig), quote=FALSE, sep="\n",
              row.names = FALSE,
              col.names = FALSE,
              append = T)
  comp_file <- read.table(comp)*1000
  comp_file[is.na(comp_file)] <- 0
  write.table(comp_file,paste0(OUTC.d,comp_wig), quote=FALSE, sep="\n",
              row.names = FALSE,
              col.names = FALSE,
              append = T)
  print(chr)
}
# ----

# 2. STEP  :  REANNOTATE COMPARTMENT. SOME CHROMOSOME MAY HAVE A COMP = +1 OTHERS -1 ----
# USE TSS ACTIVITY TO DO SO
IN.d <- paste0(OUT.d,"Compartments_WIG/",OUT.NM,"/")
OUTC.d <- create(paste0(OUT.d,"/PROCESSING_SCRIPTS/COMPARTMENT_COMPUTATION/Compartments_REANNOTATED/",OUT.NM,"/"))
act_tss.gr <- readRDS(paste0(OUT.d,"/PUBLISHED/actTSS_Kc167_dm3.GR.RData"))

for(comp in dir(IN.d,full.names = T)){
  res <- gsub(".*res([0-9]+)_.*","\\1",basename(comp))
  NM <- paste0(gsub("(.*)_(compartment_|compartment_chr)(.+)_res.*","\\1",basename(comp)))
  
  compAB.gr <- import(comp)
  compAB.gr <- addSeqinfo(compAB.gr)
  compAB.gr <- trim(compAB.gr)
  print(compAB.gr)
  compAB.grl <- split(compAB.gr,seqnames(compAB.gr))
  compA.gr <- NULL
  compB.gr <- NULL
  
  compAB_chrX.gr <- compAB.grl[[1]]
  new_COMP.grl <- lapply(compAB.grl,
         function(compAB_chrX.gr){
           print(compAB_chrX.gr)
           compAB_chrX.gr$comp <- ifelse(compAB_chrX.gr$score>=0,"A","B")
           compA_chrX.gr <- subset(compAB_chrX.gr,comp=="A")
           compB_chrX.gr <- subset(compAB_chrX.gr,comp=="B")
           chr <- as.character(unique(seqnames(compAB_chrX.gr)))
           act_tss_chrX.gr <- subset(act_tss.gr,seqnames==chr)
           #Check if this is a correct annotation of A/B domains with active genes enrichment
           act_domA.grl <- myOverlaps("actTSS","A",act_tss_chrX.gr,compA_chrX.gr,500,F)
           act_domB.grl <- myOverlaps("actTSS","B",act_tss_chrX.gr,compB_chrX.gr,500,F)
           
           a_mat1 <- length(act_domA.grl$gr1xgr2)
           b_mat1 <- length(act_domB.grl$gr1xgr2)
           a_mat2 <- length(act_domA.grl$gr1nogr2)
           b_mat2 <- length(act_domB.grl$gr1nogr2)
           f_mat <- matrix(c(tst_feat=a_mat1,tst_noFeat=a_mat2,ctrl_feat=b_mat1,ctrl_noFeat=b_mat2),nrow=2,ncol=2)
           
           ft <- fisher.test(f_mat)
           ft.pv <- ft$p.value
           ft.fc <- ft$estimate
           
           if(ft.fc < 1){
             compA_chrX.gr$comp <- "B"
             compA_chrX.gr$score <- compA_chrX.gr$score * -1
             compB_chrX.gr$comp <- "A"
             compB_chrX.gr$score <- compB_chrX.gr$score * -1
           }
           c(compA_chrX.gr,compB_chrX.gr)
         })

  
  new_COMP.gr <- sort(Reduce("c",new_COMP.grl))
  export.bw(new_COMP.gr,paste0(OUTC.d,NM,"_ab_comp_res",res,"_S2_dm3.gr.bw"))
}
# ----
