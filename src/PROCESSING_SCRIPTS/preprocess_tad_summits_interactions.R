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
source(paste0(OUT.d,"LIBS/mrc_lib.R"))
source(paste0(OUT.d,"/LIBS/data_lib.R")) # load data
TMP.d <- create(paste0(SPATH,"/TMP/"))


# PREPROCESS HIC WT----------------------------------------------------------
bin41=41
hicName <- "ramirezS2WT_1000"
cons="wang_tad"
mybin=bin41
hic_res=1000
name=paste0("K27_BORDER_",mybin,"_",hic_res,"_",hicName,"_",cons)

m.dir <- create(paste0(OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/"))
o.dir <- create(paste0(m.dir,"toDump/constraint_",cons,"/msize_",get("mybin"),"/res_",hic_res,"/"))


# Set a TAD as active or inactive as he overlaps or not with a K27 domain
seqlevelsStyle(wang.tad.gr) <- "Ensembl"
mydom.grl <- split(wang.tad.gr,seqnames(wang.tad.gr))
wang.tad.noXtrm.dm3.gr <- Reduce("c",(lapply(mydom.grl,function(X){X[-c(1,length(X))]})))

folp <- findOverlapPairs(wang.tad.noXtrm.dm3.gr,Het.gr)
wang.tad.K27.gr <- subsetByOverlaps(wang.tad.noXtrm.dm3.gr,Het.gr)
wang.tad.noK27.gr <- subsetByOverlaps(wang.tad.noXtrm.dm3.gr,Het.gr,invert=T)

pint <- pintersect(folp)
folp.df <- as.data.frame(folp)
folp.df$inter <- width(pint)
folp.df$nme <- paste0(folp.df$first.start,"_",folp.df$first.end)
dens.df <- folp.df %>% dplyr::select(nme,inter) %>% group_by(nme) %>% summarise(tot=sum(inter))

wang.tad.K27.gr$nm <- paste0(start(wang.tad.K27.gr),"_",end(wang.tad.K27.gr))

wang.tad.K27.gr$K27wdth <- 0 
wang.tad.K27.gr[match(dens.df$nme,wang.tad.K27.gr$nm)]$K27wdth <- dens.df$tot 
wang.tad.K27.gr$K27dens <- round(wang.tad.K27.gr$K27wdth*100/width(wang.tad.K27.gr))





wang.tad.noK27.gr$nm <- paste0(start(wang.tad.noK27.gr),"_",end(wang.tad.noK27.gr))
wang.tad.noK27.gr$K27wdth <- 0
wang.tad.noK27.gr$K27dens <- 0

wang.tad.gr <- sort(c(wang.tad.noK27.gr,wang.tad.K27.gr))

sapply(1:length(wang.tad.gr),function(i){
  cur.gr <- wang.tad.gr[i]
  # print(i)
  startBIN <- floor((start(cur.gr))/hic_res)
  endBIN <- floor((end(cur.gr))/hic_res)
  chr <- seqnames(cur.gr)
  
  s1 <- (as.integer(startBIN)-((get("mybin")-1)/2))*hic_res 
  e1 <- (as.integer(startBIN)+((get("mybin")-1)/2))*hic_res
  s2 <- (as.integer(endBIN)-((get("mybin")-1)/2))*hic_res
  e2 <- (as.integer(endBIN)+((get("mybin")-1)/2))*hic_res
  
  resu <- cbind(as.character(chr),s1,e1,as.character(chr),s2,e2)
  write.table(resu,paste0(o.dir,"int=startWANGtad_anc=endWANGtad_dom=",i,"_K27LV=",cur.gr$K27dens,"_TADwdth=",width(cur.gr),".bed"), quote=FALSE, sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              append = F)
})




create(paste0(OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/config/"))
conf <- paste0(OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/config/apa.cfg")
if(file.exists(conf)) unlink(conf)
d_apa <- create(paste0(OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/Matrices/",hicName,"/constraint_",cons,"/msize_",get("mybin"),"/res_",hic_res,"/"))
write.table(paste0("#!/bin/bash"),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('mydir="',o.dir,'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
# Global list variable, contain .hic data (either http or local) that I have used so far.
# can add the one you need in lirary lrc_lib.R
write.table(paste0('myhic="',mydothiclist[[hicName]],'"'),
            file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('dumpdir="',d_apa,'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0("res=",hic_res),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0("bin=",get("mybin")),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('hicName="',hicName,'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('constraint="',cons,'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('xpname="',name,'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
write.table(paste0('tadpath="',OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/config/wangTAD.rds",'"'),file=conf,quote = F,row.names = F,col.names = F,append = T)
system(paste0("chmod +x ",conf))

saveRDS(wang.tad.gr,paste0(OUT.d,"/DATA/PROCESSED/Create_APA_TAD_Domain_summits/",name,"/config/wangTAD.rds"))





#   -----------------------------------------------------------------------
