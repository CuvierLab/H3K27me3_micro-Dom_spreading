All Figure_X_mDOM.R scripts can be launch as such as long as you keep the organization of all directories as it is, especially the path to /DATA/.. .
Also, the code need rstudioapi package to run (base package under Rstudio) otherwise headers have to be modified in order to retrieve the correct pathways (see dirname(rstudioapi::getSourceEditorContext()$path))

All libraries and homemade functions needed are loaded in /LIBS/lrc_lib.R
All the data needed are loaded in /LIBS/data_lib.R. Some may have been preprocessed with scripts under /PROCESSING_SCRIPTS/.
/LIBS/mrc_lib.R is what has been used by the team to create a summary matrix called matrecap or mrc_tss.gr (GenomicRanges gene centered) that is used quite often along the manuscript.

I regret I did not had the time to comment all the functions used in details, so feel free to ask when some parts of code need highlights (a.heurteau@gmail.com or cuvier@ibcg.biotoul.fr)

[![DOI](https://zenodo.org/badge/259596873.svg)](https://zenodo.org/badge/latestdoi/259596873)
