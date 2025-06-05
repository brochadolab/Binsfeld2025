require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(corrplot)
require(gtools)
require(ggplot2)
require(corrplot)
require(gtools)
require(caTools)
library(chron)
require(fmsb)
require(beeswarm)
require(scales)
require(MASS)
require(ggfortify)
require(igraph)

#============ Set directories ==============

here_path = ""
source(paste0(here_path,"R-scripts/Aux_functions.R",collapse=NULL))

#Read Map & General PlateDatabase
file_id = paste0(here_path,"R-scripts/MasterPlateDatabase.txt",collapse=NULL)
MasterPlateDatabase = read.table(file_id,header = T, sep = "\t")

file_id = paste0(here_path,"R-scripts/Map.txt",collapse=NULL)
Map = read.table(file_id,header = T, sep = "\t")

file_id = paste0(here_path,"R-scripts/ColorScheme2.txt",collapse=NULL)
ColorScheme = read.table(file_id,header = T, sep="\t",skip = 0,comment.char = "*")

Strains = unique(MasterPlateDatabase$Strain)
Promoters = unique(MasterPlateDatabase$Promoter)

#https://paletton.com/#uid=7000u0kr4uOgSFulZxNuKpSy7ks
Prom_clrs = c("#D0471C", "#D0801C", "#195A85", "#139252")
Prom_clrs = as.vector(t(ColorScheme[1,]))


for(s in 1:length(Strains)){

#============ Choose the strain Do not change between different running parts of the pipeline ==============

strain = Strains[s]

#============ Run ScreenPipeline_1.R for the chosen strain ==============

### Use time cutoff low/high 0/8 for WT and 1/9 for deletion strains due to different inocculum size
if(strain == "WT")
  {
  time_cutoff_low = 0
  time_cutoff_high = 8
  } else 
  {
  time_cutoff_low = 1
  time_cutoff_high = 9 
  }

source(paste0(here_path,"R-scripts/ScreenPipeline_1_v2.R",collapse=NULL))

#Write output files
file_id = paste0(Out_dir,"All_OD_AUC.txt")
write.table(All_OD_AUC,file = file_id,sep="\t",row.names = T)

file_id = paste0(Out_dir,"All_Lux_AUC.txt")
write.table(All_Lux_AUC,file = file_id,sep="\t",row.names = T)

#============ Run ScreenPipeline_2.R for the chosen strain ==============

here_path_strain = paste0(here_path,strain,"/")
Out_dir = paste0(here_path_strain,"Out/")
load(paste0(Out_dir,"ScreenPipeline_1.RData"))

Prom_clrs= colorRampPalette(Prom_clrs)(length(Promoters))
source(paste0(here_path,"R-scripts/ScreenPipeline_2.R",collapse=NULL))

#Write output files
file_id = paste0(Out_dir,"ReplicateCorrelation.txt")
write.table(rep_cors,file = file_id,sep="\t",row.names = F)

#============ Run ScreenPipeline_3.R for the chosen strain ==============

here_path_strain = paste0(here_path,strain,"/")
Out_dir = paste0(here_path_strain,"Out/")
load(paste0(Out_dir,"ScreenPipeline_2.RData"))

Prom_clrs= colorRampPalette(Prom_clrs)(length(Promoters))

source(paste0(here_path,"R-scripts/ScreenPipeline_3_v2.R",collapse=NULL))

#Write output files
file_id = paste0(Out_dir,"Scores.txt")
write.table(Scores,file = file_id,sep="\t",row.names = T)

file_id = paste0(Out_dir,"Z_scores.txt")
write.table(Z_scores,file = file_id,sep="\t",row.names = T)

file_id = paste0(Out_dir,"Fit_Parameters.txt")
write.table(Results_list$Fitted_parameters,file = file_id,sep="\t",row.names = T)

file_id = paste0(Out_dir,"Fit_Rqs.txt")
write.table(Results_list$Fitted_PrsnRsqr,file = file_id,sep="\t",row.names = T)

#============ Run ScreenPipeline_4.R for the chosen strain ==============

here_path_strain = paste0(here_path,strain,"/")
Out_dir = paste0(here_path_strain,"Out/")
load(paste0(Out_dir,"ScreenPipeline_3.RData"))

#test different colors
Prom_clrs=unlist(c(ColorScheme[2,],ColorScheme[4,]))

Z_score_cutval = 1

source(paste0(here_path,"R-scripts/ScreenPipeline_4_v2_Ana.R",collapse=NULL))

file_id = paste0(Out_dir,"Cytoscape_network.txt")
write.table(Cytoscape_network,file = file_id,sep="\t",row.names = F,quote = F)
file_id = paste0(Out_dir,"Node_attributes.txt")
write.table(Node_attributes,file = file_id,sep="\t",row.names = F,quote = F)
file_id = paste0(Out_dir,"Means.txt")
write.table(Means,file = file_id,sep="\t",row.names = T,quote = F)
file_id = paste0(Out_dir,"P_vals.txt")
write.table(P_vals,file = file_id,sep="\t",row.names = T,quote = F)
file_id = paste0(Out_dir,"Means_cidal.txt")
write.table(Means_cidal,file = file_id,sep="\t",row.names = T,quote = F)
file_id = paste0(Out_dir,"Means_static.txt")
write.table(Means_static,file = file_id,sep="\t",row.names = T,quote = F)
file_id = paste0(Out_dir,"Min_growth.txt")
write.table(Min_growth,file = file_id,sep="\t",row.names = T,quote = F)

}


