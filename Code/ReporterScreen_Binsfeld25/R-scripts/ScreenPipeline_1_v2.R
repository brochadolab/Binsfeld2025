
#============ Set directories ==============

here_path_strain = paste0(here_path,strain,"/")

#create Results_list
Results_list = list()

#============ Re-set directories to load plates and run pipeline ==============

Load_dir = paste0(here_path,"R-scripts/") 
Out_dir = paste0(here_path_strain,"Out/")

clr= Prom_clrs

#PlateDatabase
PlateDatabase = MasterPlateDatabase[grep_exact(strain,MasterPlateDatabase$Strain),]
Plates = PlateDatabase$Plate_Nr

#====== Read all OD and lux data into their dedicated lists ======
file_id = paste0(Load_dir,"All_OD_data.txt",collapse=NULL)
All_OD_data = read.table(file_id,header = T, sep = "\t")

file_id = paste0(Load_dir,"All_Lux_data.txt",collapse=NULL)
All_Lux_data = read.table(file_id,header = T, sep = "\t")

strain_Lux_data = All_Lux_data[grep(strain,All_Lux_data$Plate),]
strain_OD_data = All_OD_data[grep(strain,All_Lux_data$Plate),]

Lux_data = list()
OD_data = list()

for(i in 1:length(PlateDatabase$Plate_Nr))
{
  Lux_data[[i]] = subset(strain_Lux_data,strain_Lux_data$Plate == PlateDatabase$Plate_Nr[i])[,-1]
  OD_data[[i]] = subset(strain_OD_data,strain_OD_data$Plate == PlateDatabase$Plate_Nr[i])[,-1]
}

names(Lux_data) = PlateDatabase$Plate_Nr
names(OD_data) = PlateDatabase$Plate_Nr

#=================== Plot ===============
Out_dir = paste0(here_path_strain,"Out/")

#====== Do plots automatically for OD of all plates ============
uniqdrugs= as.character(unique(Map$Drug))

file_id = paste0(Out_dir,"Screen_pipeline1_ODraw_all_plates.pdf",collapse=NULL)
#Changed the y-axis limits to 0-1.2 instead of 0-1 - Ana 170622
pdf(file_id, useDingbats = F)
for(p in 1:length(OD_data))
{
  OD_plate = OD_data[[p]] 
  Plate_id = names(OD_data)[p]
  
  par(mfrow=c(8,12),mar=c(1,0.5,1,0.5),cex=0.32)
  
  uniqdrugs= as.character(unique(Map$Drug))
  
  for(drug in uniqdrugs)
  {
    drug_wells = as.character(Map$Well[grep_exact(drug, Map$Drug)])
    drug_ConcMock = as.character(Map$ConcMock[match(drug_wells, Map$Well)])
    drug_data = OD_plate[match(drug_wells ,names(OD_plate))]
    
    plot(OD_plate[[1]],drug_data[[1]], ylim= c(0,1.2), col=clr[1], pch=19, main=drug,axes=F,frame=T)
    for(i in 2:length(drug_data)) 
    {points(OD_plate[[1]],drug_data[[i]], col= clr[i], pch=19)}
    legend("topleft", legend = paste(drug_wells, drug_ConcMock), col = clr, pch=19, bty="n")
  }
  text(x = 9,y = 0.1,labels = Plate_id,cex = 2)
}
dev.off()

#====== Do plots automatically for Lux of all plates ======

file_id = paste0(Out_dir,"Screen_pipeline1_Lux_all_plates.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)
for(p in 1:length(Lux_data))
{
  Lux_plate = Lux_data[[p]] 
  Plate_id = names(Lux_data)[p]
  
  uniqdrugs= as.character(unique(Map$Drug))
  
  par(mfrow=c(8,12),mar=c(1,0.5,1,0.5),cex=0.32)
  for(drug in uniqdrugs)
  {
    drug_wells = as.character(Map$Well[grep_exact(drug, Map$Drug)])
    drug_ConcMock = as.character(Map$ConcMock[match(drug_wells, Map$Well)])
    drug_data = Lux_plate[match(drug_wells ,names(Lux_plate))]
    
    plot(Lux_plate[[1]],drug_data[[1]], ylim= c(0,PlateDatabase$y_limits[p]), col=clr[1], pch=19, main=drug,axes=F,frame=T,type="b")
    for(i in 2:length(drug_data)) 
    {points(Lux_plate[[1]],drug_data[[i]], col= clr[i], pch=19,type="b")}
    if(time_cutoff_low == 1)
    {points(c(time_cutoff_low,time_cutoff_low),c(0,1000000),type="l")}
    points(c(time_cutoff_high,time_cutoff_high),c(0,1000000),type="l")
    legend("topleft", legend = paste(drug_wells, drug_ConcMock), col = clr, pch=19, bty="n")
  }
  text(x = 9,y = 100,labels = Plate_id,cex = 2)
}
dev.off()

#====== Background correction of OD for all plates plus plots ======

file_id = paste0(Out_dir,"Screen_pipeline1_ODcorrected_all_plates.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)

OD_data_corrected = list()

for(p in 1:length(OD_data))
{
  OD_plate = OD_data[[p]] 
  Plate_id = names(OD_data)[p]
  
  bg_data= OD_plate[4,]
  bg_mean= t(as.data.frame(sapply(bg_data,FUN=mean)))
  OD_plate_corrected = as.data.frame(t(as.data.frame(t(OD_plate)) - bg_mean))
  OD_plate_corrected[1] = OD_plate[1]
  OD_data_corrected[[length(OD_data_corrected)+1]] = OD_plate_corrected
  
  par(mfrow=c(8,12),mar=c(1,0.5,1,0.5),cex=0.32)
  
  uniqdrugs= as.character(unique(Map$Drug))
  
  for(drug in uniqdrugs)
  {
    drug_wells = as.character(Map$Well[grep_exact(drug, Map$Drug)])
    drug_ConcMock = as.character(Map$ConcMock[match(drug_wells, Map$Well)])
    drug_data = OD_plate_corrected[match(drug_wells ,names(OD_plate_corrected))]
    
    plot(OD_plate_corrected[[1]],drug_data[[1]], ylim= c(0,1), col=clr[1], pch=19, main=drug,axes=F,frame=T)
    for(i in 2:length(drug_data)) 
    {points(OD_plate_corrected[[1]],drug_data[[i]], col= clr[i], pch=19)}
    if(time_cutoff_low == 1)
      {points(c(time_cutoff_low,time_cutoff_low),c(0,1000000),type="l")}
    points(c(time_cutoff_high,time_cutoff_high),c(0,1000000),type="l")
    legend("topleft", legend = paste(drug_wells, drug_ConcMock), col = clr, pch=19, bty="n")
  }
  text(x = 9,y = 0.1,labels = Plate_id,cex = 2)
}
dev.off()

names(OD_data_corrected) = names(OD_data)

#====== Calculate AUC of corrected ODs for all plates ========

for(p in 1:length(OD_data_corrected))
{
  OD_plate_corrected = OD_data_corrected[[p]]
  OD_plate_corrected = OD_plate_corrected[OD_plate_corrected[[1]]<time_cutoff_high,]
  OD_plate_corrected = OD_plate_corrected[OD_plate_corrected[[1]]>=time_cutoff_low,]
  Plate_id = names(OD_data_corrected)[p]
  
  OD_AUC = sapply(OD_plate_corrected, FUN= trapz, x = OD_plate_corrected[[1]])
  OD_AUC = as.data.frame(OD_AUC[-c(1)])
  
  names(OD_AUC) = Plate_id
  
  if(p==1)
  {All_OD_AUC = OD_AUC} else
  {All_OD_AUC = as.data.frame(cbind(All_OD_AUC, OD_AUC))}
}

#====== Calculate AUC of Lux for all plates ==== 

for(p in 1:length(Lux_data))
{
  Lux_plate = Lux_data[[p]]
  Lux_plate = Lux_plate[Lux_plate[[1]]<time_cutoff_high,]
  Lux_plate = Lux_plate[Lux_plate[[1]]>=time_cutoff_low,]
  Plate_id = names(Lux_data)[p]

  Lux_AUC = sapply(Lux_plate, FUN= trapz, x = Lux_plate[[1]])
  Lux_AUC = as.data.frame(Lux_AUC[-c(1)])
  
  names(Lux_AUC) = Plate_id
  
  if(p==1)
  {All_Lux_AUC = Lux_AUC} else
  {All_Lux_AUC = as.data.frame(cbind(All_Lux_AUC, Lux_AUC))}
}

#====== Clean up and Save environment ==== 
#keep: uniqdrugs
rm(p,Plate_id,bg_data,bg_mean, clr,All_files,file_id,i,
   drug,drug_ConcMock,drug_data,drug_wells,plate_data,
   Lux_AUC,Lux_data,Lux_plate,
   OD_AUC,OD_data,OD_data_corrected,OD_plate,OD_plate_corrected
)

Results_list[[length(Results_list)+1]] = All_Lux_AUC
names(Results_list)[length(Results_list)] = "All_Lux_AUC"

Results_list[[length(Results_list)+1]] = All_OD_AUC
names(Results_list)[length(Results_list)] = "All_OD_AUC"

save.image(paste0(Out_dir,"ScreenPipeline_1.RData"))



