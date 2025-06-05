require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(corrplot)
require(gtools)
require(caTools)
library(ggplot2)

#============ Create asessory functions ===========
get_conc_grad <- function(c,steps=7)
{
  step = c/steps
  v=c((steps+1):1)
  v=v*step
  v = c-v
  v = round(c(v[2:(steps+1)],c),5)
  return(v)
}

get_conc_gradB <- function(c)
{
  step = c/11
  v=c(12:1)
  v=v*step
  v = c-v
  v = round(c(v[2:12],c),5)
  return(v)
}

plot_enhanced_CB <- function(plot_data,Isoboles,drug_hrz,drug_vert,
                             plotting=T,CB_title="",contour_color="grey50",no_growth_limit=0.000005,
                             grad_color ="grey")
{
  #remove non-growing concs
  if(length(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit))>0)
  {plot_data = plot_data[1:(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit)[1])]}
  plot_data = as.data.frame(t(plot_data))
  if(length(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit))>0)
  {plot_data = plot_data[1:(grep(F,sapply(plot_data,FUN=sum,na.rm=T)>= no_growth_limit)[1])]}
  #plot_data = as.data.frame(t(plot_data))
  
  ctlns <- contourLines(x = as.numeric(row.names(plot_data)),
                        y = as.numeric(names(plot_data)),
                        z = as.matrix(plot_data), levels=Isoboles)
  
  grad_color = colorRampPalette (c("white",grad_color,grad_color))
  #plot
  if(plotting==1)
  {
    try(filled.contour(x = as.numeric(row.names(plot_data)),
                       y = as.numeric(names(plot_data)),
                       z = as.matrix(plot_data),
                       col=grad_color(35),
                       main=CB_title,
                       ylab = drug_vert,xlab=drug_hrz,
                       ylim = c(0,max(as.numeric(names(plot_data)))),
                       xlim = c(0,max(as.numeric(row.names(plot_data)))),
                       levels=seq(0,1.1,by=0.05),
                       plot.axes={axis(1); axis(2);
                         sapply(1:length(Isoboles), function(x) lines(ctlns[[x]][[2]], ctlns[[x]][[3]],
                                                                      lwd=2,col=contour_color))}
    ),silent = T)
  } #if(plotting==1)
  
  return(ctlns)
}

#Like grep, but returns only exact match
grep_exact <- function(x,v)
{
  match_vec = grep(x, v, fixed = TRUE)
  uniq_match = c()
  
  #To make sure that well_index points to the exact one and only one well
  if (length(match_vec) > 0)
  {
    for(i in 1:length(match_vec))
    {
      match_vec_temp = match_vec[i]
      if (nchar(match_vec_temp) > 0 && x == v[match_vec_temp])
      {uniq_match = c(uniq_match, match_vec_temp)}
    }
  }
  return(uniq_match)
}

Trape_rule <- function(dataset)
{
  areas = c()
  
  xData = as.numeric(as.vector(dataset[[1]]))
  nr_points = length(xData)
  h = (xData[length(xData)]-xData[1])/(nr_points-1)
  
  for(i in 2:length(dataset))
  {
    yData = as.numeric(as.vector(dataset[[i]]))
    
    t = (yData[1] + yData[length(yData)])/2
    for(j in 2:(length(yData)-1))
    {
      t = t + yData[j]
    }
    
    Tn = t*h
    areas = c(areas,Tn)
  }
  names(areas) = names(dataset)[2:length(dataset)]
  return(areas)
  
}

get_feature <- function(plate_ids,feature,Plate_database)
{
  All_features = names(Plate_database)
  features = c()
  if(length(grep_exact(feature,All_features))>0)
  {
    column = grep_exact(feature,All_features)
    for(p in 1:length(plate_ids))
    {
      plate_id = plate_ids[p]
      features = c(features,as.character(Plate_database[grep_exact(plate_id,Plate_database[[1]]),column]))
    }
  }
  
  return(features)
}


# ============= Redefine directories ==============
here_path = ""
Load_dir = paste0(here_path,"")
Out_dir = paste0(here_path,"")

growth_window = c(0.005,1.2)
Bg_correct = T #FALSE
AUC = T

#Read data
file_id = paste0(here_path,"CBs_clean_Raw_data.txt")
Raw_data = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#========= extract data for compatible format & create a plate database ==========

CB_ids = unique(Raw_data$CB_id)

strains = Raw_data$strain[match(CB_ids,Raw_data$CB_id)]
drugs1 = Raw_data$Drug.1[match(CB_ids,Raw_data$CB_id)]
drugs2 = Raw_data$Drug.2[match(CB_ids,Raw_data$CB_id)]
replicates = Raw_data$replicate[match(CB_ids,Raw_data$CB_id)]

# ========== Calculate growth per CB =========
for(c in 1:length(CB_ids))
{
  cb_id = CB_ids[c]
  cb_data = Raw_data[grep_exact(cb_id,Raw_data$CB_id),]
  
  cb_time = t(cb_data[1,grep("time.in",names(cb_data))])
  cb_od = t(cb_data[grep("OD",names(cb_data))])
  
  dataset = as.data.frame(cbind(cb_time,cb_od))
  names(dataset)[1] = "time"
  
  #correct for baseline OD of medium
  corr_dataset=as.data.frame(t(as.data.frame(t(dataset[2:length(dataset)]))-t(dataset[3,2:length(dataset)])))
  corr_dataset=as.data.frame(cbind(dataset[[1]],corr_dataset))
  names(corr_dataset)[1] = "Time"
  dataset = corr_dataset
  rm(corr_dataset)
  
  #Calculate AUC until 8h
  growth = Trape_rule(as.data.frame(dataset[1:(length(grep(T,dataset[[1]]<8))),]))
  growth = as.data.frame(growth)
  names(growth) = cb_id
  row.names(growth) = cb_data$well
  
  if(c==1)
  {Growth=growth} else
  {Growth = as.data.frame (cbind (Growth,growth))}
  
  
}
Growth[Growth<0]=0

# =========== Get data and plot checkerboard ========

for(c in 1:length(CB_ids))
{
  cb_id = CB_ids[c]
  cb_raw_data = Raw_data[grep_exact(cb_id,Raw_data$CB_id),]
  cb_data = Growth[c]
  
  strain = cb_raw_data$strain[1]
  drug1 = cb_raw_data$Drug.1[1]
  drug2 = cb_raw_data$Drug.2[1]
  
  cb_data = as.data.frame(cbind(cb_raw_data$Drug.1.conc...µg.ml.,
                                cb_raw_data$Drug.2.conc...µg.ml.,
                                cb_data))
  
  data = as.data.frame(t(matrix(cb_data[[3]],ncol=8,nrow=12)))
  drug1_conc = as.data.frame(t(matrix(cb_data[[1]],ncol=8,nrow=12)))
  drug2_conc = as.data.frame(t(matrix(cb_data[[2]],ncol=8,nrow=12)))
  
  names(data) = drug1_conc[1,]
  row.names(data) = drug2_conc[[1]]
  data = data [rev(seq(1,12))]
  #calculate fitness by dividing all by non-treated control
  data = data/data[8,1]
  
  #calculate expected fitness for bliss model
  exp_fit_data = data
  for(i in 2:length(exp_fit_data))
  {exp_fit_data[i] = data[1]*data[8,i]}
  
  #Keep data
  if(c==1)
  {
    CB_Growth = list()
    CB_ExpFit = list()
  }
  CB_Growth[[length(CB_Growth)+1]] = data
  CB_ExpFit[[length(CB_ExpFit)+1]] = exp_fit_data
}
names(CB_Growth) = CB_ids
names(CB_ExpFit) = CB_ids

# ========= Extract isobles =======

Isoboles= seq(0.4,0.6,by=0.05)
for(c in 1:length(CB_ids))
{
  cb_id = CB_ids[c]
  cb_raw_data = Raw_data[grep_exact(cb_id,Raw_data$CB_id),]
  cb_data = Growth[c]
  
  strain = cb_raw_data$strain[1]
  drug1 = cb_raw_data$Drug.1[1]
  drug2 = cb_raw_data$Drug.2[1]
  replicate = cb_raw_data$replicate[1]
  
  #get fitness
  data = CB_Growth[[grep_exact(cb_id,names(CB_Growth))]]
  
  #get expected fitness
  exp_fit_data = CB_ExpFit[[grep_exact(cb_id,names(CB_ExpFit))]]
  
  #Plot
  grad_color = colorRampPalette (c("white","chocolate1","chocolate1"))
  #plot fitness
  plot_data = data
  plot_data = plot_data[rev(1:length(plot_data[[1]])),]
  CB_title = paste0(strain," rep",replicate," Bliss interaction")
  interac_ctlns = plot_enhanced_CB(plot_data,Isoboles,drug1,drug2,plotting = F)
  
  #plot expected fitness
  plot_data = exp_fit_data
  plot_data = plot_data[rev(1:length(plot_data[[1]])),]
  blissModel_ctlns = plot_enhanced_CB(plot_data,Isoboles,drug1,2,plotting = F)
  rm(plot_data)
  
  #Store the data
  if(c==1)
  {
    All_blissModel_ctlns = c()
    All_interac_ctlns = c()
  }
  
  All_blissModel_ctlns[[length(All_blissModel_ctlns)+1]] = blissModel_ctlns
  All_interac_ctlns[[length(All_interac_ctlns)+1]] = interac_ctlns
}
names(All_blissModel_ctlns) = CB_ids
names(All_interac_ctlns) = CB_ids

#=========== Plot isoboles =============
pdf(paste0(here_path,"Checkerboards_Binsfeld2025.pdf"),useDingbats = F)
par(mfrow = c(2,2))
for(c in 1:length(CB_ids))
{
  cb_id = CB_ids[c]
  cb_raw_data = Raw_data[grep_exact(cb_id,Raw_data$CB_id),]
  cb_data = Growth[c]
  
  strain = cb_raw_data$strain[1]
  drug1 = cb_raw_data$Drug.1[1]
  drug2 = cb_raw_data$Drug.2[1]
  replicate = cb_raw_data$replicate[1]
  
  #Get a gradient for isobole colors
  Isob_cols = colorRampPalette (c("white","chocolate1","chocolate1"))(length(Isoboles)+1)
  Isob_cols = Isob_cols[-1]
  
  Isob_cols_bliss = colorRampPalette (c("white","grey50","grey50"))(length(Isoboles)+1)
  Isob_cols_bliss = Isob_cols_bliss[-1]
  
  #get fitness
  data = CB_Growth[[grep_exact(cb_id,names(CB_Growth))]]
  
  #get expected fitness
  exp_fit_data = CB_ExpFit[[grep_exact(cb_id,names(CB_ExpFit))]]
  
  #Plot
  grad_color = colorRampPalette (c("white","chocolate1","chocolate1"))
  #plot fitness
  plot_data = data
  plot_data = plot_data[rev(1:length(plot_data[[1]])),]
  max_plot_x = max(as.numeric(names(plot_data)))
  max_plot_y = max(as.numeric(row.names(plot_data)))
  
  CB_title = paste0(strain," rep",replicate)
  
  #Get the contours
  blissModel_ctlns = All_blissModel_ctlns[[c]]
  interac_ctlns = All_interac_ctlns[[c]]
  
  plot(blissModel_ctlns[[1]][[2]],blissModel_ctlns[[1]][[3]],
       main=CB_title,
       ylab = drug2,xlab=drug1,
       ylim = c(0,max_plot_y), xlim=c(0,max_plot_x),
       lwd=1,col=Isob_cols_bliss[1],type="l")
  try(sapply(1:length(Isoboles), function(x) lines(blissModel_ctlns[[x]][[2]], blissModel_ctlns[[x]][[3]],
                                                   lwd=1,col=Isob_cols_bliss[x])),silent = T)
  try(sapply(1:length(Isoboles), function(x) lines(interac_ctlns[[x]][[2]], interac_ctlns[[x]][[3]],
                                                   lwd=1,col=Isob_cols[x])),silent = T)
  legend("topright",legend = (Isoboles),col=Isob_cols,lwd=1,title = "Isoboles (fitness)",bty="n")
  
}

dev.off()
