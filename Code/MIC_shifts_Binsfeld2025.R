require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(corrplot)
require(gtools)
require(caTools)
library(ggplot2)
require(drc)

#============ Create accessory functions ===========
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
file_id = paste0(here_path,"MICShifts_clean_Raw_data.txt")
Raw_data = read.table (file_id, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#========= extract data for compatible format & create a plate database ==========

CB_ids = unique(Raw_data$CB_id)

# ========== Calculate growth per MIC =========
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
  {
    if(length(growth[[1]])==length(Growth[[1]]))
    {Growth = as.data.frame (cbind(Growth,growth))} else
    {
      add_growth = rep(NA,length(Growth[[1]]))
      add_growth[match(row.names(growth),row.names(Growth))] = growth[[1]]
      add_growth = as.data.frame(add_growth)
      names(add_growth) = names(growth)
      row.names(add_growth) = row.names(Growth)
      Growth = as.data.frame(cbind(Growth,add_growth))
    }
  
  }
}
Growth[Growth<0]=0

#========== Plot all MIC shifts ========

pdf(paste0(Out_dir,"MICShifts_Binsfeld2025.pdf"),useDingbats = F)

Strains = unique(Raw_data$strain)

par(mfrow=c(2,2))
for(b in 1:length(Strains))
{
  strain = Strains[b]
  strain_CBs = unique(Raw_data$CB_id[grep_exact(strain,Raw_data$strain)])
  
  cbs_info = Raw_data[match(strain_CBs,Raw_data$CB_id),c(1:8)]
  combos = paste0(drug1 = cbs_info$Drug.1,"_",drug2 = cbs_info$Drug.2)
  cbs_info = as.data.frame(cbind(cbs_info,combos))
  
  for(c in 1:length(unique(combos)))
  {
    combo = unique(combos)[c]
    drug1 = cbs_info$Drug.1[match(combo,cbs_info$combo)]
    drug2 = cbs_info$Drug.2[match(combo,cbs_info$combo)]
    
    cb_ids = cbs_info$CB_id[grep_exact(combo,cbs_info$combo)]
    cbs_data = Growth[match(cb_ids,names(Growth))]
    
    cbs_data = as.data.frame(cbind(rep(names(cbs_data),each=length(cbs_data[[1]])),
                                    paste0(rep(names(cbs_data),each=length(cbs_data[[1]])),"_",rep(row.names(cbs_data),length(cbs_data))),
                                         as.vector(as.matrix(cbs_data))))
    
    cbs_data = as.data.frame(cbind(cbs_data,Raw_data$Drug.1.conc...µg.ml.[match(cbs_data[[2]],paste0(Raw_data$CB_id,"_",Raw_data$well))],
                                   Raw_data$Drug.2.conc...µg.ml.[match(cbs_data[[2]],paste0(Raw_data$CB_id,"_",Raw_data$well))]))
    names(cbs_data) = c("cb_id","info","growth","conc_drug1","conc_drug2")
    
    #calculate fitness by diving by "no antibitic, drug 1" data (n=1), individually per cb (replicate)
    fitness = rep(NA,length(cbs_data[[1]]))
    for (i in 1:length(cb_ids))
    {
      cb_id = cb_ids[i]
      cb_data = cbs_data[grep_exact(cb_id,cbs_data$cb_id),]
      no_antib_growth = mean(as.numeric(cb_data$growth[cb_data$conc_drug1==0]),na.rm=T)
      fitness[grep_exact(cb_id,cbs_data$cb_id)] = as.numeric(cbs_data$growth[grep_exact(cb_id,cbs_data$cb_id)])/no_antib_growth
    }
    cbs_data = as.data.frame(cbind(cbs_data,fitness))
    
    #Add replicate and category caff+/- for plotting
    cbs_data = as.data.frame(cbind(cbs_data,cbs_info$replicate[match(cbs_data$cb_id,cbs_info$CB_id)]))
    
    category = rep("Caff-",length(cbs_data[[1]]))
    category[grep(T,cbs_data$conc_drug2>0)] = "Caff+"
    
    drug1_column = rep(drug1,length(cbs_data[[1]]))
    drug2_column = rep(drug2,length(cbs_data[[1]]))
    strain_column = rep(strain,length(cbs_data[[1]]))
    
    cbs_data=as.data.frame(cbind(cbs_data,category,drug1_column,drug2_column,strain_column))
    names(cbs_data) = c("cb_id","info","growth","conc_drug1","conc_drug2","fitness","replicate","category","drug1","drug2","strain")
    
    plot_pch = cbs_data$replicate
    plot_pch[cbs_data$replicate==1] = 19
    plot_pch[cbs_data$replicate==2] = 17
    plot_pch[cbs_data$replicate==3] = 15
    plot_pch[cbs_data$replicate==4] = 18
    cbs_data=as.data.frame(cbind(cbs_data,plot_pch))
    
    #Fit MICs and calculate IC50
    plot(cbs_data$conc_drug1[cbs_data$category=="Caff-"],cbs_data$fitness[cbs_data$category=="Caff-"],
         pch=cbs_data$plot_pch[cbs_data$category=="Caff-"],col="grey",ylim=c(0,1.2),xlim=c(0,max(cbs_data$conc_drug1,na.rm=T)+max(cbs_data$conc_drug1,na.rm=T)*0.2),
         main=paste0(strain),ylab="fitness",xlab=drug1)
    points(cbs_data$conc_drug1[cbs_data$category=="Caff+"],cbs_data$fitness[cbs_data$category=="Caff+"],
         pch=cbs_data$plot_pch[cbs_data$category=="Caff+"])
    
    mL <- drm(cbs_data$fitness[cbs_data$category=="Caff-"] ~ cbs_data$conc_drug1[cbs_data$category=="Caff-"], fct = L.3(), type = "continuous")
    new_concs = seq(0,max(cbs_data$conc_drug1,na.rm=T)+max(cbs_data$conc_drug1,na.rm=T)*0.3,by=0.0001)
    predict_data = predict(mL,newdata = as.data.frame(new_concs,by=0.1))
    IC50 = round(new_concs[which.min(abs(predict_data-0.5))],4)
    points(new_concs,predict_data,type="l",col="grey")
    abline(v=IC50,col="grey",lty=2)
    text(x=max(cbs_data$conc_drug1,na.rm=T)-max(cbs_data$conc_drug1,na.rm=T)/2,y=1.1,adj=0,paste0("IC50=",IC50),col="grey")
    
    mL <- drm(cbs_data$fitness[cbs_data$category=="Caff+"] ~ cbs_data$conc_drug1[cbs_data$category=="Caff+"], fct = L.3(), type = "continuous")
    new_concs = seq(0,max(cbs_data$conc_drug1,na.rm=T)+max(cbs_data$conc_drug1,na.rm=T)*0.3,by=0.0001)
    predict_data = predict(mL,newdata = as.data.frame(new_concs,by=0.1))
    IC50_caff = round(new_concs[which.min(abs(predict_data-0.5))],3)
    points(new_concs,predict_data,type="l")
    abline(v=IC50_caff,lty=2)
    text(x=max(cbs_data$conc_drug1,na.rm=T)-max(cbs_data$conc_drug1,na.rm=T)/2,y=1,adj=0,paste0("IC50caff=",IC50_caff))
    
    legend("topright",legend = paste0("rep ",seq(1,length(unique(plot_pch)))) ,pch=unique(plot_pch),)
    
    if(b==1 && c==1)
    {Source_data = cbs_data} else
    {Source_data = as.data.frame(rbind(Source_data,cbs_data))}
    
  }
}

dev.off()

#Write out source data file
file_id = paste0(here_path,"MICShifts_Source_data.txt")
write.table(Source_data,file = file_id,quote = F,row.names = F,sep="\t")

