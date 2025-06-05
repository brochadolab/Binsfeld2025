#List of functions defined here


#======================== Basics ====================

#This function takes an object (x) and a vector (v) and returns the first position of a word contain x in v, when x is found in v. 
grep_uniq <- function(x,v)
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
  return(uniq_match[1])
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

#Gets quadrant "quad" from the dataset. Dataset = 384 well plate data. quad = 1,2,3 or 4
quadr <- function(quad,dataset)  #NEW
{
  all_wells = names(dataset)[2:length(dataset)]
  sel_wells = c()

  nr_rows = 16
  nr_cols = length(all_wells)/nr_rows
  c = 1
  for(r in 1:nr_rows)
  {
    while(c <= (nr_cols*r))
    {
      if(quad == 1 && odd(r) && odd(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 2 && odd(r) && even(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 3 && even(r) && odd(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      if(quad == 4 && even(r) && even(c))
      {sel_wells = c(sel_wells,all_wells[c])}
      
      c = c+1
    }
      
  }
  
  red_dataset = Red_dataset(dataset,sel_wells)
  return(red_dataset)
}

Red_dataset <- function(dataset, sel_wells)
{
  red_dataset = as.vector(dataset[1])
  for (well in 1:length(sel_wells))
  {
    well_name <- sel_wells[well]
    well_index <- grep_uniq(well_name,names(dataset))
    red_dataset <- cbind(red_dataset,dataset[well_index])
  }
  return (red_dataset)
}  

#Function returns a vector of 2 integers, ready to be the dimentions of a par
dimPar <- function(x,y)
{
  #x is the number of plots to be made
  #y is the number of the columns of the par
  if(x <= y)
  {
    nr_cols = x
    nr_rows = 1
  } else
  {
    if(x > y)
    {
      nr_cols = y
      nr_rows = ceiling(x/y)
    }
  }
  
  if (x == 96)
  {
    nr_cols = 12
    nr_rows = 8
  }
  
  c(nr_rows,nr_cols)
}

panel.empty <- function(x,y)
{
  xx_coord = max(x)/2
  yy_coord = max(y)/2
  text(xx_coord, yy_coord, "white", col = "white")
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="complete"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#Check: http://www.math.ist.utl.pt/~calves/courses/integra/capiii32.html
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

NonLin_Rate <- function(dataset,min_growth=0.08,max_lag=9,growth_window=c(0,0.45),time_window=c(0,10),plotting=T,time_cutoff=0)
{
  xData = dataset[[1]]
  
  if(plotting==T)
  {par(mfrow=c(8,12),mar=c(1,1,1,1))}
  
  for(w in 2:length(dataset))
  {
    yData = dataset[[w]]
    if(plotting==T)
    {
      plot(xData,yData,ylim=growth_window,xlim=time_window,col="black",pch=19,cex=0.7,
           frame=T,axes=F,ylab="",xlab="")
      text(x=2,y=(growth_window[2]-0.05),labels=names(dataset)[w])
      points(c(time_cutoff,time_cutoff),growth_window,type="l",lty=2)
    } #if(plotting==T)
    
    #Check for growth
    if(max(yData[grep(F,xData > max_lag)]) <= min_growth)
    {parameters=c(0,0,NA,0,NA,NA,0,0)} else
    {
      nonlin_fit <- SummarizeGrowth(xData, yData)
      
      #nonlin_fit <- gcFitModel(xData, yData, gcID = "undefined",control=grofit.control(model.type=c("gompertz")))
      
      if(is.na(nonlin_fit$vals[[1]])==FALSE)
      {
        fitted_curve=predict(nonlin_fit$model)
        fitted_point=0
        fitted_area=0
        if(time_cutoff>0)
        {
          fitted_point = fitted_curve[length(grep(T,xData<time_cutoff))+1]
          fitted_area = nonlin_fit$vals$auc_e
        }
        
        if(plotting==T)
        {points(xData,fitted_curve,col="red3",type="l",lwd=2)}#ylim=growth_window,xlim=time_window)} #if(plotting==T)
        
        R = cor(yData,fitted_curve)
        R2 = summary(lm(yData~fitted_curve))$r.squared
        parameters=c(nonlin_fit$vals$r,nonlin_fit$vals$k,NA,nonlin_fit$vals$auc_l,
                     R2,R,fitted_point,fitted_area)
      } else
      {parameters=c(0,0,0,0,0,2,0,0)}
    }
    
    parameters = as.data.frame(parameters)
    row.names(parameters) = c("Growth_rate_h-1","plateauOD","lag-phase","AUCintegral","Rsq","PearsonCor","fitted_point","AUCtrap")
    names(parameters) = names(dataset)[w]
    if(w==2)
    {Parameters = as.data.frame(t(parameters))} else
    {Parameters = as.data.frame(rbind(Parameters,t(parameters)))}
  }
  
  return(Parameters)   
}

#======================== For Christoph ==================== 14.02.2020

plot_96wells <- function(quad_data,here_ylim,time_window,stat_phase_start,here_color="grey",Map="Map",y_lab="OD595",x_lab="Time (h)")
{
  par(mfrow=c(8,12),mar=c(0.5,0.5,0.5,0.5))
  for(w in 2:length(quad_data))
  {
    well = names(quad_data)[w]
    label = as.character(Map$Label[match(well,Map$Well)])
    
    #Map[match(well,Map$Well),]
    here_color = as.character(Map$Color_drug[match(well,Map$Well)])
    plot(quad_data[[1]],quad_data[[w]],ylim=here_ylim,xlim=time_window, ylab=y_lab,xlab=x_lab,
         col=here_color,pch=19,cex=0.7, type="b",axes = F,frame=T,lwd=2)
    points(c(stat_phase_start,stat_phase_start),c(-1,10),type="l",col="grey")
    text(x=2,y=(here_ylim[2]-here_ylim[2]*0.1),labels=label,cex=0.7,adj=0)
  }
}

BiotekReform <- function(file_id) #30.06.2020
{
  file_content = readLines(file_id)
  
  results = list()
  
  if(21 %in% grep(96,file_content))
  {plate_format = "mtp96"} else
  {plate_format = "mtp384"}
  
  time0_i = grep_exact("Results",file_content)[1] + 2
  
  if(length(grep("Read 2:595",file_content))>0)
  {
    time_series_i = grep("Read 2:595",file_content)[1]+2
    data = file_content[time_series_i]
    flag1=T
    flag2=T
    
    while(flag1==T)
    {
      time_series_i = time_series_i+1
      data = file_content[[time_series_i]][1]
      data = unlist(strsplit(data,split='\t'))
      time = data[1]
      time = round(24 * as.numeric(times(time)),3)
      
      data = data[-c(1,2)]
      data = as.numeric(gsub(",", ".", gsub("\\.", "", data)))
      
      data = c(time,data)
      
      if(flag2==T)
      {
        time_series = data
        flag2 = F
      } else
      {
        time_series = as.data.frame(rbind(time_series,data))
      }
      
      if(file_content[[time_series_i+1]][1]=="" || length(grep("0:00:00",file_content[[time_series_i+1]][1]))>0)
      {flag1=F}
      
    }
    rm(flag1,flag2)
    time_series[1] = time_series[1]-time_series[1,1]
    
    if(plate_format == "mtp384")
    {header = paste0(rep(LETTERS[1:16],each=24),seq(1:24))} else
    {header = paste0(rep(LETTERS[1:8],each=12),seq(1:12))}
    
    names(time_series) = c("Time",header)
    row.names(time_series) = NULL
    
    time_series_OD = time_series
    results[[length(results)+1]] = time_series_OD
    names(results)[length(results)] = "OD"
    rm(time_series)
    #write.table(time_series ,file = paste0(Out_dir,Plates[p],".txt"),sep="\t",quote =F,row.names = F)
  }
  
  
  if(length(grep("Read 3:Lum",file_content))>0)
  {
    
    #get times series for OD
    time_series_i = grep("Read 3:Lum",file_content)[1]+2
    flag1=T
    flag2=T
    while(flag1==T)
    {
      time_series_i = time_series_i+1
      data = file_content[[time_series_i]]
      data = unlist(strsplit(data,split='\t'))
      time = data[1]
      time = round(24 * as.numeric(times(time)),3)
      
      data = data[-c(1,2)]
      data = as.numeric(gsub(",", ".", gsub("\\.", "", data)))
      
      data = c(time,data)
      
      if(flag2==T)
      {
        time_series = as.data.frame(t(data))
        flag2 = F
      } else
      {
        if(length(data)==length(time_series[1,]))
          time_series = as.data.frame(rbind(time_series,data))
      }
      
      if(file_content[[time_series_i+1]][1]=="" || length(grep("0:00:00",file_content[[time_series_i+1]][1]))>0)
      {flag1=F}
      
    }
    rm(flag1,flag2)
    time_series[1] = time_series[1]-time_series[1,1]
    
    if(plate_format == "mtp384")
    {header = paste0(rep(LETTERS[1:16],each=24),seq(1:24))} else
    {header = paste0(rep(LETTERS[1:8],each=12),seq(1:12))}
    names(time_series) = c("Time",header)
    row.names(time_series) = NULL
    
    time_series_Lux = time_series
    results[[length(results)+1]] = time_series_Lux
    names(results)[length(results)] = "Lux"
    rm(time_series)
    #write.table(time_series ,file = paste0(Out_dir,Plates[p],"_Lux.txt"),sep="\t",quote =F,row.names = F)
    
  }
  
  return(results)
  
}


#By Ana i June2022: Took from Martin and updated to extract OD & Lux or Fluo!
Ext_SAMIplates <- function(file_path, output_dir, invert=T ,Lux=F,Fluo=F) #ARB 09.12.2021
{
  dataset <- read.table(file_path, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, skip = 4)
  
  #Remove NA columns in case of 96 well plates recorded as 384
  if(length(dataset)>120 && sum(is.na(dataset[length(dataset)])) == length(dataset[[1]]))
  {
    dataset = dataset[1:(3+96)]
  }
  
  header =c("Family","Plate","Time",paste0(rep(LETTERS[1:8],each=12),seq(1,12)))
  if(length(dataset)>120)
  {header=c("Family","Plate","Time",paste0(rep(LETTERS[1:16],each=24),seq(1,24)))}
  
  names(dataset) = header        
  
  family_vec <- as.numeric(as.vector(dataset$Family))
  nr_families <- max(family_vec)
  
  families_OD_list <- list()
  for (family in 1:nr_families)
  {
    family_OD_data = dataset[grep_exact(family,family_vec),]
    
    if(length(grep(T,duplicated(family_OD_data))>0))
    {family_OD_data = family_OD_data[grep(F,duplicated(family_OD_data)),]}
    
    families_OD_list[[length(families_OD_list)+1]] <- family_OD_data
  }
  
  #Covert Date and time columns to time-intervals in hours and write the files
  for(family in 1:nr_families)
  {
    family_data = families_OD_list[[family]]
    Date_time <- as.character(as.vector(family_data$Time))
    
    dtm <- strptime(Date_time, format = "%m/%d/%Y %H:%M:%S", tz = "CET")
    Time_h = vector(mode = "numeric", length = length(dtm))  
    for (i in 2:length(Time_h))
    {
      j = i-1
      diff = dtm[i]-dtm[j]
      Time_h[i] = Time_h[j] + as.numeric(diff, units="hours")
    }
    
    family_data$Time = Time_h
    PlateID = as.character(family_data$Plate[1])
    
    family_data = family_data[-c(1,2)]
    
    if(invert==T)
    {
      if(length(family_data)==385)
      {orgnl_order = paste0(rep(LETTERS[1:16], each=24),rep(seq(1,24),16))} else
      {orgnl_order = paste0(rep(LETTERS[1:8], each=12),rep(seq(1,12),8))}
      new_wells = rev(orgnl_order)
      names(family_data) = c(names(family_data)[1],new_wells)
      family_data = family_data[c(1,match(orgnl_order,names(family_data)))]
    }
    
    if(Lux==T)
    {
      family_data_OD = family_data[seq(1,length(family_data[[1]]),by=2),]
      family_data_Lux = family_data[seq(2,length(family_data[[1]]),by=2),]
      
      file_id = paste0(Out_dir,PlateID,"-OD.txt")
      write.table(family_data_OD,file = file_id, append = F, sep = "\t",row.names = F,quote = F)
      file_id = paste0(Out_dir,PlateID,"-Lux.txt")
      write.table(family_data_Lux,file = file_id, append = F, sep = "\t",row.names = F,quote = F)
    }
    
    if(Fluo==T)
    {
      family_data_OD = family_data[seq(1,length(family_data[[1]]),by=2),]
      family_data_Fluo = family_data[seq(2,length(family_data[[1]]),by=2),]
      
      file_id = paste0(Out_dir,PlateID,"-OD.txt")
      write.table(family_data_OD,file = file_id, append = F, sep = "\t",row.names = F,quote = F)
      file_id = paste0(Out_dir,PlateID,"-Fluo.txt")
      write.table(family_data_Fluo,file = file_id, append = F, sep = "\t",row.names = F,quote = F)
    }
    
    if(Lux==F && Fluo==F)
    {
      family_data_OD = family_data
      
      file_id = paste0(Out_dir,PlateID,"-OD.txt")
      write.table(family_data_OD,file = file_id, append = F, sep = "\t",row.names = F,quote = F)
    }
    
  }
}


#=========== New =========
my_ggplot_converter <- function(a,add_classifier="") #added 21.11.22 by Ana
{
  x = as.vector(as.matrix(a))
  b = rep(names(a),each = length(a[[1]]))
  c = rep(row.names(a),length(a))
  
  if(add_classifier!="")
  {x = as.data.frame(cbind(x,b,c,rep(add_classifier,length(x[[1]]))))} else
  {x = as.data.frame(cbind(x,b,c))}
  return(x)
}