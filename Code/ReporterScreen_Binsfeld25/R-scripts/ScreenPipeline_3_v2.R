
#============ Open output pdf ==============

file_id = paste0(Out_dir,"Screen_pipeline3.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)

#============ Fit to average EVC - fixed by Ana 05.11.22 ==============

par(mfrow=c(2,2))

reference_EVC = All_Lux_by_OD_AUC[grep("EVC",names(All_Lux_by_OD_AUC))]
reference_EVC = rowMeans(reference_EVC)

for(i in 1:length(All_Lux_by_OD_AUC))
{
  #===== plot Promo-I against EVC =====
  
  plot(reference_EVC, All_Lux_by_OD_AUC[[i]], ylab=paste("AUC LUX/OD", Plates[i]), xlab= paste("Average AUC LUX/OD EVC"), 
       col= bxp_clrs[i], main= names(All_Lux_by_OD_AUC[i]))
  
  #===== Calculate Pearson correlation =====
  
  r=cor(reference_EVC, All_Lux_by_OD_AUC[[i]], use = "complete.obs")
  #r=round(r,3)
  
  #===== Apply linear fit, plot fit, calculate RSS, TSS and Rsqr =====
  
  fit=lm(All_Lux_by_OD_AUC[[i]]~reference_EVC)
  abline(fit, col= "black", lwd=1)
  fit_RSS=sum(fit$residuals^2)
  
  #Get all non NA wells accounted for the fit
  fit_wells = row.names(All_Lux_by_OD_AUC)[as.numeric(names(fit$residuals))]
  fit_wells_data = All_Lux_by_OD_AUC[[i]][match(fit_wells,row.names(All_Lux_by_OD_AUC))]
  
  #calculate TSS and Rsqr
  TSS= sum((fit_wells_data-mean(fit_wells_data, na.rm = T))^2, na.rm = T)
  Rsqr= 1-fit_RSS/TSS
  
  #===== Apply robust huber regression, plot fit, calculate RSS, TSS and Rsqr =====
  
  fitr=rlm(All_Lux_by_OD_AUC[[i]]~reference_EVC)
  abline(fitr, col= "black", lty= 2, lwd=1)
  fitr_RSS = fitr$residuals^2
  fitr_RSS = fitr$w*fitr_RSS
  fitr_RSS = sum(fitr_RSS)
  Rsqrr= 1-fitr_RSS/TSS
  
  #===== Trim top and lowest 10% and plot, calculate RSS, TSS and Rsqr =====
  
  cutoff_high = quantile(fitr$residuals,na.rm=T,probs = 0.90)
  cutoff_low = quantile(fitr$residuals,na.rm=T,probs = 0.1)
  
  cut_names_high = as.numeric(names(fitr$fitted.values[fitr$residuals<cutoff_high]))
  cut_names_low = as.numeric(names(fitr$fitted.values[fitr$residuals>cutoff_low]))
  cut_names= as.numeric(cut_names_high[cut_names_high %in% cut_names_low])
  No_outlrs_wells = row.names(All_Lux_by_OD_AUC)[cut_names]
  
  points(reference_EVC[match(No_outlrs_wells,row.names(All_Lux_by_OD_AUC))],
         All_Lux_by_OD_AUC[[i]][match(No_outlrs_wells,row.names(All_Lux_by_OD_AUC))],
         col= bxp_clrs[i], pch=19)
  
  #Recalculate RSS, TSS & Rsqr of fitr, excluding outliers
  fitr_w = fitr$w[match(cut_names,names(fitr$residuals))]
  fitr_r = fitr$residuals[match(cut_names,names(fitr$residuals))]
  RSS_cut = fitr_r^2
  RSS_cut = fitr_w*RSS_cut
  RSS_cut = sum(RSS_cut)
  Rsqrr_cut=1-RSS_cut/TSS
  
  legend("topleft", legend = c(paste0("Pearson = ", round(r,2)), paste0("lin_reg,"," Rsqr = ", round(Rsqr,3)), 
                               paste0("hub_reg,"," Rsqr = ", round(Rsqrr,3),", Rsqr_cut = ", round(Rsqrr_cut,3))), bty= "n", col = c("white", "black", "black"), 
         lty = c(0,1,2), cex = 0.6)
  
  
  #===== Report data ======
  
  #Fitting parameters
  fit_p = as.data.frame(c(fitr$coefficients[[2]],fitr$coefficients[[1]]))
  row.names(fit_p) = c("slope","intercept")
  names(fit_p) = names(All_Lux_by_OD_AUC)[[i]]
  if(i==1)
  {Fitted_parameters=fit_p} else
  {Fitted_parameters=as.data.frame(cbind(Fitted_parameters,fit_p))}
  
  #Pearson & Rsqrs
  fit_pr = as.data.frame(c(r,Rsqr,Rsqrr,Rsqrr_cut))
  row.names(fit_pr) = c("Pearson_cor","Rsqr_lin_reg","Rsqr_rob_reg","Rsqr_rob_reg_cut")
  names(fit_pr) = names(All_Lux_by_OD_AUC)[[i]]
  if(i==1)
  {Fitted_PrsnRsqr=fit_pr} else
  {Fitted_PrsnRsqr=as.data.frame(cbind(Fitted_PrsnRsqr,fit_pr))}
  
  #Residuals
  scores = rep(NA,length(All_Lux_by_OD_AUC[[1]]))
  scores[match(fit_wells,row.names(All_Lux_by_OD_AUC))] = as.numeric(fitr$residuals)
  scores = as.data.frame(scores)
  row.names(scores) = row.names(All_Lux_by_OD_AUC)
  names(scores) = names(All_Lux_by_OD_AUC)[i]
  if(i==1)
  {Scores=scores} else
  {Scores = as.data.frame(cbind(Scores,scores))}
  
  #fit.w
  weights = rep(NA,length(All_Lux_by_OD_AUC[[1]]))
  weights[match(fit_wells,row.names(All_Lux_by_OD_AUC))] = as.numeric(fitr$w)
  weights = as.data.frame(weights)
  row.names(weights) = row.names(All_Lux_by_OD_AUC)
  names(weights) = names(All_Lux_by_OD_AUC)[i]
  if(i==1)
  {Weights=weights} else
  {Weights = as.data.frame(cbind(Weights,weights))}
} 

Results_list2 = list(Scores, Weights, Fitted_parameters, Fitted_PrsnRsqr)
names(Results_list2) = c("Scores","Weights", "Fitted_parameters", "Fitted_PrsnRsqr")
Results_list = append(Results_list,Results_list2)
rm(Results_list2)

#===== Rearrange Scores to be grouped by drug and concentration =====

new_names = paste0(Map$Drug[match(row.names(Scores),Map$Well)],
                   Map$ConcMock[match(row.names(Scores),Map$Well)])
row.names(Scores) = new_names
new_order = row.names(Scores)[order(row.names(Scores))]
Scores = Scores[match(new_order,row.names(Scores)),]

Results_list$Scores = Scores

#============== Summary plot fitting quality ============

par(mfrow=c(1,1))
here_colors = colorRampPalette(c("grey",Prom_clrs[3]))(4)
barplot(as.matrix(Results_list$Fitted_PrsnRsqr), beside = T,las=2,col=here_colors,
        main="Overview fit quality", ylab="Correlation coefficient or Coefficient of determination (R2)")
legend("topright",legend = row.names(Results_list$Fitted_PrsnRsqr),fill=here_colors,bty="n")
rm(here_colors)

#============== Plot Scores vs growth for checking relationship of promoter activation and growth. Ana added on 06.11.22 =======


hits_cutoff = 2*sapply(Scores,FUN=sd,na.rm=T)

#for (i in 1:length(Scores))
#{
#  hits = row.names(Scores)[grep(T,Scores[[i]]>hits_cutoff[[i]])]
#  
#  if(i==1)
#  {Hits=list()}
#  Hits[[length(Hits)+1]] = hits
#  
#}
#names(Hits) = names(Scores)
#
#Results_list[[length(Results_list)+1]] = Hits
#names(Results_list)[length(Results_list)] = "Hits"
#
#Results_list[[length(Results_list)+1]] = hits_cutoff
#names(Results_list)[length(Results_list)] = "hits_cutoff"


#===== Assign hits based on stdev: Alternative 2 ======

#alternative hits_cutoff by pooling replicates
for(i in 1:length(Promoters))
{
  promoter = Promoters[i]
  promoter_scores = Scores[grep(promoter,names(Scores))]
  promoter_sd = sd(as.matrix(promoter_scores),na.rm=T)
  
  hits_cutoff[grep(promoter,names(hits_cutoff))] = 2*promoter_sd
}

#Assign hits
for(i in 1:length(Scores))
{
  hits = row.names(Scores)[grep(T,Scores[[i]]>hits_cutoff[[i]])]
  
  if(i==1)
  {Hits=list()}
  Hits[[length(Hits)+1]] = hits
  
}
names(Hits) = names(Scores)

Results_list[[length(Results_list)+1]] = Hits
names(Results_list)[length(Results_list)] = "Hits2"

Results_list[[length(Results_list)+1]] = hits_cutoff
names(Results_list)[length(Results_list)] = "hits_cutoff2"

#===== plot Scores & hits per plate alternative 2 =====

par(mfrow=c(2,1))
for(i in 1:length(Promoters))
{
  promoter = Promoters[i]
  promoter_color = Prom_clrs[i]
  promoter_scores = Scores[grep(promoter,names(Scores))]
  promoter_cutoffs = Results_list$hits_cutoff2[grep(promoter,names(Results_list$hits_cutoff2))]
  promoter_hits = Results_list$Hits2[grep(promoter,names(Results_list$Hits2))]
  
  Hits = Results_list$Hits2
  hits_cutoff = Results_list$hits_cutoff2
  
  y_up = max(as.matrix(promoter_scores),na.rm=T)*1.1
  y_down = min(as.matrix(promoter_scores),na.rm=T)*1.1
  for(j in 1:length(promoter_scores))
  {
    if(j==1)
    {
      plot(seq(1,length(promoter_scores[[j]])),promoter_scores[[j]],
           ylab="Scores",xlab="Drugs", main = promoter,cex=0.5,
           col=promoter_color, ylim=c(y_down,y_up))
      points(c(-10,1000),c(promoter_cutoffs[[j]],promoter_cutoffs[[j]]),type="l",col="black",lwd=2,lty=2)
      points(c(-10,1000),c(median(promoter_scores[[j]],na.rm=T),median(promoter_scores[[j]],na.rm=T)),type="l",col="black",lwd=2,lty=2)
    } else
    {
      points(seq(1,length(promoter_scores[[j]])),promoter_scores[[j]],pch=19,
             col=promoter_color,cex=0.5)
      points(c(-10,1000),c(promoter_cutoffs[[j]],promoter_cutoffs[[j]]),type="l",col="black",lwd=2)
    }
    
    hits_data = promoter_scores[match(promoter_hits[[j]],row.names(promoter_scores)),j]
    #points(match(promoter_hits[[j]],row.names(promoter_scores)),hits_data,col="grey",pch=19)
    text(match(promoter_hits[[j]],row.names(promoter_scores)),(hits_data+promoter_cutoffs[[j]]*0.1),
         labels = promoter_hits[[j]],cex=0.5)
    
  }
  
}

#========= Calculate Z-scores from individual plates ===========

mean_vector = sapply((Scores),FUN=median,na.rm=T)
sd_vector = sapply((Scores),FUN=sd,na.rm=T)
Z_scores = as.data.frame(t(as.data.frame((t(Scores)-mean_vector)/sd_vector)))
boxplot(Scores,las=2, main="Non-scaled Scores",col=bxp_clrs, ylab="Scores")
boxplot(Z_scores,las=2, main="Scores scaled by Z-score",col=bxp_clrs,ylab="Z-scores")

Results_list[[length(Results_list)+1]] = Z_scores
names(Results_list)[length(Results_list)] = "Z_scores"

#===== plot Heatmap for Z-scores =====
na_rows = lapply(lapply(as.data.frame(t(Scores)),FUN=is.na),FUN=sum)
na_rows = na_rows[grep(T,na_rows >= (length(Scores))-4)]
heatmap_data = Z_scores[-match(names(na_rows),row.names(Scores)),]
heat_colors = colorRampPalette(c(ColorScheme[5,1],"white",ColorScheme[5,4]))

#remove EVC
heatmap_data = heatmap_data[-grep("EVC",names(heatmap_data))]

par(mfrow=c(1,1))
heatmap.2(as.matrix(heatmap_data),trace="none",na.rm = T,margins = c(6,8),col = heat_colors,
          cexRow = 0.1, Colv = T, Rowv = T,  main="Scaled by Z-score", symbreaks = T)

heatmap.2(cov(heatmap_data,use = "complete.obs"),trace="none",na.rm = T,margins = c(6,8),col = heat_colors,
          Colv = T, Rowv = T,  main="Cov matrix,Scaled by Z-score", symbreaks = T)

for(p in 2:length(Promoters))
{
  promoter = Promoters[p]
  promoter_scores = heatmap_data[grep(promoter,names(heatmap_data))]
  promoter_scores = c(promoter_scores[[1]], promoter_scores[[2]])
  
  if(p==2)
  {heatmap_data2 = data.frame(promoter_scores)}
  else
  {heatmap_data2 [p-1] = data.frame(cbind(promoter_scores))}
}
names(heatmap_data2) = Promoters[2:8]
heatmap_drugs = rep(row.names(heatmap_data),2)

heatmap.2(cov(heatmap_data2,use = "complete.obs"),trace="none",na.rm = T,margins = c(6,8),col = heat_colors,
          Colv = T, Rowv = T,  main="Cov matrix,Scaled by Z-score", symbreaks = T)

# ======== plot Z-scores vs growth =======
All_OD_AUC_drugs = paste0(Map$Drug[match(row.names(All_OD_AUC),Map$Well)],Map$ConcMock[match(row.names(All_OD_AUC),Map$Well)])

par(mfrow=c(2,2))
for(p in 1:length(Promoters))
{
promoter = Promoters[p]
promoter_color = Prom_clrs[p]
promoter_scores = Z_scores[grep(promoter,names(Z_scores))]
promoter_scores = promoter_scores[match(All_OD_AUC_drugs,row.names(promoter_scores)),]
promoter_scores = c(promoter_scores[[1]], promoter_scores[[2]])
promoter_OD = All_OD_AUC[grep(promoter,names(All_OD_AUC))]
promoter_OD = c(promoter_OD[[1]], promoter_OD[[2]])

plot(promoter_OD,promoter_scores,
       xlab="All_OD_AUC",ylab="Z-scores", xlim = c(0.1,3.9),
       pch=1,col=promoter_color, cex=0.6, #ylim=c(-2000,10000),
       main=promoter)
  x =rlm(promoter_scores ~ promoter_OD)
  abline(x)
  y = cor.test(promoter_OD, promoter_scores, use = "complete.obs", method = "pearson")
  legend("topright",legend=c(paste0("Rfit slope = ",round(x$coefficients[2],3)), paste0("Pearson cor = ", round(y$estimate,3)),
                             if(y$p.value<0.0001)
                             { paste0("p-val = ", formatC(y$p.value, format = "e", digits = 2)) }
                             else
                             { paste0("p-val = ", round(y$p.value, 3)) } 
                             ),bty="n")
  
  if(p==1)
  {Rfitslopes = round(x$coefficients[2],3)} 
  else
  {Rfitslopes = as.data.frame(rbind(Rfitslopes,round(x$coefficients[2],3)))}
  
  if(p==1)
  {cors = round(y$estimate,3)} 
  else
  {cors = as.data.frame(rbind(cors,round(y$estimate,3)))}
  
}

colnames(cors) = "Pearson correlation"
rownames(cors) = Promoters
cors[,2]=rownames(cors)
Pcors = cors[order(cors$"Pearson correlation", decreasing = F),]
Pcors = Pcors[1]
par(mfrow=c(1,1))
barplot(Pcors[[1]], names.arg = rownames(Pcors), las=2, main = "Pearson correlation", ylab = "Pearson correlation", col= Prom_clrs, ylim = c(-0.8,0.4))

file_id = paste0(Out_dir,"Pcors.txt")
write.table(Pcors,file = file_id,sep="\t",row.names = T)

#==== Close pdf & save environment Clean up ====
dev.off()

Results_list = Results_list[-grep("Hits",names(Results_list))]
Results_list = Results_list[-grep("hits",names(Results_list))]

#Clean up
rm(Hits,hits_cutoff,hits_data, hits,j,i,p,r,x,y,z,
   promoter,promoter_color,promoter_cutoffs,promoter_hits,promoter_sd,promoter_scores,
   fit,fit_p,fit_pr,fitr,fitr_r,fitr_RSS,fitr_w,fit_RSS,
   na_rows, scores, weights,Weights, Fitted_PrsnRsqr,
   cut_names,cut_names_high,cut_names_low,cutoff_high,cutoff_low,file_id,fit_wells,fit_wells_data,
   mean_vector,new_names,new_order,No_outlrs_wells,reference_EVC,Rsqr,Rsqrr,Rsqrr_cut,RSS_cut,sd_vector,TSS,y_down,y_up, 
   All_OD_AUC_drugs, promoter_OD, Rfitslopes, pAM, pmarRAB_names, pmarRAB_points, pmarRAB_scores, 
   pacrAB_names, pacrAB_points, pacrAB_scores, points_a, points_b, points_a_names, points_b_names, nlsdf, promoter_points, Pcors, cors
)

save.image(paste0(Out_dir,"ScreenPipeline_3.RData"))


#Anas version of lines 321-327
#plot(All_OD_AUC[[p]],Scores[match(All_OD_AUC_drugs,row.names(Scores)),p],
#     xlab="All_OD_AUC",ylab="Residuals",
#     pch=1,col=bxp_clrs[p], cex=0.6, #ylim=c(-2000,10000),
#     main=names(All_OD_AUC)[p])
#x =rlm(Scores[match(All_OD_AUC_drugs,row.names(Scores)),p] ~ All_OD_AUC[[p]])
#abline(x)
#legend("topright",legend=paste0("Rfit slope=",round(x$coefficients[2],3)),bty="n")







