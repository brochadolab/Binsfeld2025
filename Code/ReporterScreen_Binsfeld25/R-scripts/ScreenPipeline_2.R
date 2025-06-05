
#============ Open output pdf ==============

file_id = paste0(Out_dir,"Screen_pipeline2.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)

#============ Remove non-growers ============
median_All_OD_AUC = sapply(All_OD_AUC,FUN=median)
dummy = as.data.frame(t(as.data.frame(t(All_OD_AUC))/median_All_OD_AUC))
All_OD_AUC[dummy<0.1] = NA
All_Lux_AUC[dummy<0.1] = NA

for(i in 1:length(All_OD_AUC))
{
  for(s in 1:sum(is.na(All_OD_AUC[i])))
  {
    wells = rownames(subset(All_OD_AUC[i], is.na(All_OD_AUC[i])))
    Background = rep(Plates[i],length(wells))
    dummy = cbind(Map[match(wells, Map$Well),c(2,3,5)],Background)
    
    
  }
  if(i == 1){removed_wells = dummy} else {removed_wells = rbind(removed_wells,dummy)}
}

file_id = paste0(Out_dir,"removed_wells.txt")
write.table(removed_wells,file = file_id,sep="\t",row.names = T)

rm(dummy, removed_wells, wells, Background)


bxp_clrs=rep(Prom_clrs, each=2)
my_pch= sample(c(15:19), length(Prom_clrs),replace = T)

# ========= Plot replicate correlation =========

#replicate correlation
rep_cors_OD = c()
rep_cors_Lux = c()
rep_cors_LuxbyOD = c()
rep_cors_promoters = c()
for(p in 1:length(Promoters))
{
  promoter = Promoters[p]
  promoter_OD = All_OD_AUC[grep(promoter,names(All_OD_AUC))]
  promoter_Lux = All_Lux_AUC[grep(promoter,names(All_Lux_AUC))]
  promoter_ratio = All_Lux_AUC/All_OD_AUC
  promoter_ratio = promoter_ratio[grep(promoter,names(All_Lux_AUC))]
  
  x = cor(promoter_OD,use="complete.obs")
  rep_cor = x[col(x)==(row(x)-1)]
  rep_cors_OD = c(rep_cors_OD,rep_cor)
  rep_cors_promoters = c(rep_cors_promoters,rep(promoter, length(rep_cor)))
  
  x = cor((promoter_Lux),use="complete.obs")
  rep_cor = x[col(x)==(row(x)-1)]
  rep_cors_Lux = c(rep_cors_Lux,rep_cor)
  
  x = cor((promoter_ratio),use="complete.obs")
  rep_cor = x[col(x)==(row(x)-1)]
  rep_cors_LuxbyOD = c(rep_cors_LuxbyOD,rep_cor)
}
rep_cors = as.data.frame(cbind(rep_cors_promoters,rep_cors_OD,rep_cors_Lux,rep_cors_LuxbyOD))
rep_cors$rep_cors_OD = as.numeric(as.matrix(rep_cors$rep_cors_OD))
rep_cors$rep_cors_Lux = as.numeric(as.matrix(rep_cors$rep_cors_Lux))
rep_cors$rep_cors_LuxbyOD = as.numeric(as.matrix(rep_cors$rep_cors_LuxbyOD))

par(mfrow=c(2,2))
boxplot(rep_cors[c(2,3,4)],ylim=c(0,1),xlim=c(0,3.5),pch=1,col=alpha("black", alpha=0), border = alpha("black", alpha=0),
         ylab="Pearson correlation coefficients", names=c("OD","LUX","LUX/OD"),
         main="Replicate correlation")

points(runif(8,0.9,1.1),rep_cors[[2]],lwd=1, pch=21, bg= Prom_clrs)
points(runif(8,1.9,2.1),rep_cors[[3]],lwd=1, pch=21, bg= Prom_clrs)
points(runif(8,2.9,3.1),rep_cors[[4]],lwd=1, pch=21, bg= Prom_clrs)

points(c(0.9,1.1),c(median(rep_cors$rep_cors_OD),median(rep_cors$rep_cors_OD)),type="l",lwd=2)
points(c(1.9,2.1),c(median(rep_cors$rep_cors_Lux),median(rep_cors$rep_cors_Lux)),type="l",lwd=2)
points(c(2.9,3.1),c(median(rep_cors$rep_cors_LuxbyOD),median(rep_cors$rep_cors_LuxbyOD)),type="l",lwd=2)

legend("topleft", legend = Promoters, pch = 19, col = Prom_clrs, bty = "n", cex = 0.75)

plot(All_Lux_AUC[[1]],All_Lux_AUC[[2]],col=Prom_clrs[5],main="EVC lux replicates", xlab = colnames(All_Lux_AUC[1]), ylab = colnames(All_Lux_AUC[2]))
plot(All_Lux_AUC[[1]],All_Lux_AUC[[2]],main="EVC lux replicates",col="white", xlab = colnames(All_Lux_AUC[1]), ylab = colnames(All_Lux_AUC[2]))
text(All_Lux_AUC[[1]],All_Lux_AUC[[2]],col=Prom_clrs[5],labels=row.names(All_Lux_AUC))

plot(All_Lux_AUC[[1]],All_Lux_AUC[[2]],main="EVC lux replicates, zoom low signal wells",col="white",ylim=c(1,2000),xlim=c(1,2000), xlab = colnames(All_Lux_AUC[1]), ylab = colnames(All_Lux_AUC[2]))
text(All_Lux_AUC[[1]],All_Lux_AUC[[2]],col=Prom_clrs[5],labels=row.names(All_Lux_AUC))

rm(x)
Results_list[[length(Results_list)+1]] = rep_cors
names(Results_list)[length(Results_list)] = "Replicate_correlation"

# ========= Normalize lux by OD =========

All_Lux_by_OD_AUC = All_Lux_AUC/All_OD_AUC

par(mfrow=c(1,1))
boxplot((All_Lux_by_OD_AUC),main="ratioLux/OD for all plates",las=2,border="grey",ylab= "AUC Lux/OD", col = bxp_clrs,pch = my_pch, cex = 1.5)


#OD

boxplot(All_OD_AUC,main="AUC OD for all plates",las=2,border="black",ylab= "AUC OD",col = bxp_clrs,pch = 21, cex = 1, bg=bxp_clrs, ylim=c(0,4.5))

sel_wells = as.character(as.matrix(Map$Well[grep("Water",Map$Drug)]))
for(p in 1:length(All_OD_AUC))
{
  wells_data = All_OD_AUC[match(sel_wells,row.names(All_OD_AUC)),p]
  x_data = rep(p,length(wells_data))
  points(x_data,wells_data,pch=19,col="black")
}
legend("topleft",legend = paste0("Water wells n=",length(sel_wells)),pch=19,col="black",bty="n")


#LUX

boxplot(log10(All_Lux_AUC),main="AUC Lux for all plates",las=2,border="black",ylab= "log10 (AUC Lux)",col = bxp_clrs,pch = 21, cex = 1, bg=bxp_clrs)

sel_wells = as.character(as.matrix(Map$Well[grep("Water",Map$Drug)]))
for(p in 1:length(All_Lux_AUC))
{
  wells_data = log10(All_Lux_AUC)[match(sel_wells,row.names(All_Lux_AUC)),p]
  x_data = rep(p,length(wells_data))
  points(x_data,wells_data,pch=19,col="black")
}
legend("topleft",legend = paste0("Water wells n=",length(sel_wells)),pch=19,col="black",bty="n")


#LUX/OD

boxplot(log10(All_Lux_by_OD_AUC),main="AUC Lux/OD for all plates",las=2,border="black",ylab= "log10(AUC Lux/OD)",col = bxp_clrs,pch = 21, cex = 1, bg=bxp_clrs)

sel_wells = as.character(as.matrix(Map$Well[grep("Water",Map$Drug)]))
for(p in 1:length(All_Lux_by_OD_AUC))
{
  wells_data = log10(All_Lux_by_OD_AUC)[match(sel_wells,row.names(All_Lux_by_OD_AUC)),p]
  x_data = rep(p,length(wells_data))
  points(x_data,wells_data,pch=19,col="black")
}
legend("topleft",legend = paste0("Water wells n=",length(sel_wells)),pch=19,col="black",bty="n")


#=========== Normalize data to median AUCLux/OD =============

median_Lux_by_OD = sapply(All_Lux_by_OD_AUC,FUN=median,na.rm=T)
normalized_Lux_by_OD = as.data.frame(t(as.data.frame(t(All_Lux_by_OD_AUC))/median_Lux_by_OD))

boxplot((normalized_Lux_by_OD),main="normalized_Lux_by_OD",las=2,border="grey",ylim=c(0,2),ylab= "median norm [AUC Lux/OD]",col = bxp_clrs,pch = my_pch, cex = 1.5)
boxplot((normalized_Lux_by_OD),main="normalized_Lux_by_OD",las=2,border="grey",ylab= "median norm [AUC Lux/OD]",col = bxp_clrs,pch = my_pch, cex = 1.5)

#=========== Plot heatmap =============

na_rows = lapply(lapply(as.data.frame(t(normalized_Lux_by_OD)),FUN=is.na),FUN=sum)
na_rows = na_rows[grep(T,na_rows >= (length(normalized_Lux_by_OD))-4)]

heatmap_data = normalized_Lux_by_OD[-match(names(na_rows),row.names(normalized_Lux_by_OD)),]
row.names(heatmap_data) = paste0(Map$Drug[match(row.names(heatmap_data),Map$Well)],"_",Map$ConcMock[match(row.names(heatmap_data),Map$Well)])

heat_colors = colorRampPalette(c("#F48D00","white","#0A5A9D"))
heatmap.2(as.matrix(log(heatmap_data)),trace="none",na.rm = T,margins = c(6,8),col = heat_colors, cexRow = 0.1)

dev.off()

#====== Clean up and Save environment ====
#keep: bxp_clrs,uniqdrugs
rm(na_rows,normalized_Lux_by_OD,prom,prom_cor,prom_data,i,ind, my_pch,p,
   med_All_Lux_by_OD_AUC,median_All_OD_AUC,median_Lux_by_OD,
   promoter,promoter_Lux,promoter_OD,promoter_ratio,bg_cor,
   rep_cor,rep_cors_Lux, rep_cors_LuxbyOD,rep_cors_OD,rep_cors_promoters,sel_wells,
   test_wells,test_wells_conc,test_wells_data,wells_data,x_data
)

save.image(paste0(Out_dir,"ScreenPipeline_2.RData"))
