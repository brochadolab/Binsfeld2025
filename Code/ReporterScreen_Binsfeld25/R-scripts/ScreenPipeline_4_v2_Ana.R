#============ Open output pdf ==============

file_id = paste0(Out_dir,"Screen_pipeline4.pdf",collapse=NULL)
pdf(file_id, useDingbats = F)

#========= Calculate p-vals based on a Wilcox test =================
Drugs = sort(unique(Map$Drug))
Scores = Z_scores

for(p in 1:length(Promoters))
{
  promoter = Promoters[p]
  promoter_scores = Scores[grep(promoter,names(Scores))]
  
  water_scores = promoter_scores[grep("Water",row.names(promoter_scores)),]
  
  p_vals = c()
  means = c()
  for(d in 1:length(Drugs))
  {
    drug = Drugs[d]
    drug_scores = promoter_scores[grep(drug,row.names(promoter_scores)),]
    p_val = wilcox.test(as.vector(as.matrix(drug_scores)),as.vector(as.matrix(water_scores)))
    p_val = p_val$p.value
    p_vals = c(p_vals,p_val)
    drug_mean = mean(as.vector(as.matrix(drug_scores)),na.rm=T)
    means = c(means,drug_mean)
  }
  
  p_vals = p.adjust(p_vals,method = "BH")
  p_vals = as.data.frame(p_vals)
  row.names(p_vals) = Drugs
  names(p_vals) = promoter
  
  means = as.data.frame(means)
  row.names(means) = Drugs
  names(means) = promoter
  
  if(p==1)
  {
    P_vals = p_vals
    Means = means
  } else
  {
    P_vals = as.data.frame(cbind(P_vals,p_vals))
    Means = as.data.frame(cbind(Means,means))
  }
}

#========= Sensititvity analysis on Z-score =================

cutval_sens = seq(0,6,by=0.05)
hits_sens = c()
for(c in 1:length(cutval_sens))
{
  cutval = cutval_sens[c]
  
  cutval_hits = 0
  for(i in 1:length(Promoters))
  {
    hits_p = c(row.names(P_vals)[grep(T,P_vals[i] < 0.05)])
    hits_m = c(row.names(Means)[grep(T,abs(Means[i]) > cutval)]) 
    hits = intersect(hits_p, hits_m)
    cutval_hits = cutval_hits+length(hits)
  }
  hits_sens = c(hits_sens,cutval_hits)
}
hit_rate_sens = hits_sens/(length(Drugs[-grep("Water",Drugs)])*length(Promoters[-grep("EVC",Promoters)])) 
plot(cutval_sens,hit_rate_sens,cex=0.6,col=Prom_clrs[1],pch=19,main="Sensitivity analysis")
abline(v=Z_score_cutval)

Sens_analysis = as.data.frame(cbind(cutval_sens,hits_sens,hit_rate_sens))
names(Sens_analysis) = c("Cutoff_Z-score","nr_hits","hit_rate")
rm(cutval_sens,hits_sens,hit_rate_sens)

#========= Fix Z-score cutoff and get hits list =================
cutval = rep(Z_score_cutval,length(Means))
hitlist = list()
for(i in 1:length(Promoters))
{
  hits_p = c("Water", row.names(P_vals)[grep(T,P_vals[i] < 0.05)])
  hits_m = c("Water", row.names(Means)[grep(T,abs(Means[i]) > cutval[i])]) 
  hits = intersect(hits_p, hits_m)
  promoter = Promoters[i]
  promoter_scores = Scores[grep(promoter,names(Scores))]
  
  for(h in 1:length(hits))
  {
    hit = hits[h]
    hit_scores = promoter_scores[grep(hit,row.names(promoter_scores)),]
    hit_scores = as.vector(as.matrix(hit_scores))
    if(h==1)
    {hits_scores_promoter = list()}
    hits_scores_promoter[[length(hits_scores_promoter)+1]] = hit_scores
    names(hits_scores_promoter)[length(hits_scores_promoter)] = hit
  }
  
  hitlist[[length(hitlist)+1]] = hits_scores_promoter
}
names(hitlist) = Promoters

#==== plot volcano plots of pval and means per promoter ====

par(mfrow=c(2,2))
for(i in 1:length(Promoters))
{
  plot(Means[[i]],-log10(P_vals[[i]]), main = Promoters[i],ylab="-log10(p-val)",xlab="Median(scores)",col=Prom_clrs[i])
  points(Means[match(names(hitlist[[i]]),row.names(Means)),i],-log10(P_vals[match(names(hitlist[[i]]),row.names(P_vals)),i]),pch=19,col=Prom_clrs[i])
  
  if(length(hitlist[[i]])>1)
  {text(Means[match(names(hitlist[[i]]),row.names(Means)),i],-log10(P_vals[match(names(hitlist[[i]]),row.names(P_vals)),i]),labels = names(hitlist[[i]]))}
  
  abline(v = cutval[i])
  abline(v = -cutval[i])
  abline(h = -log10(0.05))
}

#==== plot volcano plots of pval and means all together ====
par(mfrow=c(1,1))
here_pch = c(1,17,19,17,17,19,19,17)
here_cols = c(Prom_clrs[1], Prom_clrs[5], Prom_clrs[3], Prom_clrs[8], Prom_clrs[7], Prom_clrs[2], Prom_clrs[4], Prom_clrs[6])
for(i in 2:length(Promoters)) #omit EVC
{
  if(i==2)
  {plot(Means[[i]],-log10(P_vals[[i]]),col="darkgrey", cex = 2,
        main = "Volcano plot",ylab="-log10(p-val)",xlab="Mean(Z-scores)",
        ylim=c(0,5),xlim=c(c(-7,7)))} 
  else
  {points(Means[[i]],-log10(P_vals[[i]]),col="darkgrey", cex = 2)}
  
  points(Means[match(names(hitlist[[i]]),row.names(Means)),i],-log10(P_vals[match(names(hitlist[[i]]),row.names(P_vals)),i]),col=here_cols[i],pch=here_pch[i], cex = 2)
  
  #if(length(hitlist[[i]])>1)
  #{text(Means[match(names(hitlist[[i]]),row.names(Means)),i],-log10(P_vals[match(names(hitlist[[i]]),row.names(P_vals)),i]),labels = names(hitlist[[i]]))}
}
abline(v = Z_score_cutval)
abline(v = -Z_score_cutval)
abline(h = -log10(0.05))
legend("topleft",legend = Promoters[-1],col=here_cols[-1],bty="n",pch=here_pch[-1],y.intersp = 0.5, cex = 2)
rm(here_pch)


#====== Get format to plot network ====

Edges = c()
Edge_weight = c()
for(i in 1:length(hitlist))
{
  prom = names(hitlist[i])
  hits = names(hitlist[[i]])
  
  hits_means = unlist(lapply(hitlist[[i]],FUN=mean, na.rm=T))
  
  #edge_weight = Weighted_average_scores[match(hits,row.names(Weighted_average_scores)),
                                        #match(prom,names(Weighted_average_scores))]
  edge_weight = hits_means 
  Edge_weight = c(Edge_weight,edge_weight)
  Edges = c(Edges,c(rbind(rep(prom,length(hits)),hits)))
}

#Create network
net = make_graph(Edges,directed=F)

#Give node attributes: drug or promoter
nodes = vertex.attributes(net)$name
node_type = rep("drug",length(nodes))
node_type[grep(T,nodes %in% Promoters)]="promoter"
#graph_attr(net, name="type") <- node_type
V(net)$type <- node_type

V(net)[V(net)$type=="drug"]$color=Prom_clrs[1]
V(net)[V(net)$type=="promoter"]$color=Prom_clrs[4]

#Give edge attributes: width reflect weighted average Z-score or p-val
E(net)$width = Edge_weight

#====== Get cytoscape file ====

for(i in 1:length(hitlist))
{
  promoter = names(hitlist[i])
  hits = names(hitlist[[i]])
  hits_means = unlist(lapply(hitlist[[i]],FUN=mean, na.rm=T))
  direction = rep("activation",length(hits))
  direction[grep(T,hits_means<0)] = "down-regulation"
  
  promoter_hits = as.data.frame(cbind(rep(promoter,length(hits)),
                      direction,hits,
                      hits_means))
  
  if(i==1)
  {Cytoscape_network = promoter_hits} else
  {Cytoscape_network = as.data.frame(rbind(Cytoscape_network,promoter_hits))}
  
}
names(Cytoscape_network) = c("Promoter","Direction","Drug","Mean_scores")
Cytoscape_network = Cytoscape_network[-grep("Water",Cytoscape_network$Drug),]
row.names(Cytoscape_network) = NULL

#prepare Node_attributes file
Node_attributes = as.data.frame(rbind(
  cbind(unique(Cytoscape_network$Drug),
        rep("Drug",length(unique(Cytoscape_network$Drug))),
        table(Cytoscape_network$Drug)),
  cbind(unique(Cytoscape_network$Promoter),
        rep("Regulator",length(unique(Cytoscape_network$Drug))),
        table(Cytoscape_network$Promoter))
))
row.names(Node_attributes) = NULL
names(Node_attributes) = c("Node","Type","Degree")
Node_attributes$Type[grep("pacrAB",Node_attributes$Node)] = "Regulon"
Node_attributes$Type[grep("micF",Node_attributes$Node)] = "Regulon"
Node_attributes$Type[grep("ompF",Node_attributes$Node)] = "Regulon"
Node_attributes$Type[grep("tolC",Node_attributes$Node)] = "Regulon"

# ============ Get general stats plots ============

hit_rate = length(Cytoscape_network[[1]])/(length(Drugs[-grep("Water",Drugs)])*length(Promoters[-grep("EVC",Promoters)]))
total_drugs_hit = unique(Cytoscape_network$Drug)

hits_per_promoter = rev(sort(table(Cytoscape_network$Promoter)))
hits_per_drug = rev(sort(table(Cytoscape_network$Drug)))
hits_per_direction = rev(sort(table(Cytoscape_network$Direction)))

hits_per_promoter_direction = matrix(ncol=length(unique(Cytoscape_network$Promoter)), nrow=length(unique(Cytoscape_network$Direction)))
colnames(hits_per_promoter_direction) = unique(Cytoscape_network$Promoter)
rownames(hits_per_promoter_direction) = c("down-regulation","activation")

for(i in 1:length(unique(Cytoscape_network$Promoter)))
{
hits_per_promoter_direction [2,i] = sum(Cytoscape_network$Direction[grep(unique(Cytoscape_network$Promoter)[i], Cytoscape_network$Promoter)]=="activation")
hits_per_promoter_direction [1,i] = sum(Cytoscape_network$Direction[grep(unique(Cytoscape_network$Promoter)[i], Cytoscape_network$Promoter)]=="down-regulation")
}
hits_per_promoter_direction = hits_per_promoter_direction[,order(colSums(hits_per_promoter_direction),decreasing=T)]

hits_per_drug_direction = matrix(ncol=length(unique(Cytoscape_network$Drug)), nrow=length(unique(Cytoscape_network$Direction)))
colnames(hits_per_drug_direction) = unique(Cytoscape_network$Drug)
rownames(hits_per_drug_direction) = c("down-regulation","activation")

for(i in 1:length(unique(Cytoscape_network$Drug)))
{
  hits_per_drug_direction [2,i] = sum(Cytoscape_network$Direction[grep(unique(Cytoscape_network$Drug)[i], Cytoscape_network$Drug)]=="activation")
  hits_per_drug_direction [1,i] = sum(Cytoscape_network$Direction[grep(unique(Cytoscape_network$Drug)[i], Cytoscape_network$Drug)]=="down-regulation")
}
hits_per_drug_direction = hits_per_drug_direction[,order(colSums(hits_per_drug_direction),decreasing=T)]


barplot(hits_per_promoter_direction,col=Prom_clrs[c(1,4)],main="Hits per promoter",las=2,ylim=c(0,15))

par(mar=c(10,5,5,5))
barplot(hits_per_drug_direction,las=2,col=Prom_clrs[c(1,4)],main="Hits per drug")

ABcount = Cytoscape_network
ABcount$Category = Map$Category[match(ABcount$Drug, Map$Drug)]
barplot(c(sum(ABcount$Category=="Antibiotic"), sum(ABcount$Category!="Antibiotic")), col = Prom_clrs[4], 
        names.arg = c("Antibiotic", "Non-Antibiotic"), ylab = "No. of interactions")
barplot(c(length(subset(ABcount, ABcount$Category=="Antibiotic" & ABcount$Direction=="down-regulation")$Drug),
          length(subset(ABcount, ABcount$Category!="Antibiotic" & ABcount$Direction=="down-regulation")$Drug)), add = T, col = Prom_clrs[1])
legend("topright", legend = paste0(c("activation", "down-regulation")), col = Prom_clrs[c(4,1)], bty = "n", cex = 2, pch = 15)

if(strain=="WT")
{
Novcount = Cytoscape_network
Novcount$Novelty = rep("Yes", length(Novcount$Drug))
Novcount$Novelty[c(18,19,20,24,32,49,50)] = "No"

barplot(c(sum(Novcount$Novelty=="Yes"), sum(Novcount$Novelty!="Yes")), col = Prom_clrs[4], 
        names.arg = c("Novel", "Known"), ylab = "No. of interactions")
barplot(c(length(subset(Novcount, Novcount$Novelty=="Yes" & Novcount$Direction=="down-regulation")$Drug),
          length(subset(Novcount, Novcount$Novelty!="Yes" & Novcount$Direction=="down-regulation")$Drug)), add = T, col = Prom_clrs[1])
legend("topright", legend = paste0(c("activation", "down-regulation")), col = Prom_clrs[c(4,1)], bty = "n", cex = 2, pch = 15)
}

pie(hits_per_direction,col = Prom_clrs[c(4,1)])
text(0,-1,labels = paste0("n=",length(Cytoscape_network[[1]])))

Cytoscape_network$Mean_scores = as.numeric(as.matrix(Cytoscape_network$Mean_scores))
 
p1 <- ggplot(Cytoscape_network,aes(x=Promoter,y=Mean_scores)) +
  geom_boxplot(color = Prom_clrs[4]) +
  geom_jitter(width=0.1,color = Prom_clrs[4]) + 
  
  #scale_color_manual(values = as.matrix(Prom_clrs[-1])) +
  ylim(c(-7,7)) +
  ggtitle("Mean Z-scores per promoter") + ylab("Mean Z-scores")
print(p1)

Means_cidal = Means[grep("Cidal",Map$Action[match(row.names(Means),Map$Drug)]),]
Means_static = Means[grep("Static",Map$Action[match(row.names(Means),Map$Drug)]),]

ggplot_data = as.data.frame(rbind(my_ggplot_converter(Means_cidal,add_classifier="Cidal"),
                                  my_ggplot_converter(Means_static,add_classifier="Static")))
names(ggplot_data) = c("Means","Promoter","Drug","Action")
ggplot_data$Means = as.numeric(ggplot_data$Means)

p_val = wilcox.test(as.vector(as.matrix(Means_cidal)),as.vector(as.matrix(Means_static)))
p_val = round(p_val$p.value,4)

p2 <- ggplot(ggplot_data, aes(x=Action, y=Means)) + 
  geom_violin(color=Prom_clrs[4]) + geom_boxplot(width=0.1,color=Prom_clrs[4]) +
  geom_text(aes(label = paste0("n=",sum(Action=="Cidal")), x=1, y = max(Means) + 0.5), position = position_dodge(1), vjust = 0) +
  geom_text(aes(label = paste0("n=",sum(Action=="Static")), x=2, y = max(Means) + 0.5), position = position_dodge(1), vjust = 0) +
  ggtitle(paste0("Cidal vs Static, Wilcox p-val=",p_val)) + xlab("Mode of action") + ylab("Mean Z-scores")

print(p2)

#c(paste0("n=",sum(ggplot_data$Action=="Cidal")),paste0("n=",sum(ggplot_data$Action=="Static")))

#get min fitness per drug
Min_growth = c()
for(d in 1:length(Drugs))
{
  drug_wells = Map$Well[grep(Drugs[d],Map$Drug)]
  drug_OD_data = All_OD_AUC[match(drug_wells,row.names(All_OD_AUC)),]
  drug_min_growth = quantile(as.vector(as.matrix(drug_OD_data)),probs=0.10,na.rm = T)
  
  Min_growth = c(Min_growth,drug_min_growth)
}
Min_growth = as.data.frame(Min_growth)
row.names(Min_growth) = Drugs

Min_growth = as.data.frame(cbind(Min_growth,rep("no-hit",length(Min_growth[[1]]))))
names(Min_growth)[2] = "hits"
Min_growth$hits[match(names(hits_per_drug),row.names(Min_growth))] = "hit"
Min_growth$Min_growth = as.numeric(as.matrix(Min_growth$Min_growth))
Min_growth = Min_growth[-grep("Water",row.names(Min_growth)),]

p_val = wilcox.test(as.vector(as.matrix(Min_growth$Min_growth[grep("hit", Min_growth$hits)])),as.vector(as.matrix(Min_growth$Min_growth[grep("no-hit", Min_growth$hits)])))
p_val = round(p_val$p.value,4)

p3 <- ggplot(Min_growth, aes(x=hits, y=Min_growth)) + 
  geom_violin(color=Prom_clrs[4]) + geom_boxplot(width=0.1,color=Prom_clrs[4]) +
  geom_jitter(width=0.1,color=Prom_clrs[4]) +
  geom_text(aes(label = paste0("n=",sum(hits=="hit")), x=1, y = max(Min_growth) + 0.5), position = position_dodge(1), vjust = 0) +
  geom_text(aes(label = paste0("n=",sum(hits=="no-hit")), x=2, y = max(Min_growth) + 0.5), position = position_dodge(1), vjust = 0) +
  ggtitle(paste0("Hit vs No-Hit, Wilcox p-val=",p_val))

print(p3)

#What percentage of drugs are giving hits?
pie(table(Min_growth$hits),col=as.vector(as.matrix(ColorScheme[c(1,4)])))
text(0,-1,length(Min_growth[[1]]))

#Library pie
Drug_categories = unique(Map$Drug)
Drug_categories = Map$Target[match(Drug_categories,Map$Drug)]
pievec = c(table(Drug_categories)[3], table(Drug_categories)[4], table(Drug_categories)[8], table(Drug_categories)[9], table(Drug_categories)[10], 
           table(Drug_categories)[1], table(Drug_categories)[2], table(Drug_categories)[5], table(Drug_categories)[6], table(Drug_categories)[7], table(Drug_categories)[11])
piecol = c(colorRampPalette(c(ColorScheme[2,1], ColorScheme[4,1]))(5), colorRampPalette(c(ColorScheme[2,4], ColorScheme[4,4]))(6))
pie(pievec, col = piecol, clockwise = T, init.angle = 180)

#close file
dev.off()

#================= update Results_list ===========
Results_list2 = list(Means, P_vals,Sens_analysis,
                     Cytoscape_network, hits_per_promoter,hits_per_direction,hits_per_drug,hit_rate)
names(Results_list2) = c("Means","P_vals","Sens_analysis",
                         "Cytoscape_network","hits_per_promoter","hits_per_direction", "hits_per_drug","hit_rate")
Results_list = append(Results_list,Results_list2)
rm(Results_list2)

rm(c,d,i,h,p,
   p1,net,drug,drug_mean,drug_scores,hitlist,hits_m,hits_means,hits_p,means,p_vals,p_val,
   promoter,promoter_hits,promoter_scores,water_scores,direction,cutval,
   edge_weight,Edge_weight,Edges,file_id,hit_rate_sens,hits_sens,hit,hit_rate,hit_scores,
   hit_scores,prom,total_drugs_hit, ABcount, Novcount)

save.image(paste0(Out_dir,"ScreenPipeline_4.RData"))


