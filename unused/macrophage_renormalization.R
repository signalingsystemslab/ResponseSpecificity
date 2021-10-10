# renormalization with internal control genes

# devtools::install_github("Zhanglab-SHT/ISnorm", build_vignettes = T)

library(ISnorm)
library(dbscan)
library(Seurat);library(ggpubr);library(ggplot2)
# browseVignettes('ISnorm')

macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = readRDS("output/macrophage_M1_IFNg_500genes_DBEC.rds")
macro = readRDS("output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
macro = readRDS("output/macrophage_baselines_0hr_500genes_DBEC.rds")

macro2 = readRDS("output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds")
macro = readRDS("output/macrophage_M0all_500genes_Dec2020.rds")
macro = readRDS("output/macrophage_BMDM2_WT_MM_500genes.rds") #_DBEC.rds")
macro = subset(macro, subset = (stimulus!="P3CSK"&stimulus!="IFNb"&stimulus!="Unstim"))

macro = readRDS("./output/macrophage_PMs_500genes_Dec2020.rds")

macro = readRDS("./output/macrophage_BMDM_B6.WT_NOD.AireGW_healthy_500genes_DBEC_filtered.rds")
macro = readRDS("./output/macrophage_BMDM1_NOD.AireGW_healthy_sick_500genes_DBEC_filtered.rds")
macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds")

# macro = merge(macro, macro2)
if(0){
  #caluclate gene distances 
  data = macro[["RNA"]]@counts #start from raw counts, could also start from data
  data = data.frame(data)
  gene_dis<-calculate.dis(data, detection_rate=0.90, ncore=5)
  
  #use DBscan for finding IS genes
  spike_candidate<-dbscan.pick(dis=gene_dis, ngene=(1:floor(nrow(gene_dis)/4)) *5)
  str(spike_candidate)
  spike_candidate
  
  #normalize matrix
  candidate_res<-candidate.norm(mat=data, spike_candidate=spike_candidate[-1]) #spike_candidate[-1] #take off first/first few options
  names(candidate_res)
  candidate_res$sf[1:3,1:3]
  candidate_res$inst[1:3,1:3]
  
  # plot every candidate geneset instability score
  candidate_res$inst[1:3,1:3]
  
  #choose best IS geneset
  # undebug(opt.candidate)
  ISnorm_res<-opt.candidate(mat=data, candidate_res=candidate_res, baseline_threshold = .1)
  names(ISnorm_res)
  
  #get normalized matrix
  data.ISnormalized = ISnorm_res$normalized
  
  #view new size factors for each cell
  sizefactors= ISnorm_res$size_factor
  
  macro[['ISnorm']] <- CreateAssayObject(counts = data.ISnormalized)
  if(0){ #if started from macro[["RNA]]@data, store here in [["ISnorm]] data slot
    macro <- SetAssayData(
      object = macro,
      assay = 'ISnorm',
      new.data = log2(as.matrix(data.ISnormalized)+1),
      slot = 'data'
    )
  }else( #if started from macro[["RNA]]@counts, store here
    macro <- SetAssayData(
      object = macro,
      assay = 'ISnorm',
      new.data = as.matrix(data.ISnormalized),
      slot = 'counts'
    )
  )
}

# csv.data = as.data.frame(macro[["ISnorm"]]@data)
# write.csv(csv.data, "output/macrophage_M0all_500genes_Dec2020_ISnorm.data.csv",row.names = T)

#plot--------------------------------------------------
macro@meta.data$timept = gsub("0hr","0.0hr", macro@meta.data$timept)
macro@meta.data$timept = gsub("24hr","x24hr", macro@meta.data$timept)

p1=VlnPlot(subset(macro, subset= (stimulus=="LPS"|stimulus=="Unstim")&replicate=="rep2"), assay = "ISnorm", slot = "data",
           features = "Adgre1", 
        group.by = "timept",y.max = 12)
p2=VlnPlot(subset(macro, subset= (stimulus=="LPS"|stimulus=="Unstim")&replicate=="rep2"), assay = "RNA", 
           features = "Adgre1", 
        group.by = "timept",y.max = 12)
p1|p2

#plot single timepoint
p1=VlnPlot(subset(macro, subset= (timept=="3hr"|timept=="0.0hr") &replicate=="rep1"), assay = "ISnorm", 
        features = "Ifit3", pt.size = 0.1,
        group.by = "stimulus", y.max = 9)+theme(legend.position = "None")
p2=VlnPlot(subset(macro, subset= (timept=="3hr"|timept=="0.0hr") &replicate=="rep2"&stimulus!="IFNb"), assay = "ISnorm", 
        features = "Ifit3", pt.size = 0.1,
        group.by = "stimulus",y.max = 9)+theme(legend.position = "None")
p1|p2

#plot BMDM2
gene = "Ccl5"; y.max=15
p1=VlnPlot(subset(macro, subset= type =="BMDM2_WT"|type =="BMDM2_MM"), assay = "ISnorm", 
           features = gene, split.by = "type",split.plot = F, pt.size = 0.1,
           group.by = "stimulus", y.max = y.max)
p2=VlnPlot(subset(macro, subset= type =="BMDM2_WT"|type =="BMDM2_MM"), assay = "RNA", 
           features = gene, split.by = "type",split.plot = F, pt.size = 0.1,
           group.by = "stimulus",y.max = y.max)
p1|p2

# saveRDS(macro, "output/macrophage_M0all_500genes_Dec2020.rds")
# saveRDS(macro, "output/macrophage_M0_rep2only_500genes_DBEC.rds")
# saveRDS(macro, "output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds")
# saveRDS(macro, "output/macrophage_M1_IFNg_500genes_DBEC.rds")
# saveRDS(macro, "output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
# saveRDS(macro, "output/macrophage_BMDM2_WT_MM_500genes_DBEC.rds")
# saveRDS(macro, "output/macrophage_BMDM2_WT_MM_500genes.rds")
# saveRDS(macro, "output/macrophage_BMDM2_WT_MM_500genes_LPT.rds")
# saveRDS(macro, "./output/macrophage_PMexpts_Feb2021_500genes_DBEC.rds")
# saveRDS(macro, "./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds")
# saveRDS(macro, "output/macrophage_baselines_0hr_500genes_DBEC.rds")
# saveRDS(macro, "./output/macrophage_PMs_500genes_Dec2020.rds")
# saveRDS(macro, "./output/macrophage_BMDM_B6.WT_NOD.AireGW_healthy_500genes_DBEC_filtered.rds")
# saveRDS(macro, "./output/macrophage_BMDM1_NOD.AireGW_healthy_sick_500genes_DBEC_filtered.rds")

# plot baseline markers----------------------------------------------
macro = readRDS("output/macrophage_baselines_0hr_500genes_DBEC.rds")
VlnPlot(macro, features = c("Adgre1","Rela", "Jun", "Irf3"), group.by = "type", assay = "ISnorm" , y.max = 8, ncol = 4)
VlnPlot(macro, features = c("Adgre1","Nfkbia", "Il1b"), group.by = "type", assay = "ISnorm" , y.max = 8, ncol = 4)

DimPlot(macro, reduction = 'umap', group.by = "type")
DimPlot(macro, reduction = 'umap', group.by = "stimulus")

FeaturePlot(macro, reduction = 'umap', features = c("Arg1", "Egr2", "Retnla","Chil3"))
p1=VlnPlot((macro), assay = "ISnorm", features = c("Arg1", "Egr2", "Retnla","Chil3"),  pt.size = 0.1,
           group.by = "type", y.max = 11, ncol = 4) +
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "green") 
p2=VlnPlot((macro), assay = "ISnorm", features = c("Cd86"),  pt.size = 0.1,
           group.by = "type", y.max = 8) +
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "green") 
p1|p2

macro= readRDS("output/macrophage_M0M1M2_combined_500genes_DBEC.rds")
VlnPlot(subset(macro, timept=="0.5hr"|timept=="0.25hr"), assay = "ISnorm", features = c("Nos2","Cxcl10"),  pt.size = 0.1,
        group.by = "type", y.max = 8,  ncol = 4) +
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "green") 

# plot multiple for ISnorm counts-------------------------------------------------------------
# plotted Tnf, Cxcl10, Il12b
if(1){
macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
gene = "Rhof"
ymax = 10
pt.size = 0
p1=VlnPlot(object = subset(macro, subset= timept=="0.25hr"|timept=="0.0hr"), pt.size = pt.size,
        features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p2=VlnPlot(object = subset(macro, subset= timept=="1hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p3=VlnPlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p4=VlnPlot(object = subset(macro, subset= timept=="8hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p.g1 = p1/p2/p3/p4

macro = readRDS("output/macrophage_M1_IFNg_500genes_DBEC.rds")
p1=VlnPlot(object = subset(macro, subset= timept=="0.5hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p2=VlnPlot(object = subset(macro, subset= timept=="1hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p3=VlnPlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p4=VlnPlot(object = subset(macro, subset= timept=="8hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p.g2 = p1/p2/p3/p4
p.g2


macro = readRDS("output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
p1=VlnPlot(object = subset(macro, subset= timept=="0.5hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p2=VlnPlot(object = subset(macro, subset= timept=="1hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p3=VlnPlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p4=VlnPlot(object = subset(macro, subset= timept=="8hr"|timept=="0.0hr"), pt.size = pt.size,
           features = c(gene), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p.g3 = p1/p2/p3/p4
p.g3
p.g1|p.g2|p.g3
}

#plot multiple, group by timepoint----
if(1){
  macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
  gene = "Tnf"
  ymax = 10
  pt.size = 0
  p1=VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p2=VlnPlot(object = subset(macro, subset= stimulus=="PIC"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p3=VlnPlot(object = subset(macro, subset= stimulus=="IFNb"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p4=VlnPlot(object = subset(macro, subset= stimulus=="P3CSK"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p5=VlnPlot(object = subset(macro, subset= stimulus=="CpG"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p6=VlnPlot(object = subset(macro, subset= stimulus=="TNF"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p.g1 = (p1|p2|p3)/(p4|p5|p6)
  p.g1
  
  macro = readRDS("output/macrophage_M1_IFNg_500genes_DBEC.rds")
  ymax = 10
  p1=VlnPlot(object = subset(macro, subset= (stimulus=="LPS"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p2=VlnPlot(object = subset(macro, subset= (stimulus=="PIC"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p3=VlnPlot(object = subset(macro, subset= (stimulus=="IFNb"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p4=VlnPlot(object = subset(macro, subset= (stimulus=="P3CSK"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p5=VlnPlot(object = subset(macro, subset= (stimulus=="CpG"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p6=VlnPlot(object = subset(macro, subset= (stimulus=="TNF"|timept=="0.0hr")&timept!="x24hr"), pt.size = pt.size,cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  
  p.g2 = (p1|p2|p3)/(p4|p5|p6)
  p.g2
  
  
  macro = readRDS("output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
  ymax = 10
  p1=VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p2=VlnPlot(object = subset(macro, subset= stimulus=="PIC"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p3=VlnPlot(object = subset(macro, subset= stimulus=="IFNb"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p4=VlnPlot(object = subset(macro, subset= stimulus=="P3CSK"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p5=VlnPlot(object = subset(macro, subset= stimulus=="CpG"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  p6=VlnPlot(object = subset(macro, subset= stimulus=="TNF"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 2, colour = "green") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
  
  p.g3 = (p1|p2|p3)/(p4|p5|p6)
  p.g3
  p.g1|p.g2|p.g3
}
