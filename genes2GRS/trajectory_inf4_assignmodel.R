# single cell GRN fitting, trajectory inference, 

# redoing the GRN assignments by testing the main seven GRNs on each gene, with set parameters, 
# and calculating the cost for each one for each gene, to select the best assignment
# probably some caveats to this, is that the selection will depend on the initial parameters

# no need to divide into clusters, can also use more genes, not just ones in Cheng et al clusters

# load libraries
library(deSolve);library(optimx);library(gdata); library(reshape2) ; library(ggalt); library(dplyr)
library(dynamicTreeCut); library(ggplot2);library(amap); library(readxl);library(FME)
library(RColorBrewer);library(reshape2);library(ggdendro);library(grid);library(ape);library(extrafont)
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
# devtools::install_github("KSheu/ksheu.library1")
# install.packages("devtools")
# BiocManager::install("Seurat")
# library(Hmisc)
# example(deSolve)
# vignette("FME")
# vignette("FMEdyna")
# vignette("FMEsteady")
# vignette("FMEother")
# vignette("FMEmcmc")

#read data------------------------------------------

#from normalized and downsampled matrix----
library(ksheu.library1);library(pheatmap);library(RColorBrewer)
# macro = readRDS("./output/macrophage_tutorial3_4sets_subsample_sctransform.rds")
# macro = readRDS("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
# macro = readRDS("./output/macrophage_M0all_gt80_500genes_DBEC.rds")
macro = readRDS("./output/macrophage_M0all_500genes_Dec2020.rds")
# macro = readRDS("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds")
# macro = readRDS("./output/macrophage_tutorial3_subsample_LPS_sctransform.rds") #test with the subset, no unstims
# macro.0hr = readRDS(paste0("./output/macrophage_tutorial3_subsample_0hr_sctransform.rds"))
# macro = merge(macro.0hr, macro)

data = macro[["ISnorm"]]@data #previously everything fit to RNA, but ISnorm expands the upper bound
data = macro[["RNA"]]@data 
data = data.frame(data)
meta = macro@meta.data
meta$stimulus = gsub("IFNb", "IFNB", meta$stimulus)
colnames(data) = paste0(meta$stimulus, "_",meta$timept)
my.dataframe = data.frame(cbind(label = colnames(data), t(data)))

head(my.dataframe)[1:5]
str(my.dataframe)
# my.dataframe$time = gsub(".*_","",my.dataframe$label)
frame = my.dataframe[!grepl("LPSlo", my.dataframe$label),]
frame = cbind(label = frame$label,  as.data.frame(sapply(frame[,-1], as.numeric)) )
frame$stim = gsub("_.*", "", frame$label)
frame$timept = gsub(".*_", "", frame$label)
frame$timept = (gsub("hr", "", frame$timept))
frame$replicate = as.factor(macro@meta.data$replicate[!grepl("LPSlo", macro@meta.data$stimulus)])

# plot a gene ---------------------------------------------------------
gene = "Cxcl10"
genes = colnames(frame)[-c(1, ncol(frame)-2)]

ggplot(frame[!grepl("24", frame$timept),], aes((timept), get(gene)) )+ ylab(gene)+
  geom_violin()+  geom_point(aes(color = as.factor(replicate)), position = "jitter", size = 0.2, alpha = 0.7)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "right")+
  stat_summary(fun = mean, geom='point', size = 2, colour = "black")+
  stat_summary(fun = median, geom='point', size = 2, colour = "green")+
  facet_wrap(~stim, scales = "free_x", nrow = 1)

ggplot(frame[!grepl("24", frame$timept),], aes((timept), get(gene)) )+ ylab(gene)+
  geom_violin()+  geom_point(aes(color = as.factor(replicate)), position = "jitter", size = 0.2, alpha = 0.7)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "right")+
  stat_summary(fun = mean, geom='point', size = 2, colour = "black")+
  stat_summary(fun = median, geom='point', size = 2, colour = "green")+
  facet_wrap(~stim+replicate, scales = "free_x", nrow = 1)

# input trajectories --------------------------------------------------
tf.ap1=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_AP1activationBMDM.txt")
tf.nfkb=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_NFkBactivationBMDM.txt")
tf.irf=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_IRFactivationBMDM.txt")
tf.p38=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/Table_p38activationBMDM.txt")

tf.table.m = rbind(melt(data.frame(tf.ap1, tf = "AP1")), 
                   melt(data.frame(tf.nfkb, tf = "NFkB")), 
                   melt(data.frame(tf.irf, tf = "IRF")),
                   melt(data.frame(tf.p38, tf = "p38")))
tf.table.m$variable = as.numeric(gsub("X","",tf.table.m$variable))
colnames(tf.table.m) = c("stim", "tf","time", "value")
tf.table.m$stim = gsub("CPG", "CpG", tf.table.m$stim)
tf.table.m$stim = gsub("PAM", "P3CSK", tf.table.m$stim)

ggplot(tf.table.m, aes(time, value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 1)+facet_wrap( ~stim, ncol = 7)+theme_bw(base_size = 16)
ggplot(tf.table.m, aes((time), value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 0.5, linetype= "dashed")+
  geom_xspline(aes(color = tf, group = tf), spline_shape=-0.2, size=1)+
  facet_wrap( ~stim, ncol = 7)+theme_bw(base_size = 16) +xlim(0,300)

ggplot(tf.table.m[grepl("LPS|CpG|^PIC|P3CSK|TNF|IFNB", tf.table.m$stim),], aes((time), value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 0.5, linetype= "dashed")+
  geom_xspline(aes(color = tf, group = tf), spline_shape=-0.2, size=1)+
  facet_wrap( ~stim, ncol = 3)+theme_bw(base_size = 16) +xlim(0,300)

# interpolated------------------------------------------------------
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="AP1", ]; plot(z$time, z$value)
input_lps.ap1 <-data.frame(xspline(z$time, z$value, shape=-0.3, draw=F)); ggplot(input_lps.ap1, aes(x,y)) +geom_point()+ylim(0,1)
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="NFkB", ]; plot(z$time, z$value)
input_lps.nfkb <-data.frame(xspline(z$time, z$value, shape=-0.3, draw=F)); ggplot(input_lps.nfkb, aes(x,y)) +geom_point()+ylim(0,1)
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="IRF", ]; plot(z$time, z$value)
input_lps.irf <-data.frame(xspline(z$time, z$value, shape= 0.4, draw=F)); ggplot(input_lps.irf, aes(x,y)) +geom_point()+ylim(0,1)



# load data per gene --------------------------------------------------
my.dataframe.select = (my.dataframe[!grepl("LPSlo", my.dataframe$label),])
tmp = as.data.frame(sapply(my.dataframe.select[,-1], as.numeric)) 
my.dataframe.select = cbind(label = my.dataframe.select$label, (tmp))
str(my.dataframe.select)
my.dataframe.selectmean =  aggregate(my.dataframe.select[,-1], 
                                     by=list(my.dataframe.select$label), FUN=mean) 
my.dataframe.selectmean$stim = gsub("_.*", "", my.dataframe.selectmean$Group.1)

#get the gene ----------------------------------------------------------
# df=my.dataframe.selectmean[!grepl("Unstim", my.dataframe.selectmean$Group.1), c("Group.1",gene)]
df=my.dataframe.selectmean[, c("Group.1",gene)] #keep the 0hr
# df$time = c(rep(c(30, 60, 180, 300, 480), 6),0)
df$time = c(c(15, 30, 60, 180, 300, 480), c(15, 60, 180, 480), rep(c(15, 30, 60, 180, 300, 480), 4), 0)
names(df)=c("label","species", "time")

unstim.cpg = data.frame(label = "CpG_0hr", species = df$species[nrow(df)], time = 0) #repeat Unstim for each stimulus
unstim.ifnb = data.frame(label = "IFNB_0hr", species = df$species[nrow(df)], time = 0)
unstim.lps = data.frame(label = "LPS_0hr", species = df$species[nrow(df)], time = 0)
unstim.p3c = data.frame(label = "P3CSK_0hr", species = df$species[nrow(df)], time = 0)
unstim.pic = data.frame(label = "PIC_0hr", species = df$species[nrow(df)], time = 0)
unstim.tnf = data.frame(label = "TNF_0hr", species = df$species[nrow(df)], time = 0)
df = rbind(df[-nrow(df),],  unstim.cpg, unstim.ifnb, unstim.lps, unstim.p3c, unstim.pic, unstim.tnf)
# plot data mean
df.m=melt(df[,-1],id.vars=c("time"),variable.name="label",value.name="species")
df.m$label = gsub("_.*", "", df$label)
ggplot(data=df.m, aes(x=time,y=species,color=label))+
  geom_line()+xlim(0,500)+
  geom_point(size=3)+facet_wrap(~label, nrow = 1)


normalize <- function(x){
  return(x/ (max(x)))
}
data_scaled <- df.m %>%
  # group_by(label) %>%
  mutate(mRNA = normalize(species))
ggplot(data=data_scaled, aes(x=time,y=mRNA,color=label))+
  geom_line()+ylim(0,1)+xlim(0,500)+
  geom_point(size=3)+facet_wrap(~label, nrow = 1)



# figureout which cluster for the gene--------------------------------
# not needed for this script
# clusters = read_excel("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/1-s2.0-S2405471217300121-mmc4.xlsx")
# genes = colnames(my.dataframe.selectmean)[colnames(my.dataframe.selectmean) %in% clusters$Symbol]
# clusters.subset = clusters[clusters$Symbol %in% genes,]

# differential equations----------------------------------------------
odeModel_steadystate <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    # print(paste0("current St_name ",St_name))
    # St_name = "TNF"
    # gene.clust = clusters$`Cluster ID`[clusters$Symbol==gene]
    
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    # TFA <-approxfun(TimepointsA[-length(TimepointsA)], TFA_profile[-length(TFA_profile)])
    # TFN <-approxfun(TimepointsN[-length(TimepointsN)], TFN_profile[-length(TFN_profile)])
    # TFI<-approxfun(TimepointsI[-length(TimepointsI)], TFI_profile[-length(TFI_profile)])
    
    ktA<-Pars[[1]]
    ktN<-Pars[[2]]
    ktI<-Pars[[3]]
    
    #synthesis and deg rates are set to be the same in the ODEs below
    kdg1s<-Pars[[4]]
    kdg2s<-Pars[[5]]
    kdg2l<-Pars[[6]]
    kdg7s<-Pars[[7]]
    kdg3s<-Pars[[8]]
    kdg3l<-Pars[[9]]
    kdg10l<-Pars[[10]]
    k0<-Pars[[11]]
    
    n = 1
    tau = 0
    
    #TF activity forcing function 
    ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    
    # tf_general = mean(ap1, nfkb, irf)
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0
    
    #AND gate
    # fAIN<-(1.0-k0)*((ktA*ktN*ktI*ap1*nfkb*irf)^n /(1.0+ktA*ap1+ktN*nfkb+ktI*irf+ktA*ktN*ap1*nfkb+
    #                                                  ktN*ktI*nfkb*irf+ktA*ktI*ap1*irf+
    #                                                  ktA*ktN*ktI*ap1*nfkb*irf)^n)+k0
    # mRNA ODEs
    if(1){
      dG1S <-kdg1s*fA-kdg1s*G1S #mA
      dG2S <-kdg2s*fN-kdg2s*G2S #mB
      dG2L <-kdg2l*fN-kdg2l*G2L #mC
      
      # dG10L <-kdg10l*fAIN-kdg10l*G10L #mD
      dG10L <-kdg2s*fN-kdg10l*G10L #mD, mod.2S
      
      dG7S <-kdg7s*fIN-kdg7s*G7S #mE
      dG3S <-kdg3s*fI-kdg3s*G3S #mF
      dG3L <-kdg3l*fI-kdg3l*G3L #mG
    }
    
    if(0){ #ksyn=1 for modelv3
      dG1S <-fA-kdg1s*G1S #mA
      dG2S <-fN-kdg2s*G2S #mB
      dG2L <-fN-kdg2l*G2L #mC
      
      # dG10L <-kdg10l*fAIN-kdg10l*G10L #mD
      dG10L <-fN-kdg10l*G10L #mD, mod.2S
      
      dG7S <-fIN-kdg7s*G7S #mE
      dG3S <-fI-kdg3s*G3S #mF
      dG3L <-fI-kdg3l*G3L #mG
    }
    
    
    # return the rate of change
    
    
    
    return( list(c(dG1S,dG2S,dG2L,dG10L,dG7S,dG3S,dG3L)) )
  })
}

odeModel <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    # print(paste0("current St_name ",St_name))
    # St_name = "TNF"
    # gene.clust = clusters$`Cluster ID`[clusters$Symbol==gene]
    
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    # TFA <-approxfun(TimepointsA[-length(TimepointsA)], TFA_profile[-length(TFA_profile)])
    # TFN <-approxfun(TimepointsN[-length(TimepointsN)], TFN_profile[-length(TFN_profile)])
    # TFI<-approxfun(TimepointsI[-length(TimepointsI)], TFI_profile[-length(TFI_profile)])
    
    ktA<-Pars[[1]]
    ktN<-Pars[[2]]
    ktI<-Pars[[3]]
    
    #synthesis and deg rates are set to be the same in the ODEs below
    kdg1s<-Pars[[4]]
    kdg2s<-Pars[[5]]
    kdg2l<-Pars[[6]]
    kdg7s<-Pars[[7]]
    kdg3s<-Pars[[8]]
    kdg3l<-Pars[[9]]
    kdg10l<-Pars[[10]]
    k0<-Pars[[11]]
    
    n = 1
    tau = 0
    
    
    #TF activity forcing function 
    ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    
    if (St_name=="LPS"|St_name=="P3CSK"|St_name=="CpG"){
      mRNA_life_mod = 480*p38
      # print(mRNA_life_mod)
      kdg10l<-log(2)/mRNA_life_mod
    }
    
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0
    
    #AND gate
    fAIN<-(1.0-k0)*((ktA*ktN*ktI*ap1*nfkb*irf)^n /(1.0+ktA*ap1+ktN*nfkb+ktI*irf+ktA*ktN*ap1*nfkb+
                                                     ktN*ktI*nfkb*irf+ktA*ktI*ap1*irf+
                                                     ktA*ktN*ktI*ap1*nfkb*irf)^n)+k0
    # mRNA ODEs
    if(1){
    dG1S <-kdg1s*fA-kdg1s*G1S #mA
    dG2S <-kdg2s*fN-kdg2s*G2S #mB
    dG2L <-kdg2l*fN-kdg2l*G2L #mC
    
    # dG10L <-kdg10l*fAIN-kdg10l*G10L #mD
    dG10L <-kdg2s*fN-kdg10l*G10L #mD, mod.2S
    
    dG7S <-kdg7s*fIN-kdg7s*G7S #mE
    dG3S <-kdg3s*fI-kdg3s*G3S #mF
    dG3L <-kdg3l*fI-kdg3l*G3L #mG
    }
    
    
    if(0){ #ksyn=1 for modelv4, using a p38 input curve to adjust kdg10l
      dG1S <-fA-kdg1s*G1S #mA
      dG2S <-fN-kdg2s*G2S #mB
      dG2L <-fN-kdg2l*G2L #mC
      
      # dG10L <-kdg10l*fAIN-kdg10l*G10L #mD
      dG10L <-fN-kdg10l*G10L #mD, mod.2S
      
      dG7S <-fIN-kdg7s*G7S #mE
      dG3S <-fI-kdg3s*G3S #mF
      dG3L <-fI-kdg3l*G3L #mG
    }
    
  
    # return the rate of change
    
    
    
    return( list(c(dG1S,dG2S,dG2L,dG10L,dG7S,dG3S,dG3L)) )
  })
}

# params and ss.init -------------------------------------------------
# pars = list(n= 1, k0=.1, kd = 0.5, kt=1, k_syn = 1, k_deg=log(2)/30, ap1_input= 0.02, nfkb_input=0.02, irf_input = 0.02)
# pars = c(mRNA = 0.1, ktA =1.48 , ktN = 2.0, ktI = 2.0, k0=.005, k_syn = 1, k_deg=log(2)/30, stimulus = 1)
# pars = c(mRNA = 0.05, n = 1, ktT =0.1  k0=.0005, k_syn = 1, k_deg=log(2)/30, stimulus = 1, tau = 10)
if (1){
mRNA_life<-c(30,60,240,300, 480)
kk<-4
kk_mod = 5
# TF DNA binding rates
kt1<-0.48      #AP1
kt2<-1.3      #NFkB
kt3<-1.25      #IRF

# mRNA degradation rates
kd1s<-log(2)/mRNA_life[1]     #GRN 1S 
kd2s<-log(2)/mRNA_life[1]     #GRN 2S
kd2l<-log(2)/mRNA_life[kk]     #GRN 2L

# kd10l<-log(2)/mRNA_life[kk]     #GRN 10L
kd10l<-log(2)/mRNA_life[1]      #GRN mod.2S

kd7s<-log(2)/mRNA_life[1]#[kk]     #GRN 7S
kd3s<-log(2)/mRNA_life[1]     #GRN 3S
kd3l<-log(2)/mRNA_life[kk]     #GRN 3L

# k0<-0.15;
k0<-0.005;
mRNA = 0.05

pars = c(kt1,kt2,kt3,kd1s,kd2s,kd2l,kd7s,kd3s,kd3l,kd10l,k0, mRNA, stimulus = 1)
stims = c("CpG", "IFNB","LPS","P3CSK", "PIC", "TNF")

bTF1<-0.05     #basal TF activity AP1
bTF2<-0.05     # --- NFkB
bTF3<-0.05     # --- IRF
f1o<-(1.0-k0)*((kt1*bTF1)/(1.0+kt1*bTF1))+k0
f2o<-(1.0-k0)*((kt2*bTF2)/(1.0+kt2*bTF2))+k0
f3o<-(1.0-k0)*((kt3*bTF3)/(1.0+kt3*bTF3))+k0
f7o<-(1.0-k0)*((kt2*bTF2+kt3*bTF3+kt2*kt3*bTF2*bTF3)/(1.0+kt2*bTF2+kt3*bTF3+kt2*kt3*bTF2*bTF3))+k0
f10o<-(1.0-k0)*((kt1*kt2*kt3*bTF1*bTF2*bTF3)/(1.0+kt1*bTF1+kt2*bTF2+kt3*bTF3+kt1*kt2*bTF1*bTF2+
                                                kt2*kt3*bTF2*bTF3+kt1*kt3*bTF1*bTF3+
                                                kt1*kt2*kt3*bTF1*bTF2*bTF3))+k0


# steady.state<- c(G1S = f1o,G2S = f2o,G2L = f2o,G10L = f10o, G7S = f7o,G3S = f3o,G3L = f3o)
steady.state<- c(G1S = f1o,G2S = f2o,G2L = f2o,G10L = f2o, G7S = f7o,G3S = f3o,G3L = f3o)

# steady.state <- runsteady(y = c(steady.state), #init
#                      fun = odeModel_steadystate,
#                      parms = pars, time = c(0,1e5))$y


# times = df$time[1:5]
times = c(0, 15, 30, 60, 180, 300, 480)
# times = seq(0, 480, length=500)
}
# wrapper_steady state -----------------------------------------------------------------
solve_model_steadystate <- function(pars, times = seq(-3000, 0, length=500)) {
  # times = seq(0, 480, length=500)
  cnt = 1
  for (i in stims){
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # print(pars[c("stimulus")])
    
    pars[12] <- list(steady.state)
    
    ss.init <- unlist(pars[12])
    params <- pars[-12]
    out <- ode(ss.init, times, odeModel_steadystate, params)
    out.frame = as.data.frame(out)
    
    out.frame$label = i
    
    if(cnt==1){
      collect_out <- out.frame
      
    }else{
      collect_out <-rbind(collect_out, out.frame)
    }
    cnt=cnt+1
    
  }
  
  return(collect_out)
}

solve.ss = solve_model_steadystate(pars)
steady.state = (solve.ss[nrow(solve.ss),2:8])


# wrapper solve model v3-----------------------------------------------------------------
solve_model <- function(pars, times = seq(0, 480, length=500)) {
  # times = seq(0, 480, length=500)
  cnt = 1
  for (i in stims){
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # print(pars[c("stimulus")])
    
    pars[12] <- list(steady.state)
    
    ss.init <- unlist(pars[12])
    params <- pars[-12]
    out <- ode(ss.init, times, odeModel, params)
    out.frame = as.data.frame(out)
    
    out.frame$label = i
    
    if(cnt==1){
      collect_out <- out.frame
      
    }else{
      collect_out <-rbind(collect_out, out.frame)
    }
    cnt=cnt+1
    
  }
  
  return(collect_out)
}

# cost function ------------------------------------------------------
cost = function(pars, times = c(0, 15, 30, 60, 180, 300, 480), grs_name = "mRNA") {
  
  cnt = 1
  for (i in stims){
    # print(i)

    #the data points
    data = data.frame(data_scaled[data_scaled$label==i, c(1,4)]) #fit to normalized
    # data = df.m[df.m$label==i, c(1,3)] #fit to unnormalized
    colnames(data) = c("time", grs_name)
    
    # pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # # print(pars[c("stimulus")])
    # 
    # pars[12] <- list(steady.state)
    # 
    # ss.init <- unlist(pars[12])
    # params <- pars[-12]
    # out <- ode(ss.init, times, odeModel, params)
    # out.frame = as.data.frame(out)
   
    out.frame = solve2[solve2$time >-100, ]
    out.frame = out.frame[out.frame$label == i,]
    
    if(cnt==1){
      Cost = modCost(out.frame, data, weight = "none") # try weight = "std" or "mean"
      
      # cost_resids = cost$residuals
      # cost_resids$name = i
      # collect_resids = cost_resids
      
    }else{
      Cost = modCost(out.frame, data, cost = Cost) # try weight = "std" or "mean"
      
      # cost_resids = cost$residuals
      # cost_resids$name = i
      # collect_resids <-rbind(collect_resids, cost_resids)
    }
    cnt=cnt+1
    
  }
  
  return ((Cost))
}

# find solutions and costs-----------------------------------

solve = solve_model(pars)
solve = rbind(solve.ss, solve)
solve.truncate = solve[solve$time >-100, ]
G1S.scale = max(solve.truncate$G1S)
G2S.scale = max(solve.truncate$G2S)
G2L.scale = max(solve.truncate$G2L)
G10L.scale = max(solve.truncate$G10L)
G7S.scale = max(solve.truncate$G7S)
G3S.scale = max(solve.truncate$G3S)
G3L.scale = max(solve.truncate$G3L)

#max normalized
solve2 = solve
solve2$G1S = solve2$G1S/G1S.scale
solve2$G2S = solve2$G2S/G2S.scale
solve2$G2L = solve2$G2L/G2L.scale
solve2$G10L = solve2$G10L/G10L.scale
solve2$G7S = solve2$G7S/G7S.scale
solve2$G3S = solve2$G3S/G3S.scale
solve2$G3L = solve2$G3L/G3L.scale
# write.table(solve2[solve2$time >-100, ], "./fit_mean/max_normalized_default_params.txt",quote=F,row.names = F,sep="\t")

solve2.m = melt(solve2[,-1])
solve2.m$time = solve2$time
ggplot(solve2.m, aes(time, value))+geom_line(aes(color = label))+
  facet_grid(rows = vars(variable), cols = vars(label), scales = "free_y" )+
  xlim(-100, 500)

#unnormalized
solve.m = melt(solve.truncate[,-1])
solve.m$time = solve.truncate$time
ggplot(solve.m, aes(time, value))+geom_line(aes(color = label))+
  facet_grid(rows = vars(variable), cols = vars(label), scales = "free_y" )+
  xlim(-100, 500)



# get cost for each GRS_-----------------------------------------------------
collect_cost = data.frame()
for (i in c("G1S", "G2S", "G2L","G10L", "G7S", "G3S", "G3L")){
  cost_gene = cost(pars, grs_name = i)$model
  collect_cost = c(collect_cost, cost_gene)
}

# do for each gene to select best GRS---------------------------------------------------
collect_cost_all = data.frame()
for (i in genes[1:500]){
  skip_to_next <- FALSE
  tryCatch( 
    {
      print(i)
      gene = i
      
      #get the gene --------
      # df=my.dataframe.selectmean[!grepl("Unstim", my.dataframe.selectmean$Group.1), c("Group.1",gene)]
      # df$time = rep(c(30, 60, 180, 300, 480), 6)
      # names(df)=c("label","species", "time")
      
      df=my.dataframe.selectmean[, c("Group.1",gene)] #keep the 0hr
      # df$time = c(rep(c(30, 60, 180, 300, 480), 6),0)
      df$time = c(c(15, 30, 60, 180, 300, 480), c(15, 60, 180, 480), rep(c(15, 30, 60, 180, 300, 480), 4), 0)
      names(df)=c("label","species", "time")
      
      unstim.cpg = data.frame(label = "CpG_0hr", species = df$species[nrow(df)], time = 0) #repeat Unstim for each stimulus
      unstim.ifnb = data.frame(label = "IFNB_0hr", species = df$species[nrow(df)], time = 0)
      unstim.lps = data.frame(label = "LPS_0hr", species = df$species[nrow(df)], time = 0)
      unstim.p3c = data.frame(label = "P3CSK_0hr", species = df$species[nrow(df)], time = 0)
      unstim.pic = data.frame(label = "PIC_0hr", species = df$species[nrow(df)], time = 0)
      unstim.tnf = data.frame(label = "TNF_0hr", species = df$species[nrow(df)], time = 0)
      df = rbind(df[-nrow(df),],  unstim.cpg, unstim.ifnb, unstim.lps, unstim.p3c, unstim.pic, unstim.tnf)
      
      
      # plot data mean
      # df.m=melt(df[,-1],id.vars=c("time"),variable.name="label",value.name="mRNA")
      # df.m$label = gsub("_.*", "", df$label)
      df.m=melt(df[,-1],id.vars=c("time"),variable.name="label",value.name="species")
      df.m$label = gsub("_.*", "", df$label)
      
      
      normalize <- function(x){
        return(x/ (max(x)))
      }
      data_scaled <- df.m %>%
        # group_by(label) %>%
        mutate(mRNA = normalize(species))
      
      collect_cost = data.frame()
      for (s in c("G1S", "G2S", "G2L","G10L", "G7S", "G3S", "G3L")){
        cost_gene = cost(pars, grs_name = s)$model
        collect_cost = data.frame(c(collect_cost, cost_gene) )
      }
      colnames(collect_cost) = c("G1S", "G2S", "G2L","G10L", "G7S", "G3S", "G3L")
      
      collect_cost_all = rbind(collect_cost_all, collect_cost)
      
    }, 
    
    error = function(e) { 
      message(paste0("Error for ", gene))
      skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }    
  
}

#for RNA@data-----------
# collect_cost_all$gene = genes[!grepl("^Il22$|^Il2$|^Klrg1$", genes)]
# collect_cost_all$best = colnames(collect_cost_all)[max.col(collect_cost_all[,-c(8,9)] *-1,ties.method="first")]
# rownames(collect_cost_all) = collect_cost_all$gene
# pie(table(collect_cost_all$best))
# pheatmap(collect_cost_all[, -c(8:9)], scale = "row", cluster_cols = F, clustering_method = "ward.D2",
#          colorRampPalette((brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103) )
# write.table(collect_cost_all, "./fit_mean/reassign_genes_to_grns_unweightedcost_modv3.txt",sep = "\t",quote=F,row.names = F)
 
# collect_cost_all$gene = genes[!grepl("^Ccr6$|^Cd1d2$|^Cd27$|^Cd3d$|^Cd79a$|^Chil3$|^Cpa3$|^Cxcl13$|^Dpp4$|^Eomes$|^Foxp3$|^Gzmb$|^Ifna1$|^Ighd$|^Ighg2c$|^Iglc3$|^Il22$|^Il2$|^Il4$|^Klrg1$|^Pdcd1$|^Tbx21$|^Tigit.BD.custom|^Tnfrsf17$|^Trat1$", genes)]
# collect_cost_all$best = colnames(collect_cost_all)[max.col(collect_cost_all[,-c(8,9)] *-1,ties.method="first")]
# rownames(collect_cost_all) = collect_cost_all$gene
# pie(table(collect_cost_all$best))
# pheatmap(collect_cost_all[, -c(8:9)], scale = "column", cluster_cols = F)
# pheatmap(collect_cost_all[, -c(8:9)], scale = "row", cluster_cols = F)
# # write.table(collect_cost_all, "./fit_mean/reassign_genes_to_grns_meanweightedcost.txt",sep = "\t",quote=F,row.names = F)

# *****For old gene list ******--------------------------------------------------------
# collect_cost_all$gene = genes[!grepl("^Il22$|^Il2$|^Klrg1$", genes)]
# collect_cost_all$best = colnames(collect_cost_all)[max.col(collect_cost_all[,-c(8,9)] *-1,ties.method="first")]
# rownames(collect_cost_all) = collect_cost_all$gene
# pie(table(collect_cost_all$best))
# pheatmap(collect_cost_all[, -c(8:9)], scale = "row", cluster_cols = F, clustering_method = "ward.D2",
#          colorRampPalette((brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103) )
# # write.table(collect_cost_all, "./fit_mean/reassign_genes_to_grns_unweightedcost_modv3_p38input.txt",sep = "\t",quote=F,row.names = F)
# collect_cost_all = read.delim("./fit_mean/reassign_genes_to_grns_unweightedcost_modv3_p38input.txt")

# For new 500 gene list Dec 2020--------------------------------------------------------
collect_cost_all$gene = genes[1:500]
collect_cost_all$best = colnames(collect_cost_all)[max.col(collect_cost_all[,-c(8)] *-1,ties.method="first")]
rownames(collect_cost_all) = collect_cost_all$gene
pie(table(collect_cost_all$best))
pheatmap(collect_cost_all[, -c(8:9)], scale = "row", cluster_cols = F, clustering_method = "ward.D2",
         colorRampPalette((brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103) )
# write.table(collect_cost_all, "./fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_modv3_p38input_500genes_Dec2020.txt",sep = "\t",quote=F,row.names = F)
# write.table(collect_cost_all, "./fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_modv3_p38input_ksynfix_500genes_Dec2020.txt",sep = "\t",quote=F,row.names = F)
# write.table(collect_cost_all, "./fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_modv3_p38input_ksynfix_n3_500genes_Dec2020.txt",sep = "\t",quote=F,row.names = F)
#for ISnorm@data--------
# write.table(collect_cost_all, "./fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_modv3_p38input_ksynfix_500genes_Dec2020_ISnorm.txt",sep = "\t",quote=F,row.names = F)

collect_cost_all = read.delim("./fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_modv3_p38input_ksynfix_500genes_Dec2020_ISnorm.txt", sep = ",")

############################################################################################
# new OdeModel equation, (overwritting the one above)---------------------------------------
# differential equations----------------------------------------------
odeModel_steadystate <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    # print(paste0("current St_name ",St_name))
    # St_name = "TNF"
    gene.clust = collect_cost_all$best[collect_cost_all$gene==gene] #quick check using existing
    # gene.clust = clusters$`Cluster ID`[clusters$Symbol==gene]
    
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    # TFA <-approxfun(TimepointsA[-length(TimepointsA)], TFA_profile[-length(TFA_profile)])
    # TFN <-approxfun(TimepointsN[-length(TimepointsN)], TFN_profile[-length(TFN_profile)])
    # TFI<-approxfun(TimepointsI[-length(TimepointsI)], TFI_profile[-length(TFI_profile)])
    
    n = Pars[c("n")]
    ktA<-Pars[c("ktA")]
    ktN<-Pars[c("ktN")]
    ktI<-Pars[c("ktI")]
    k0 <-Pars[c("k0")]
    k_deg = Pars[c("k_deg")]
    k_syn = k_deg #Pars[c("k_syn")]
    # tau = Pars[c("tau")]
    tau = 0
    
    #TF activity forcing function 
    ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    
    # tf_general = mean(ap1, nfkb, irf)
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0

    
    # mRNA ODEs
    
    if (gene.clust =="G1S"){ #keeping all original names, but using modv3,p38input
      dmRNA <-k_syn*fA - k_deg*mRNA#mA G1S
    }
    if (gene.clust =="G2S"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mB G2S
    }
    if (gene.clust =="G2L"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mC G2L
    }
    if (gene.clust =="G10L"){
      dmRNA <-k_syn*fN - k_deg*mRNA  #mD G10L
    }
    if (gene.clust =="G7S"){
      dmRNA <-k_syn*fIN - k_deg*mRNA #mE G7S
    }
    if (gene.clust =="G3S"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mF G3S
    }
    if (gene.clust =="G3L"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mG G3L
    }
    
    return(list(c(dmRNA) ))
  })
}

odeModel <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    # print(paste0("current St_name ",St_name))
    # St_name = "TNF"
    gene.clust = collect_cost_all$best[collect_cost_all$gene==gene]
    # gene.clust = clusters$`Cluster ID`[clusters$Symbol==gene]
    
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    # TFA <-approxfun(TimepointsA[-length(TimepointsA)], TFA_profile[-length(TFA_profile)])
    # TFN <-approxfun(TimepointsN[-length(TimepointsN)], TFN_profile[-length(TFN_profile)])
    # TFI<-approxfun(TimepointsI[-length(TimepointsI)], TFI_profile[-length(TFI_profile)])

    n = Pars[c("n")]
    ktA<-Pars[c("ktA")]
    ktN<-Pars[c("ktN")]
    ktI<-Pars[c("ktI")]
    k0 <-Pars[c("k0")]
    k_deg = Pars[c("k_deg")]
    k_deg_p38 = Pars[c("k_deg")]
    k_p38 = 1 #Pars[c("k_p38")]
    k_syn = k_deg #Pars[c("k_syn")]
    tau = Pars[c("tau")]
    

    #TF activity forcing function 
    ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    
    if (St_name=="LPS"|St_name=="P3CSK"|St_name=="CpG"){
      mRNA_life_mod = 480*p38*k_p38
      # print(mRNA_life_mod)
      k_deg_p38<-log(2)/mRNA_life_mod
    }
    
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0

    
    # mRNA ODEs
    
    if (gene.clust =="G1S"){ #keeping all original names, but using modv3,p38input
      dmRNA <-k_syn*fA - k_deg*mRNA#mA G1S
    }
    if (gene.clust =="G2S"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mB G2S
    }
    if (gene.clust =="G2L"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mC G2L
    }
    if (gene.clust =="G10L"){
      dmRNA <-k_syn*fN - k_deg_p38*mRNA  #mD G10L, mod. 2S
    }
    if (gene.clust =="G7S"){
      dmRNA <-k_syn*fIN - k_deg*mRNA #mE G7S
    }
    if (gene.clust =="G3S"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mF G3S
    }
    if (gene.clust =="G3L"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mG G3L
    }
    
    return(list(c(dmRNA) ))
  })
}

# params and init-------------------------------------------------------------
pars = c(mRNA = 0.05, n = 1, ktA =0.48 , ktN = 1.3, ktI = 1.25, k0=.005, k_deg=log(2)/30, stimulus = 1, tau = 10)

stims = c("CpG", "IFNB","LPS","P3CSK", "PIC", "TNF")

solve_model_steadystate <- function(pars, times = seq(-3000, 0, length=500)) {
  # times = seq(0, 480, length=500)
  cnt = 1
  for (i in stims){
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # print(pars[c("stimulus")])
    
    ss.init <- pars[c("mRNA")]
    params <- pars[c("n", "ktA", "ktN", "ktI", "k0", "k_deg", "stimulus", "tau")]
    out <- ode(ss.init, times, odeModel_steadystate, params)
    out.frame = as.data.frame(out)
    out.frame$label = i
    
    if(cnt==1){
      collect_out <- out.frame
      
    }else{
      collect_out <-rbind(collect_out, out.frame)
    }
    cnt=cnt+1
    
  }
  
  return(collect_out)
}

solve_model <- function(pars, times = seq(-500, 480, length=500)) {
  # times = seq(0, 480, length=500)
  
  solve.ss = solve_model_steadystate(pars)
  steady.state = (solve.ss[nrow(solve.ss),2])
  
  cnt = 1
  for (i in stims){
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # print(pars[c("stimulus")])
    
    pars[c("mRNA")] <- runsteady(y = c(mRNA=0), #init
                                 fun = odeModel_steadystate,
                                 parms = pars, time = c(0,1e5))$y
    
    # pars[c("mRNA")] <- steady.state
    
    ss.init <- pars[c("mRNA")]
    params <- pars[c("n", "ktA", "ktN", "ktI", "k0", "k_deg", "stimulus", "tau")]
    out <- ode(ss.init, times, odeModel, params)
    out.frame = as.data.frame(out)
    out.frame$label = i
    
    if(cnt==1){
      collect_out <- out.frame
      
    }else{
      collect_out <-rbind(collect_out, out.frame)
    }
    cnt=cnt+1
    
  }
  
  return(collect_out)
}



#new cost function - overwirtting the above-------------------------------------
cost = function(pars, times = c(0, 15, 30, 60, 180, 300, 480),  grs_name = "mRNA") {
  
  solve = solve_model(pars)
  # solve = rbind(solve.ss, solve)
  solve.truncate = solve[solve$time >-100, ]
  mRNA.scale = max(solve.truncate$mRNA)
  solve2 = solve
  solve2$mRNA = solve2$mRNA/mRNA.scale
  
  cnt = 1
  for (i in stims){
    # print(i)
    
    #the data points
    data = data.frame(data_scaled[data_scaled$label==i, c(1,4)]) #fit to normalized
    # data = df.m[df.m$label==i, c(1,3)] #fit to unnormalized
    colnames(data) = c("time", grs_name)
  
    
    out.frame = solve2[solve2$time >-100, ]
    out.frame = out.frame[out.frame$label == i,]
  
    if(cnt==1){
      Cost = modCost(out.frame, data, weight = "none") # try weight = "std" or "mean"
      
      # cost_resids = cost$residuals
      # cost_resids$name = i
      # collect_resids = cost_resids
      
    }else{
      Cost = modCost(out.frame, data, cost = Cost) # try weight = "std" or "mean"
      
      # cost_resids = cost$residuals
      # cost_resids$name = i
      # collect_resids <-rbind(collect_resids, cost_resids)
    }
    cnt=cnt+1
    
  }
  
  return ((Cost))
  # return ((Cost$model)) # for GenSA
}




# fitting-------------------------------------------------------------
fit = modFit(f = cost, p = pars, lower = 0.001,  method = "Marq") #
fit
summary(fit)

# optim simulated annealing ------------------------------------------
# fit.sann = modFit(f = cost, p = pars, lower = 0 , method = "SANN") #"Nelder-Mead"

# try simulation annealing with GenSA-----------------------------------
cost_SAwrapper = function(pars.fit, times = c(0, 15, 30, 60, 180, 300, 480),  grs_name = "mRNA") {
  pars2 = pars.fit
  cost.sa = cost(pars = pars2)
  
  # update BEST_SCORE
  if(cost.sa$model < BEST_SCORE[length(BEST_SCORE)]) {
    BEST_SCORE <<- c(BEST_SCORE, cost.sa$model)}else{
      BEST_SCORE <<- c(BEST_SCORE,BEST_SCORE[length(BEST_SCORE)])}
  
  # plt.bestScore
  NUM_ITER <<- NUM_ITER + 1
  cat('iter:',NUM_ITER,'pars:',unlist(pars2),'\n')
  cat('score:',signif(cost.sa$model,2),'Best_score',signif(BEST_SCORE[length(BEST_SCORE)],2))
  
  return ((cost.sa$model)) # for GenSA
}



library(GenSA)
set.seed(1234)
NUM_ITER <- 0; BEST_SCORE<- 9999
fit.sann <- GenSA(par = pars, fn = cost_SAwrapper,
                    lower=rep(0.0001,9), upper = c(0.05, 6, 5, 5, 5, 1, log(2)/15, 1, 100),
                    control=list(verbose=T,
                                 nb.stop.improvement=50,
                                 max.call = 1000))


# rerun the model and plot---------------------------------------------
pars.fit <- fit$par
optim <- solve_model(pars.fit)
cost.fitted = cost(pars.fit)$model
# optim <- solve_model(pars)
# cost.fitted = cost(pars)$model
ggplot(NULL, aes(time, mRNA, color = label))+
  geom_line(data = optim)+
  geom_point(data = data_scaled)+
  xlim(-100, 500)+ ylim(0,1.1)+
  facet_wrap(~label, nrow = 1)

# for unnormalized
ggplot(NULL, aes(time, mRNA, color = label))+
  geom_line(data = optim)+
  geom_point(data = df.m)+
  xlim(-100, 500)+ # ylim(0,1.2)+
  facet_wrap(~label, nrow = 1)


collect_out_norm<- optim %>%
  mutate(mRNA = normalize(mRNA))

ggplot(NULL, aes(time, mRNA, color = label))+
  geom_line(data = collect_out_norm)+
  geom_point(data = data_scaled) +
  xlim(-100, 500)+# ylim(0,1.2)+
  facet_wrap(~label, nrow = 1)


# fit for every gene and save the fit.par------------------------------
# induced_genes = read.delim("F:/scRNAseq_macro/bulk_rnaseq/Cheng2017_induced.txt")
# induced_genes = induced_genes$gene
# genes_fit = genes[genes %in% induced_genes]
genes_fit = c("Tnf", "Ccl5","Il6","Cxcl10", "Ifit3", "Cmpk2")
# genes_fit_cont = genes_fit[74:105]

count = 1
for (i in genes_fit){
  
  skip_to_next <- FALSE
  
  
  tryCatch( 
    {
      print(i)
      gene = i
      
      #get the gene ----------------------------------------------------------
      # df=my.dataframe.selectmean[!grepl("Unstim", my.dataframe.selectmean$Group.1), c("Group.1",gene)]
      df=my.dataframe.selectmean[, c("Group.1",gene)] #keep the 0hr
      # df$time = c(rep(c(30, 60, 180, 300, 480), 6),0)
      df$time = c(c(15, 30, 60, 180, 300, 480), c(15, 60, 180, 480), rep(c(15, 30, 60, 180, 300, 480), 4), 0)
      names(df)=c("label","species", "time")
      
      unstim.cpg = data.frame(label = "CpG_0hr", species = df$species[nrow(df)], time = 0) #repeat Unstim for each stimulus
      unstim.ifnb = data.frame(label = "IFNB_0hr", species = df$species[nrow(df)], time = 0)
      unstim.lps = data.frame(label = "LPS_0hr", species = df$species[nrow(df)], time = 0)
      unstim.p3c = data.frame(label = "P3CSK_0hr", species = df$species[nrow(df)], time = 0)
      unstim.pic = data.frame(label = "PIC_0hr", species = df$species[nrow(df)], time = 0)
      unstim.tnf = data.frame(label = "TNF_0hr", species = df$species[nrow(df)], time = 0)
      df = rbind(df[-nrow(df),],  unstim.cpg, unstim.ifnb, unstim.lps, unstim.p3c, unstim.pic, unstim.tnf)
      # plot data mean
      df.m=melt(df[,-1],id.vars=c("time"),variable.name="label",value.name="species")
      df.m$label = gsub("_.*", "", df$label)
      
      
      normalize <- function(x){
        return(x/ (max(x)))
      }
      data_scaled <- df.m %>%
        # group_by(label) %>%
        mutate(mRNA = normalize(species))
      
      fit = modFit(f = cost, p = pars, lower =0.001,  method = "Marq") #
      saveRDS(fit, paste0("./fit_mean_modv3_p38input/", gene, ".rds"))
      fitpar = data.frame(t(fit$par))
      fitpar$gene = i
      
      if(count == 1){
        collect.fitpars = fitpar
      }else{
        collect.fitpars = rbind(collect.fitpars, fitpar)
      }
      
      count = count+1
    }, 
    
    error = function(e) { 
      message(paste0("Error for ", gene))
      skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }    
  
}


# read the rds back in for induced genes ----
count = 1
for (i in genes_fit){
  tryCatch( 
  {
      print(i)
      gene = i
  
  fit = readRDS(paste0("./fit_mean_modv3_p38input/", gene, ".rds"))
  fitpar = data.frame(t(fit$par))
  fitpar$gene = i
  
  if(count == 1){
    collect.fitpars = fitpar
  }else{
    collect.fitpars = rbind(collect.fitpars, fitpar)
  }
  count = count+1
  }, 
  
  error = function(e) { 
    message(paste0("Error for ", gene))
    skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }   
  
}




# write.table(collect.fitpars, "./fit_mean_modv3_p38input/collect.fitpars_modelv3_p38input_ksynfix_500genes_Dec2020.txt",quote = F, sep = "\t", row.names = F)
collect.fitpars.m = melt(collect.fitpars, id.vars = "gene")
ggplot(collect.fitpars.m[!grepl("mRNA|stim|kt", collect.fitpars.m$variable),], aes(x = variable, y = as.numeric(value))) +
  geom_boxplot(outlier.shape = NA)+geom_point(position = "jitter")+facet_wrap(~ variable, drop = TRUE, scales = "free")
ggplot(collect.fitpars.m[!(grepl("mRNA|stim|kt", collect.fitpars.m$variable)|grepl("Il10|Socs1|Cxcl11|Nkg7", collect.fitpars.m$gene)),
                         ], aes(x = variable, y = as.numeric(value))) +
  geom_boxplot(outlier.shape = NA)+geom_point(position = "jitter")+
  facet_wrap(~ variable, drop = TRUE, scales = "free")
ggplot(collect.fitpars.m[!(grepl("mRNA|stim|k_|k0|tau|n", collect.fitpars.m$variable)|grepl("Il10|Socs1|Cxcl11|Nkg7|Ccr7", collect.fitpars.m$gene)),
                         ], aes(x = variable, y = as.numeric(value))) +
  geom_boxplot(outlier.shape = NA)+geom_point(position = "jitter")+
  facet_wrap(~ variable, drop = TRUE, scales = "free_x")

# solve the model for all genes with fitted params-------------------------------------
collect.fitpars = read.delim("./fit_mean_modv3_p38input/collect.fitpars_modelv3_p38input_ksynfix_500genes_Dec2020.txt")
for (i in seq(1:nrow(collect.fitpars))){
  
  print(collect.fitpars$gene[i])
  gene = collect.fitpars$gene[i]
  pars.fitted = unlist(collect.fitpars[i, -ncol(collect.fitpars)])
  optim <- solve_model(pars.fitted)
  optim$gene = gene
  
  cost.fitted = data.frame(cost = cost(pars.fitted)$model)
  cost.fitted$gene = gene
  
  if(i == 1){
    collect.optim = optim
    collect.cost = cost.fitted
  }else{
    collect.optim = rbind(collect.optim, optim)
    collect.cost = rbind(collect.cost, cost.fitted)
  }
  
}

# plot cost and simulation as heatmap-----------------------------------------------------------
rownames(collect.cost) = collect.cost$gene
hist(collect.cost$cost, breaks = 100)
pheatmap(collect.cost[,1,drop=FALSE] , cluster_rows = F, cluster_cols = F, show_rownames = T)
# write.table(collect.cost, "./fit_mean_modv3_p38input/collect.cost_fittedmodelv3_p38input_ksynfix_500genes_Dec2020.txt", quote = F,row.names = F, sep = "\t")

#simplot----
collect.optim$timeround = round(collect.optim$time)
tmp = collect.optim[(collect.optim$timeround %in% c(1, 30, 60, 180, 301, 480)), ] 
tmp$stim_time = paste0(tmp$label, "_", tmp$timeround)
collect.optim.cast = dcast(tmp, formula = gene~stim_time, value.var = "mRNA")

rownames(collect.optim.cast) = collect.optim.cast$gene
p=pheatmap(collect.optim.cast[, 
                              c(2,4,7,3,5,6, 
                                8,10,13,9,11,12,
                                14,16,19,15,17,18,
                                20,22,25,21,23,24,
                                26,28,31,27,29,30,
                                32,34,37,33,35,36)], scale = "row",
           gaps_col = c(6,12,18,24,30), clustering_method = "ward.D2",
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
           breaks=c(-4,seq(-1.5,1.5,length=100),4), 
           cluster_rows =T, cluster_cols = F, cutree_rows = 4,show_rownames = T)

#simplot with 15 mins----
tmp = collect.optim[(collect.optim$timeround %in% c(1, 15, 30, 60, 180, 301, 480)), ] 
tmp$stim_time = as.factor(paste0(tmp$label, "_", tmp$timeround))
collect.optim.cast = dcast(tmp, formula = gene~stim_time, value.var = "mRNA")

rownames(collect.optim.cast) = collect.optim.cast$gene
p=pheatmap(collect.optim.cast[,
                              c(as.numeric(tmp$stim_time)+1)[1:42]] , scale = "row",
           gaps_col = c(7,14,21,28,35), clustering_method = "ward.D2",
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
           breaks=c(-4,seq(-1.5,1.5,length=100),4), 
           cluster_rows =T, cluster_cols = F, cutree_rows = 4,show_rownames = T)


#plot actual data in same order
mat = as.data.frame(t(my.dataframe.selectmean[, -1]))
mat = as.data.frame(sapply(mat,as.numeric))
colnames(mat) = my.dataframe.selectmean$Group.1
rownames(mat) = colnames(my.dataframe.selectmean)[-1]
mat = mat[rownames(mat) %in% collect.optim.cast$gene,]
# mat = mat[!grepl("Il10", rownames(mat)),]
pheatmap(mat, scale = "row", cluster_rows = F, cluster_cols = F,
         gaps_col = c(5,10,15,20,25), clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1.5,1.5,length=100),4))

gene_order = rownames(collect.optim.cast)[p$tree_row$order]
pheatmap(mat[gene_order, ], scale = "row", cluster_rows = F, cluster_cols = F,
         # gaps_col = c(5,10,15,20,25),
         clustering_method = "ward.D2", cutree_rows = 3,
         # gaps_row = c(16,28,53),
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1.5,1.5,length=100),4))
pheatmap(mat[gene_order, c(35,1:6, 35,7:10, 35,11:16, 35,17:22, 35,23:28, 35,29:34)], scale = "row", cluster_rows = F, cluster_cols = F,
         gaps_col = c(7,12,19,26,33), clustering_method = "ward.D2", cutree_rows = 3,
         # gaps_row = c(150, 223,386),
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1.5,1.5,length=100),4))

collect.cost2 = collect.cost
pheatmap(collect.cost2[gene_order,1,drop=FALSE] , 
         # gaps_row = c(40, 77,86),
         cluster_rows = F, cluster_cols = F, show_rownames = T)

# plot single gene simulation ---------------------------------------------
mat.m = melt(cbind(gene = rownames(mat), mat), id.vars = "gene")
mat.m$time = gsub("..*_", "",mat.m$variable)
mat.m$time = as.numeric(gsub("hr", "",mat.m$time))*60
mat.m$label = gsub("_..*", "", mat.m$variable)
colnames(mat.m) = c("gene", "variable","mRNA", "time", "label")

ggplot(mat.m[mat.m$gene=="Nos2",], aes(time, mRNA/max(mRNA), color = label))+
  geom_line(data = collect.optim[collect.optim$gene=="Il6",] )+
  geom_point() +
  xlim(-100, 500)+# ylim(0,1.2)+
  facet_wrap(~label, nrow = 1)
