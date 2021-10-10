# calculating distances between single cell distributions

# info theory scripts
# install.packages("backports")
# install.packages("devtools")
# install.packages("Rtools")

# library(devtools)
# install_github("sysbiosig/SLEMI")
# devtools::install_github("tylermorganwall/rayshader")
# install.packages('./infotheo/mi-by-decoding-master/R/release/MIdecoding_1.0-1.tar.gz', repos = NULL, type="source")
library(SLEMI);library(ggplot2);library(ggrepel);library(ksheu.library1);library(Seurat);library(pheatmap);library(R.matlab)
library(plotly);library(rayshader)
# help(SLEMI)
# example(MIdecoding)

setwd("F:/scRNAseq_macro/scRNAseq_macro/")
# signaling input data--------------------------
x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_LPS100ng_756_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_ikbamut_10ngTNF_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_P3CSK4100ng_547_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_ikbamut_10ngTNF_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_TNF10ng_754_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_ikbamut_10ngTNF_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_LPS100ng_756_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_TNF10ng_754_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_CpG100nM_610_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_PIC50ug_566_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_CpG100nM_610_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_P3CSK4100ng_547_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_TNF10ng_754_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_P3CSK4100ng_547_smoothed.mat")

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_TNF10ng_754_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_PIC50ug_566_smoothed.mat")

names = c("LPS100ng_756","TNF10ng_754","PIC50ug_566","P3CSK4100ng_547")
collect = data.frame()
for (j in seq(1:(length(names)-1) )){
  for (k in seq(j+1,length(names)) ){
    
    stim1 = names[j]
    stim2 = names[k]
  print(paste0(stim1,"&",stim2))
  x.1 = readMat(paste0("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_",stim1,"_smoothed.mat"))
  x.2 = readMat(paste0("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_",stim2,"_smoothed.mat"))

if(1){
  mat.1 = data.frame(x.1$data); 
  colnames(mat.1) = paste0("X",seq(0, (ncol(mat.1)-1)*5, 5)) #timepoints every 5mins
  order = rownames(mat.1)[order(apply(mat.1[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
  mat.1 = mat.1[match(order, rownames(mat.1)), ]
  mat.1 = cbind(label = stim1, mat.1)
  pheatmap(mat.1[,-1], cluster_cols = F, cluster_rows = F, show_rownames = F,
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
           breaks=c(0,seq(0.1,6.9,length=100),7))
  
  mat.2 = data.frame(x.2$data);
  colnames( mat.2) = paste0("X",seq(0, (ncol( mat.2)-1)*5, 5)) #timepoints every 5mins
  order = rownames( mat.2)[order(apply( mat.2[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
  mat.2 =  mat.2[match(order, rownames( mat.2)), ]
  mat.2 = cbind(label = stim2,  mat.2)
  pheatmap(mat.2[,-1], cluster_cols = F, cluster_rows = F, show_rownames = F,
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
           breaks=c(0,seq(0.1,4.9,length=100),5))
}

for (i in 3:98){
  
  my.dataframe = rbind(mat.1, mat.2)
  my.dataframe$label = as.factor(my.dataframe$label)
  my.dataframe = na.omit(my.dataframe[,c(1,i)]) #for timepoint, for timeseries use [,c(1,3:i)]; [,c(1,i)]
  output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                          # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                          testing = F, boot_num = 10, boot_prob = .8 ,testing_cores = 4) #without0hr
  
  tmp = data.frame(time = (i-1)*5, pair = paste0(stim1,"_",stim2), cc = output_capacity$cc)
  collect = rbind(collect, tmp)
}
ggplot(collect[grepl(paste0(stim1,"..*",stim2), collect$pair),], aes(time/60, cc))+geom_point()+
  theme_bw(base_size = 12)+ylim(c(0,1))+xlab("Time (h)")+ylab("channel capacity")+ggtitle(paste0(stim1,"VS",stim2))
# ggsave(paste0("F://scRNAseq_macro/scRNAseq_macro/images/NFkBsignaling", stim1,"VS",stim2, "timepoint.png"))
  }
}
# write.table(collect, "F://scRNAseq_macro/scRNAseq_macro/infotheo/SLEMI_NFkBsignaling_pairwise_timeseries.txt",quote=F,sep="\t",row.names = F)
# write.table(collect, "F://scRNAseq_macro/scRNAseq_macro/infotheo/SLEMI_NFkBsignaling_pairwise_timepoint.txt",quote=F,sep="\t",row.names = F)

collect = read.delim("F://scRNAseq_macro/scRNAseq_macro/infotheo/SLEMI_NFkBsignaling_pairwise_timeseries.txt")
ggplot(collect[grepl("LPS..*PIC", collect$pair),], aes(time/60, cc))+geom_point()+
  theme_bw(base_size = 12)+ylim(c(0,1))+xlim(0,5)+ xlab("Time (h)")+ylab("max MI (bits)")+
  geom_vline(xintercept = 3, color = "gold", size = 1,linetype = "dashed")

# signaling codon vs gene expression------
# plot features----
expt.info = readxl::read_excel("F://enhancer_dynamics/nfkb_trajectories_20190302/expt_info.xlsx")
features = read.delim("F://enhancer_dynamics/features.txt", sep = ",")
ggplot(features[grepl("610|547|756|566|754|783", features$ID), ], aes(as.factor(ID),off_times))+geom_boxplot(outlier.shape =NA)+
  geom_point(position = "jitter", aes(color = as.factor(ID)))+theme_bw(base_size = 18)
ggplot(features[grepl("610|547|756|566|754|783", features$ID), ], aes(as.factor(ID),oscpower))+geom_boxplot(outlier.shape =NA)+
  geom_point(position = "jitter", aes(color = as.factor(ID)))+theme_bw(base_size = 18)
ggplot(features[grepl("610|547|756|566|754|783", features$ID), ], aes(as.factor(ID),derivatives_002))+geom_boxplot(outlier.shape =NA)+
  geom_point(position = "jitter", aes(color = as.factor(ID)))+theme_bw(base_size = 18)

feature.names = colnames(features)[-c(1,2)]
features.wanted = features[grepl("610|547|756|566|754|783", features$ID), ]

features.wanted = readxl::read_excel("C://Users/kathe/Downloads/data_table_KS.xlsx", sheet = 3)

#take medians of features to codewords----
k <- which(is.na(features.wanted), arr.ind=TRUE)
features.wanted.median = (features.wanted)
features.wanted.median[k] <-0
# features.wanted.median[k] <- colMeans(features.wanted.median[,-1], na.rm=TRUE)[features.wanted.median[,1]]

features.wanted.median$speed = rowMeans(features.wanted.median[,c("derivatives_002", "max_pk1_speed", "pk1_time")])
features.wanted.median$peakAmp = rowMeans(features.wanted.median[,c("fold_change", "peak2peak", "pk1_amp")])
features.wanted.median$oscVSnonosc = (features.wanted.median[,c("oscpower")])
features.wanted.median$totalactivity = (features.wanted.median[,c("max_integral")])
features.wanted.median$duration = rowMeans(features.wanted.median[,c("duration_002", "num_peaks")])
features.wanted.median$earlyLate = (features.wanted.median[,c("time2HalfMaxIntegral")])

features.codons = features.wanted.median[, c("ID","Name","speed","peakAmp","oscVSnonosc","totalactivity",
                                      "duration","earlyLate")]
features.codons$Name= ifelse(grepl("P3C", features.codons$Name), "P3CSK",
                             ifelse(grepl("LPS", features.codons$Name), "LPS",
                                    ifelse(grepl("PIC", features.codons$Name), "PIC",
                                           ifelse(grepl("TNF", features.codons$Name), "TNF",
                                                  ifelse(grepl("CpG", features.codons$Name), "CpG","Unstim")
                                           ))))

# do the same for IKbaMM and WT--------------
features.wanted = readxl::read_excel("C://Users/kathe/Downloads/data_table_KS.xlsx", sheet = 2) #WT zscore
features.wanted = readxl::read_excel("C://Users/kathe/Downloads/data_table_KS.xlsx", sheet = 4) #MM zscore

#take medians of features to codewords----
features.wanted.median = (features.wanted)

features.wanted.median$speed = rowMeans(features.wanted.median[,c("derivatives_002", "max_pk1_speed", "pk1_time")])
features.wanted.median$peakAmp = rowMeans(features.wanted.median[,c("fold_change",  "pk1_amp")])
features.wanted.median$oscVSnonosc = rowMeans(features.wanted.median[,c("oscpower", "peak2peak","num_peaks")])
features.wanted.median$totalactivity = rowMeans(features.wanted.median[,c("max_integral")])
features.wanted.median$duration = rowMeans(features.wanted.median[,c("duration_002" )])
features.wanted.median$earlyLate = rowMeans(features.wanted.median[,c("time2HalfMaxIntegral")])

features.codons = features.wanted.median[, c("Ligand","speed","peakAmp","oscVSnonosc","totalactivity",
                                             "duration","earlyLate")]
features.codons$Name= ifelse(grepl("P3C", features.codons$Ligand), "P3CSK",
                             ifelse(grepl("LPS", features.codons$Ligand), "LPS",
                                    ifelse(grepl("PolyIC", features.codons$Ligand), "PIC",
                                           ifelse(grepl("TNF", features.codons$Ligand), "TNF",
                                                  ifelse(grepl("CpG", features.codons$Ligand), "CpG","Unstim")
                                           ))))
#calc LPT MI for codons---- 
collect = data.frame()

for (codon in c("speed","peakAmp","oscVSnonosc","totalactivity","duration","earlyLate")){
  
  skip_to_next <- FALSE
  print(codon)
  
  tryCatch(
    {
      print(codon)
      data = features.codons
      # data = data[grep(paste0(list[k],"|",list[j]), data$Name),]
      
      #--------------------------------mi using SLEMI -----------------------
      
      output_capacity <- capacity_logreg_main(data, signal = "Name", 
                                              response = codon, #c("speed","peakAmp","oscVSnonosc","totalactivity","duration","earlyLate"),#codon,
                                              testing=T, boot_prob = 0.5, boot_num = 30, testing_cores = 4)
      sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
      
      
      tmp = data.frame(time = "MM", cc = output_capacity$cc, sd = sd, codon = codon) #"all")#
      collect = rbind(collect, tmp)
    }, error = function(e) { skip_to_next <<- TRUE; print("error")})
  
  if(skip_to_next) { next }  
}
# write.table(collect, "./infotheo/SLEMI_3stim_NFkBsignalingcodons_IkBaWT_zscore.txt",quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_3stim_NFkBsignalingcodons_IkBaMM_zscore.txt",quote = F, sep = "\t",row.names = F)


#calc pairwise for codons----  
collect = data.frame()
list = c("CpG","P3CSK","LPS", "TNF", "PIC", "Unstim") #
list = c("LPS", "TNF", "PIC") #

for (j in seq(1:(length(list)-1) )){
  for (k in seq(j+1,length(list)) ){
    
    for (codon in c("speed","peakAmp","oscVSnonosc","totalactivity","duration","earlyLate")){
      
      skip_to_next <- FALSE
      print(codon)
      
      tryCatch(
        {
          
          
          print(list[j])
          print(list[k])
          
          data = features.codons
          data = data[grep(paste0(list[k],"|",list[j]), data$Name),]
          
          #--------------------------------mi using SLEMI -----------------------
          
          output_capacity <- capacity_logreg_main(data, signal = "Name", 
                                                  response = codon, #c("speed","peakAmp","oscVSnonosc","totalactivity","duration","earlyLate"),#codon,
                                                  testing=T, boot_prob = 0.5, boot_num = 30, testing_cores = 4)
          sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
          
          
          tmp = data.frame(time = "MM", stim1 = list[j], stim2 = list[k], cc = output_capacity$cc, sd = sd, codon = codon) #"all")#
          collect = rbind(collect, tmp)
          
          # output_mi  <- mi_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
          #                              paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/mi_", i),
          #                              pinput=rep(1/6,6))
          
        }, error = function(e) { skip_to_next <<- TRUE; print("error")})
      
      if(skip_to_next) { next }  
      
      
     }
  }
}
collect$pair = paste0(collect$stim1,"_", collect$stim2)
# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons.txt",quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_all6.txt",quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaMM.txt",quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_all6_IkBaMM.txt",quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaWT_zscore.txt",quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaMM_zscore.txt",quote = F, sep = "\t",row.names = F)

collect$pair = factor(collect$pair, levels= c("CpG_LPS","CpG_P3CSK", "P3CSK_LPS",
                                              "CpG_TNF","CpG_Unstim","LPS_TNF","LPS_Unstim", "P3CSK_TNF","P3CSK_Unstim",
                                              "CpG_PIC","LPS_PIC","P3CSK_PIC",
                                              "TNF_PIC", "PIC_Unstim",
                                              "TNF_Unstim"))
collect$pair = factor(collect$pair, levels= c("CpG_Unstim", 
                                              "LPS_Unstim", "P3CSK_Unstim","PIC_Unstim","TNF_Unstim",
                                              "CpG_LPS","CpG_P3CSK","CpG_PIC","CpG_TNF",
                                              "P3CSK_LPS","LPS_PIC","LPS_TNF",
                                              "P3CSK_PIC","P3CSK_TNF",
                                              "TNF_PIC"))
collect$color = as.factor(ifelse(grepl("Unstim",collect$pair), "Unstim", 
                                 ifelse(grepl("CpG",collect$pair), "CpG", 
                                        ifelse(grepl("LPS",collect$pair), "LPS", 
                                               ifelse(grepl("P3CSK",collect$pair), "P3CSK", 
                                                      ifelse(grepl("PIC",collect$pair), "PIC",
                                                             ifelse(grepl("TNF",collect$pair), "TNF",0)
                                                      )))) ))
colors_list = c("Unstim" = "gray", "CpG"="#F8766D", "IFNb"="#B79F00","LPS"= "#00BA38","P3CSK"= "#00BFC4","PIC"= "#619CFF","TNF"= "#F564E3")
ggplot(collect[!grepl("IFNb", collect$pair),], aes((cc)))+geom_density(aes(group = pair, fill = color), alpha = 0.3)+
  scale_fill_manual(values = colors_list)+
  facet_wrap(~pair, ncol = 4, scales = "free_y")+theme_bw()

ggplot(collect[!grepl("blah", collect$pair),], aes(pair, cc))+
  geom_bar(aes(fill = color), stat = "identity", position = "nudge",alpha = 0.7)+
  scale_fill_manual(values = colors_list)+
  facet_wrap(~codon, ncol = 3, scales = "fixed")+theme_bw()+ylab("max MI")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(collect[!grepl("blah", collect$pair),], aes(codon, cc))+
  geom_bar(aes(fill = color), stat = "identity", position = "nudge",alpha = 0.7)+
  scale_fill_manual(values = colors_list)+
  facet_wrap(~pair, ncol = 5, scales = "fixed")+theme_bw()+ylab("max MI")
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

if(1){
  collect = rbind(read.delim("./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaWT_zscore.txt"))
                  # read.delim( "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_all6.txt"))
  collect.cast_codon = dcast(collect, pair~codon, value.var = "cc")
  # collect_genes = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt")
  collect_genes = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_WT_ISnorm_allgenes.txt")
}
if(1){
  collect = rbind(read.delim("./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaMM_zscore.txt"))
                  # read.delim( "./infotheo/SLEMI_pairwise_NFkBsignalingcodons_all6_IkBaMM.txt"))
  collect.cast_codon = dcast(collect, pair~codon, value.var = "cc")
  collect_genes = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_MM_ISnorm_allgenes.txt")
}

if(0){
  collect = rbind(read.delim("./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaWT_zscore.txt"),
                  read.delim("./infotheo/SLEMI_pairwise_NFkBsignalingcodons_IkBaMM_zscore.txt"))
                  
  ggplot(collect, aes(factor(time, levels= c("WT","MM")), cc))+ geom_bar(stat= "identity", aes(fill = time))+
    geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
    facet_grid(pair~codon)+theme_bw(base_size = 14)+ylab('max MI (bits)')+xlab("")
  
  collect.cast_codon = dcast(collect, pair~codon, value.var = "cc")
  collect_genes = rbind(read.delim("./infotheo/channel_capacity_pairwise_BMDM2_WT_ISnorm_allgenes.txt"),
                        read.delim("./infotheo/channel_capacity_pairwise_BMDM2_MM_ISnorm_allgenes.txt"))
  collect_genes$pair = paste0(collect_genes$stim1,"_", collect_genes$stim2)
  rownames(collect_genes) = collect_genes$pair
  ggplot(collect_genes[grepl("Socs3",collect_genes$gene),], aes(time, cc))+
    geom_bar(stat= "identity", aes(fill = time))+facet_grid(~pair)
  
}
collect_genes$pair = paste0(collect_genes$stim1,"_", collect_genes$stim2)
library(reshape2)
collect.cast = dcast(collect_genes, pair~gene, value.var = "cc")
collect.all = cbind(collect.cast_codon, collect.cast[match(collect.cast_codon$pair, (collect.cast$pair)),-1])
collect.all = collect.all[!grepl("CpG|Unstim",collect.all$pair),]

rownames(collect.all) = collect.all$pair
ggplot(collect.all, aes(pair, Cxcl10))+geom_bar(stat= "identity")
ggplot(collect.all, aes(oscVSnonosc, collect.all[,"Cxcl10"], label = pair))+geom_point(size = 5)+theme_bw()+geom_text_repel()
ggplot(collect.all, aes(totalactivity, collect.all[,"Tlr2"], label = pair))+geom_point(size = 5)+theme_bw()+geom_text_repel()
ggplot(collect.all, aes(speed, collect.all[,"Tlr2"], label = pair))+geom_point(size = 5)+theme_bw()+geom_text_repel()

corrs.all = data.frame()
for (i in 8:499){
  
  skip_to_next <- FALSE

  
  tryCatch(
    {
  print(colnames(collect.all)[i])
  corrs = data.frame(
    gene = colnames(collect.all)[i],
  osc = cor.test(collect.all$oscVSnonosc, collect.all[,i], method = "pearson")$estimate,
  speed =cor.test(collect.all$speed, collect.all[,i], method = "pearson")$estimate,
  peakAmp=cor.test(collect.all$peakAmp, collect.all[,i], method = "pearson")$estimate,
  earlyLate=cor.test(collect.all$earlyLate, collect.all[,i], method = "pearson")$estimate,
  duration=cor.test(collect.all$duration, collect.all[,i], method = "pearson")$estimate,
  totalactivity=cor.test(collect.all$totalactivity, collect.all[,i], method = "pearson")$estimate,
  # all6codons=cor.test(collect.all$all, collect.all[,i], method = "pearson")$estimate,
  osc.pval = cor.test(collect.all$oscVSnonosc, collect.all[,i], method = "pearson")$p.value,
  speed.pval =cor.test(collect.all$speed, collect.all[,i], method = "pearson")$p.value,
  peakAmp.pval=cor.test(collect.all$peakAmp, collect.all[,i], method = "pearson")$p.value,
  earlyLate.pval=cor.test(collect.all$earlyLate, collect.all[,i], method = "pearson")$p.value,
  duration.pval=cor.test(collect.all$duration, collect.all[,i], method = "pearson")$p.value,
  totalactivity.pval=cor.test(collect.all$totalactivity, collect.all[,i], method = "pearson")$p.value
  # all6codons.pval=cor.test(collect.all$all, collect.all[,i], method = "pearson")$p.value
  
  )
  corrs.all = rbind(corrs.all, corrs)
    }, error = function(e) { skip_to_next <<- TRUE; print("error")})
  
  if(skip_to_next) { next }  
}

rownames(corrs.all) = corrs.all$gene
clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
corrs.all$clusters = clusters$manual_fix[match(rownames(corrs.all), clusters$gene)]
corrs.all$clusters = gsub("mB","mC", corrs.all$clusters)
corrs.all$clusters = ifelse(corrs.all$clusters=="mA","AP1",
                                 ifelse(corrs.all$clusters=="mC","NFkB",
                                        ifelse(corrs.all$clusters=="mD","NFkB&p38",
                                               ifelse(corrs.all$clusters=="mE","NFkB|IRF","IRF")))) 
corrs.all$clusters = factor(corrs.all$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))                                  
pheatmap(corrs.all[corrs.all$clusters=="NFkB",c(2:7)], scale = "none", cluster_rows = T, cluster_cols = F)
         # colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         # breaks=c(-2,seq(-1,1,length=100),2)
         
# corrs.mm = corrs.all
# corrs.wt = corrs.all
corrs.wt.pval = corrs.wt[apply(corrs.wt[,c(8:13)]<1, 1, any),] 
corrs.mm.pval = corrs.mm[apply(corrs.mm[,c(8:13)]<1, 1, any),] 
corrs.all = cbind(wt= corrs.wt.pval[,c(2:7)], mm=corrs.mm.pval[match(corrs.wt.pval$gene, corrs.mm.pval$gene),c(2:7,14)])
# write.table(corrs.all, "./infotheo/SLEMI_pairwise_signaling2geneexpression_WTMM.corrs.all.txt", row.names = T, quote = F, sep="\t")

corrs.wt.pval = corrs.wt[apply(corrs.wt[,c(8:13)]<0.25, 1, any),] 
corrs.mm.pval = corrs.mm[apply(corrs.mm[,c(8:13)]<0.25, 1, any),] 
corrs.all = cbind(wt= corrs.wt.pval[,c(2:7)], mm=corrs.mm.pval[match(corrs.wt.pval$gene, corrs.mm.pval$gene),c(2:7,14)])
test = corrs.all[corrs.all$mm.clusters=="NFkB",]
test = test[rowSums( abs(test[,c(1:12)]))>0.25, ] #corrs across all sum to >1.5
p=pheatmap(na.omit(test[,c(1,7,2,8,3,9,4,10,5,11,6,12)]), cluster_rows = T, cluster_cols = F, #clustering_method = "ward.D2",
         gaps_col = c(2,4,6,8,10), cutree_rows = 4)


corrs.all.m = melt(corrs.all)
ggplot(corrs.all.m[!grepl('pval', corrs.all.m$variable)&grepl("^NFkB$", corrs.all.m$clusters),], aes(clusters, value))+
  geom_boxplot()+geom_violin()+geom_point(position = "jitter")+
  facet_wrap(~variable, scales = "free")

##################################################################################
# input data each timept individually----
collect = data.frame()
list = c(500, 1000, 2000, 3000, 4000, 5000)
# for (i in c("0.25hr","1hr","3hr","8hr")){ #c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") #0.5 and 5 have no IFNb
for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i, "_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
  macro = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
  macro = subset(macro, subset = timept == i)
  
  if(0){ #include unstim
    # macro = readRDS(paste0("./output/macrophage_tutorial3_subsample_8hr_sctransform.rds"))
    macro.0hr = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_0hr_DBEC.rds"))
    macro = merge(macro.0hr, macro)
  }
  
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  # data = data[ ,grepl("LPS|TNF|PIC", colnames(data))]
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  str(my.dataframe)
  
  
  if(0){ #subsample #of cells
    for (n_cells in list){
      print(n_cells)
      set.seed(1)
      my.dataframe.subset = my.dataframe[sample(rownames(my.dataframe), n_cells),]
      output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", response = colnames(my.dataframe)[-1],
                                              # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cell_subsets/cc_with0hr_",n_cells,"_", i),
                                              output_path = NULL,
                                              testing = F)
      tmp = data.frame(time = i, num_cells = n_cells, cc = output_capacity$cc)
      collect = rbind(collect, tmp)
    }
  }
  #--------------------------------mi using SLEMI -----------------------
  
  #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
  output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                          # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                          testing = F, boot_num = 10, boot_prob = .8 ,testing_cores = 4) #without0hr
  tmp = data.frame(time = i, num_cells = "all", cc = output_capacity$cc)
  collect = rbind(collect, tmp)
  
  # output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
  #                                         paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_with0hr_", i), testing = F)

  # output_mi  <- mi_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1], 
  #                              paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/mi_", i),
  #                              pinput=rep(1/6,6))
} 

# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M0_6stim_500genes.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M0_rep2only_6stim_500genes.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M1_IFNg_6stim_500genes.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M2_IL4_gt80_6stim_500genes.txt", sep = "\t",quote = F, row.names = F)

collect = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_M0_rep2only_6stim_500genes.txt")
ggplot(collect[grepl("all", collect$num_cells),], aes(time, cc, group = num_cells))+geom_point(size=2)+
  geom_hline(yintercept = log10(6)/log10(2), linetype = "dashed")+ ylab("channel capacity")+
  geom_line(aes(group = num_cells, color = num_cells), size = 1)+theme_bw(base_size = 16) +ylim(0,2.6)

ggplot(cc.with0hr, aes(time, channel_capacity, group = set))+geom_point(size=2)+
  geom_line(aes(group = set, color = set), size = 1)+theme_bw(base_size = 16) +xlim(0,8)+ylim(1.5,3)


# based on WGCNA modules-------------------------------
collect = data.frame()
list = c(500, 1000, 2000, 3000, 4000, 5000)
wgcna.sets = read.delim("./wgcna/WGCNA_5modules_20210121.txt")
wgcna.sets = read.delim("./wgcna/WGCNA_3modules_20210410.txt")
wgcna.sets = split.data.frame(wgcna.sets,wgcna.sets$net.colors)
# for (i in c("0.25hr","1hr","3hr","8hr")){ #c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") #0.5 and 5 have no IFNb
# for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
for (i in c("PM_B6.LFD","PM_B6.old","PM_B6.HFD")){
  
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i, "_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
  # macro = subset(macro, subset = timept == i)
  
  macro = readRDS(paste0("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds"))
  macro = subset(macro, subset = type == i)
  
  if(0){ #include unstim
    # macro = readRDS(paste0("./output/macrophage_tutorial3_subsample_8hr_sctransform.rds"))
    macro.0hr = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_0hr_DBEC.rds"))
    macro = merge(macro.0hr, macro)
  }
  
  for (geneset in wgcna.sets){
    print(geneset)
    genes = geneset$gene
    data = macro[["ISnorm"]]@data
    data = data.frame(data)
    meta = macro@meta.data
    colnames(data) = paste0(meta$stimulus, "_",meta$timept)
    # data = data[ ,grepl("LPS|TNF|PIC", colnames(data))]
    data.select = data[rownames(data) %in% genes,]
    my.dataframe = cbind(label = colnames(data.select), data.frame(t(data.select)))
    my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
    rownames(my.dataframe) = seq(1:nrow(my.dataframe))
    str(my.dataframe)
    
    
    if(0){ #subsample #of cells
      for (n_cells in list){
        print(n_cells)
        set.seed(1)
        my.dataframe.subset = my.dataframe[sample(rownames(my.dataframe), n_cells),]
        output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", response = colnames(my.dataframe)[-1],
                                                # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cell_subsets/cc_with0hr_",n_cells,"_", i),
                                                output_path = NULL,
                                                testing = F)
        tmp = data.frame(time = i, num_cells = n_cells, cc = output_capacity$cc)
        collect = rbind(collect, tmp)
      }
    }
    #--------------------------------mi using SLEMI -----------------------
    
    #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
    output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                            # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                            testing = F, boot_num = 10, boot_prob = .8 ,testing_cores = 4) #without0hr
    tmp = data.frame(time = i, num_cells = paste0("all_ME",geneset$net.colors), cc = output_capacity$cc)
    collect = rbind(collect, tmp)
    
    # output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
    #                                         paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_with0hr_", i), testing = F)
    
    # output_mi  <- mi_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1], 
    #                              paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/mi_", i),
    #                              pinput=rep(1/6,6))
    
  }
}
collect = collect[!duplicated(collect),]
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M0_6stim_500genes.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M0_rep2only_6stim_WGCNAgenesets.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M1_IFNg_6stim_WGCNAgenesets.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_M2_IL4_gt80_6stim_WGCNAgenesets.txt", sep = "\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/cell_subsets/collect_cc_without0hr_PMexpt_Feb2021_rmUnstim_5stim_WGCNAgenesets.txt", sep = "\t",quote = F, row.names = F)

collect_allgenes = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_M0_rep2only_6stim_500genes.txt")
collect = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_M0_rep2only_6stim_WGCNAgenesets.txt")
collect = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_M1_IFNg_6stim_WGCNAgenesets.txt")
collect = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_M2_IL4_gt80_6stim_WGCNAgenesets.txt")
collect = rbind(collect_allgenes, collect)
ggplot(collect[grepl("all", collect$num_cells),], aes(time, cc, group = num_cells))+geom_point(size=2)+
  geom_hline(yintercept = log10(6)/log10(2), linetype = "dashed")+ ylab("channel capacity")+
  geom_line(aes(group = num_cells, color = num_cells), size = 1)+theme_bw(base_size = 16) +ylim(0,2.6)
ggplot(collect[grepl("all", collect$num_cells),], aes(time, cc, group = num_cells))+geom_point(size=2)+
  geom_hline(yintercept = log10(6)/log10(2), linetype = "dashed")+ ylab("channel capacity")+
  scale_color_brewer(palette = "Paired")+
  geom_line(aes(group = num_cells, color = num_cells), size = 2)+theme_bw(base_size = 16) +ylim(0,2.6)

#for PMs
collect = read.delim("./infotheo/cell_subsets/collect_cc_without0hr_PMexpt_Feb2021_rmUnstim_5stim_WGCNAgenesets.txt")
ggplot(collect[grepl("all", collect$num_cells),], aes(time, cc, group = num_cells))+geom_point(size=2)+
  geom_hline(yintercept = log10(6)/log10(2), linetype = "dashed")+ ylab("channel capacity")+
  geom_line(aes(group = num_cells, color = num_cells), size = 1)+theme_bw(base_size = 16) 


# search for high channel gene for 1, 2, 3,... genes----
collect = data.frame()
collect_all = data.frame()
# for (i in c("0.25hr","0.5hr","1hr","3hr","5hr","8hr")){
# for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
# for (i in c("0.25hr","1hr","3hr","8hr")){
# for (i in c("BMDM2_WT", "BMDM2_MM")){
for (i in c("PM_B6.LFD","PM_B6.old","PM_B6.HFD")){
  print(i)
  
  # macro = readRDS("./output/macrophage_M0all_500genes_Dec2020.rds")
  # macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
  # macro = readRDS(paste0("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds"))
  # macro = readRDS("./output/macrophage_M1_IFNg_500genes_DBEC.rds")
  # macro = readRDS("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
  # macro = subset(x=macro, subset = timept==i)
  # macro = subset(macro, subset= stimulus!="IFNb")
  
  # macro = readRDS("output/macrophage_BMDM2_WT_MM_500genes_DBEC.rds")
  
  macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds")
  macro = subset(macro, type == i)
  # macro = subset(x=macro, subset = stimulus=="LPS"|stimulus=="TNF"|stimulus=="PIC"|stimulus=="P3CSK"|stimulus=="IFNb")
  macro = subset(x=macro, subset = stimulus=="LPS"|stimulus=="TNF"|stimulus=="PIC")
  
  if(0){ #include unstim
    # macro = readRDS(paste0("./output/macrophage_tutorial3_subsample_8hr_sctransform.rds"))
    macro.0hr = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_0hr_DBEC.rds"))
    macro = merge(macro.0hr, macro)
  }
  
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  str(my.dataframe)
  
  #--------------------------------mi using SLEMI -----------------------
  # iterate through genes
  store = "blah"
  current_max = 0
  for (gene in colnames(my.dataframe)[-1]){
    
    skip_to_next <- FALSE
    print(gene)
    
    tryCatch(
      {
        #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
        output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = gene,
                                                output_path = NULL, testing = F)
        if (output_capacity$cc > current_max){
          store = gene
          current_max = output_capacity$cc
          print(store)
          print(current_max)
        }
        tmp0 = data.frame(time = "3hr", type = i, cc = output_capacity$cc, gene = gene)
        collect_all = rbind(collect_all, tmp0)
        
      }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }  
    
  }
  tmp = data.frame(time = "3hr", type = i, cc = current_max, gene = store)
  collect = rbind(collect, tmp)
  head(collect)
    
}

# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M0all_gt80_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect, "./infotheo/SLEMI_singlegene_collectmax_M0all_gt80_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M0_rep2only_SCT_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_BMDM2_WTMM_3stim_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_5stim_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_gt120_5stim_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_gt120_3stim_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)

collect_all0 = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_SCT_500genes.txt")
collect_all$ISnorm = collect_all0$cc[match(paste0(collect_all$gene, collect_all$time), 
                                           paste0(collect_all0$gene,collect_all0$time))]
collect_all$diff = collect_all$cc-collect_all$ISnorm
ggplot(collect_all, aes(cc, ISnorm))+geom_point(aes(color = time))


collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_gt120_5stim_ISnorm_500genes.txt")

# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M1_IFNg_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect, "./infotheo/SLEMI_singlegene_collectmax_M1_IFNg_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M1_IFNg_SCT_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_SCT_500genes.txt")
collect_all0 = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_500genes.txt")
collect_all$RNA = collect_all0$cc[match(paste0(collect_all$gene, collect_all$time), 
                                        paste0(collect_all0$gene,collect_all0$time))]
ggplot(collect_all, aes(cc, RNA))+geom_point(aes(color = time))

# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect, "./infotheo/SLEMI_singlegene_collectmax_M2_IL4_gt80_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_SCT_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt")

# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall_with0hr.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect, "./infotheo/SLEMI_singlegene_collectmax_with0hr.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_collectall.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect, "./infotheo/SLEMI_singlegene_collectmax.txt",sep = "\t", quote = F, row.names = F)


#graph it USE THIS
ggplot(collect_all[grepl("3hr", collect_all$time),], aes( cc))+ geom_density()+
  # geom_point(aes(color = time),position = position_jitter(seed = 1))+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+
  geom_vline(xintercept = 0.7, linetype="dotted", size = 1)+
  geom_vline(xintercept = (log10(6)/log10(2)), linetype="dashed", size = 1)

ggplot(collect_all[!grepl("0.5hr|^5hr|x24hr", collect_all$time),], aes(time, cc))+geom_violin()+
  geom_point(aes(color = time),position = position_jitter(seed = 1))+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+ylab("channel capacity")+
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_text_repel(subset(collect_all[!grepl("0.5hr|^5hr|x24hr", collect_all$time),], cc > 0.9), 
                  mapping = aes(label = gene), size = 4, position = position_jitter(seed = 1))

ggplot(collect_all[], aes(time, (cc)))+geom_violin()+
  geom_point(aes(color = time),position = "jitter")+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+ylab("channel capacity")+
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  facet_grid(~type)+
  geom_text_repel(subset(collect_all, cc > 0.8), 
                  mapping = aes(label = gene), size = 4)

collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all= read.delim("./infotheo/SLEMI_singlegene_collectall_BMDM2_WTMM_3stim_ISnorm_500genes.txt")
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_gt120_5stim_ISnorm_500genes.txt")
dcast = dcast(collect_all, gene~time, value.var = "cc")
dcast = dcast(collect_all, gene~type, value.var = "cc")
# wgcna = read.delim("./wgcna/WGCNA_5modules_20210121.txt")
# dcast$geneset = as.factor(wgcna$net.colors[match(dcast$gene, wgcna$gene)])
rownames(dcast) = dcast$gene
ggplot(dcast,aes(PM_B6.LFD, PM_B6.HFD))+geom_point(aes(color = geneset), size = 2)+
  # geom_text_repel(subset(dcast, PM_B6.old > 0.4),mapping = aes(label = gene), size = 4)+
  geom_abline(slope = 1, intercept = 0)+theme_bw(base_size = 18)
ggplot(na.omit(dcast), aes(geneset, PM_B6.LFD-PM_B6.old))+geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color = geneset),position = "jitter", size = 2)+
  theme_bw(base_size = 18)
pheatmap(na.omit(dcast[rowSums(dcast[,-1])>1.5,-1]), scale = "none", 
         colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
         cluster_cols = F, cluster_rows = T)

#plot Fig2 cc----
mutual = readRDS("./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
# genes = unlist(mutual[15,4])
genes = unlist(mutual[41,4])
genes = c("Ifit3","Mx1", "Mx2", "Oasl1","Gsr", "Sod2", "Gclm", "Gsta3")
genes = c(genesets[[4]],genesets[[5]],genesets[[6]])
genes = c(genesets[[3]],genesets[[9]])

# genes = unlist(mutual[67,4])
pheatmap(na.omit(dcast[genes[order(genes)],-1]), scale = "none", clustering_method = "ward.D2", 
         colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
         cluster_cols = F, cluster_rows = F)


p=pheatmap(na.omit(dcast[,-1][rowSums(dcast[,-1] > .75) > 0, ]), scale = "none", clustering_method = "ward.D2", 
         colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
         cluster_cols = F, cluster_rows = T)

macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
data = GetAssayData(object = macro, assay = "ISnorm", slot = "data")
data <- as.data.frame( as.matrix(data))
meta = macro@meta.data
meta$stimulus <- factor(meta$stimulus, levels = c("Unstim", "LPS", "PIC", "IFNb", "P3CSK", "CpG", "TNF"))
meta = meta[order(meta$stimulus), ]
meta = meta[order(meta$timept), ]
# meta = meta[grepl("LPS|P3C",meta$stimulus),]
meta = meta[grepl("3hr",meta$timept),]
count =0
for (i in c(names(table(meta$timept)))) {
  count[i] <- table(meta$timept)[names(table(meta$timept))==i]
}
colseps <- cumsum(count[-1])
col_order = rownames(meta)
colors_list = list(stimulus = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
pheatmap(data[genes,col_order], scale = "row", cluster_rows = T, cluster_cols = F, clustering_method = "ward.D2", 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         gaps_col = colseps[-6],
         breaks=c(-2,seq(-1,1,length=100),2),
         annotation_col = data.frame(macro@meta.data[,c(7,6)]), 
         annotation_colors = colors_list,
         show_colnames = F, show_rownames = T)


row_order = na.omit(dcast$gene[rowSums(dcast[,-1] > .75) > 0])[p$tree_row$order]
####USE THIS
p=pheatmap(na.omit(dcast[(dcast[,4] > .7),4, drop=F]), scale = "none", clustering_method = "ward.D2", 
           colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
           cluster_cols = F, cluster_rows = T, show_rownames = T)
row_order = na.omit(dcast$gene[(dcast[,4] > .7)])[p$tree_row$order]

pheatmap(data[row_order,col_order], scale = "row", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
            colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
            gaps_col = colseps[-6],
            breaks=c(-2,seq(-1,1,length=100),2),
            annotation_col = data.frame(macro@meta.data[,c(7,6)]), 
            annotation_colors = colors_list,
            show_colnames = F, show_rownames = T)


mat = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_MEAN.txt",row.names = 1)
mat = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_VAR.txt",row.names = 1)
mat = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_FF.txt",row.names = 1)
pheatmap(mat[row_order,c(25,1:24)], scale = "row", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         gaps_col = c(1,5,9,13,17,21),
         breaks=c(-2,seq(-1,1,length=100),2),
         show_colnames = T, show_rownames = T)
pheatmap(mat[row_order,c(25,c(1,5,9,13,17,21),c(1,5,9,13,17,21)+1, c(1,5,9,13,17,21)+2,c(1,5,9,13,17,21)+3)], 
         scale = "none", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
         colorRampPalette((brewer.pal(n = 11, name ="YlOrRd"))[2:9])(103),
         gaps_col = c(1,7,13,19),
         breaks=c(0,seq(.01,2,length=100),5),#breaks=c(0,seq(.01,3,length=100),10),#breaks=c(0,seq(1,5,length=100),10),
         show_colnames = T, show_rownames = T)

# mean vs variance for a few genes-------------------------
mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T)
mat.m0.3hr = mat.m0[grepl("3", mat.m0$time),]
ggplot(mat.m0.3hr[grepl("Tnfaip3$", mat.m0.3hr$gene)&grepl("", mat.m0.3hr$stimulus),], aes(mean, fano))+
  geom_point(aes(color=stimulus),size = 5)+theme_bw(base_size = 14)
ggplot(mat.m0.3hr[grepl("Ifit3$", mat.m0.3hr$gene)&grepl("", mat.m0.3hr$stimulus),], aes(mean, fano))+
  geom_point(aes(color=stimulus),size = 5)+theme_bw(base_size = 14)

# mean vs variance contribution to channel capacity------
p=ggplot()
# for (i in c("LPS|TNF", "P3CSK|LPS", "TNF|PIC","P3CSK|TNF","LPS|PIC","CpG|P3CSK","LPS|IFNb","TNF|IFNb")){ #,""
for (i in c("LPS|TNF")){ #,""
  print(i)

  # i = "P3CSK|TNF"
  if(0){
    collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")  #for all stim
    mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T)
    timept ="3"
    pseudobulk = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_MEAN.txt")
    
  }
  
  if(0){
    collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt")
    mat.m0 = read.table("./output/macrophage_M1_IFNg_500genes_calculatedistfeatures.txt", header = T)
    timept ="3"
  }
  
  if (1){
    collect_all = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt")
    mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T)
    timept ="3"
    pseudobulk = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_MEAN.txt")
    
  }
  if (0){
    collect_all = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes_8hr.txt")
    mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T)
    timept ="8"
    pseudobulk = read.delim("./output/macrophage_M0_rep2only_500genes_pseudobulk_MEAN.txt")
  }
   
  if(0){
    collect_all = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_WT_ISnorm_allgenes.txt")
    collect_all$time="8hr"
    mat.m0 = read.table("./output/macrophage_BMDM2_WT_500genes_calculatedistfeatures.txt", header = T)
    timept= "8"
    pseudobulk = read.delim("./output/macrophage_BMDM2_WT_500genes_pseudobulk_MEAN.txt")
    
  }
  
  if(0){
    collect_all = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_MM_ISnorm_allgenes.txt")
    collect_all$time="8hr"
    mat.m0 = read.table("./output/macrophage_BMDM2_MM_500genes_calculatedistfeatures.txt", header = T)
    timept= "8"
  }
  if(0){ #WT-MM
    collect_all = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_WT_ISnorm_allgenes.txt")
    collect_all2 = read.delim("./infotheo/channel_capacity_pairwise_BMDM2_MM_ISnorm_allgenes.txt")
    collect_all$time="8hr"
    collect_all$WT.cc = collect_all$cc
    collect_all$MM.cc = collect_all2$cc[match(paste0(collect_all$gene,collect_all$stim1,collect_all$stim2), 
                                              paste0(collect_all2$gene,collect_all2$stim1,collect_all2$stim2))]
    collect_all$cc = collect_all$cc - collect_all$MM.cc
    mat.m0 = read.table("./output/macrophage_BMDM2_MM_500genes_calculatedistfeatures.txt", header = T)
    timept= "8"
  }
  collect_all$pair = paste0(collect_all$stim1,"_", collect_all$stim2)
  pair = gsub("\\|","_",i)
  collect_all = collect_all[grepl(pair, collect_all$pair),]
  
  mat.m0=mat.m0[grepl(i,mat.m0$stimulus),]
  
  mat.m0$gene.time = paste0(mat.m0$gene,"_", mat.m0$time,"hr")
  mat.m0.agg=aggregate(mat.m0[grepl("",mat.m0$stimulus),], by = list(mat.m0$gene.time), FUN = function(x){(mean(x))}) #
  mat.m0.agg$cc = collect_all$cc[match(mat.m0.agg$Group.1, paste0(collect_all$gene, "_",collect_all$time))]
  mat.m0$dev2mean = mat.m0.agg$mean[match(paste0(mat.m0$gene,"_",mat.m0$time,"hr"), mat.m0.agg$Group.1)]
  mat.m0$dev2mean = abs(mat.m0$mean-mat.m0$dev2mean)^2
  
  pseudobulk = pseudobulk[,grepl(paste0(timept,"hr"), colnames(pseudobulk))];
  pseudobulk = pseudobulk[,grepl(i, colnames(pseudobulk))]
  pseudobulk$diff = (pseudobulk[,1]-pseudobulk[,2])
  mat.m0$diff = pseudobulk$diff[match(rownames(pseudobulk), mat.m0$gene)]
    
  # wgcna = read.delim("./wgcna/WGCNA_3modules_20210410.txt")
  # mat.m0$gceneset = as.factor(wgcna$net.colors[match(mat.m0$gene, wgcna$gene)])
  mat.m0.agg=aggregate(mat.m0[grepl("",mat.m0$time),], by = list(paste0(mat.m0$gene,"_", mat.m0$time,"hr")), FUN = function(x){(mean(x))}) #
  mat.m0.agg$cc = collect_all$cc[match(mat.m0.agg$Group.1, paste0(collect_all$gene, "_",collect_all$time))]
  # mat.m0.agg$geneset = (mat.m0$geneset[1:nrow(mat.m0.agg)])
  # cor.mean = cor.test(mat.m0.agg$cc[!grepl("^0.25$|1|8",mat.m0.agg$time)],
  #                     mat.m0.agg$dev2mean[!grepl("^0.25$|1|8",mat.m0.agg$time)], method = "pearson")$estimate
  # cor.ff = cor.test(mat.m0.agg$cc[!grepl("^0.25$|1|8",mat.m0.agg$time)],
  #                     mat.m0.agg$fano[!grepl("^0.25$|1|8",mat.m0.agg$time)], method = "pearson")$estimate
  # cor.cv = cor.test(mat.m0.agg$cc[!grepl("^0.25$|1|8",mat.m0.agg$time)],
  #                   (sqrt(mat.m0.agg$var)/mat.m0.agg$mean)[!grepl("^0.25$|1|8",mat.m0.agg$time)], method = "pearson")$estimate
  
  #########USE THIS
  
  mat.m0.agg.3hr= mat.m0.agg[grepl(timept,mat.m0.agg$time),]
  rownames(mat.m0.agg.3hr)= gsub(paste0("_",timept,"hr"), "",mat.m0.agg.3hr$Group.1)
  # p=p+
  # ggplot(mat.m0.agg.3hr, aes(log2(var),log2(cc)))+ #aes(log2(dev2mean),log2(fano)))+
  #   geom_point(aes(color = cc), size =3)+theme_bw(base_size = 14)+
  #   scale_color_gradient(low="blue", high="red")+ggtitle(i)+
  #   # facet_wrap(~time, ncol=4)+
  #   geom_text_repel(      data = mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr)) , ] ,
  #                         aes(label = c(rownames(mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr) ),] ) )),
  #                         
  #     # data = mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr)) |(log2(mat.m0.agg.3hr$fano)<(-2)) , ] ,
  #                   # aes(label = c(rownames(mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr) ),] ), rownames(mat.m0.agg.3hr[log2(mat.m0.agg.3hr$fano)<(-2),]) )),
  #                   box.padding = 1.5) +
  #   # geom_smooth(method = "lm", se=F)+stat_cor(label.y = -4,method = "pearson", p.digits = 0)+
  #   ylab("log2(mean squared deviation)")+xlab("avg. Fano factor")
  # clusters = readxl::read_excel("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/1-s2.0-S2405471217300121-mmc4.xlsx")
  # clusters = readxl::read_excel("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/manual_1-s2.0-S2405471217300121-mmc4.xlsx")
  # mat.m0.agg.3hr$clusters = clusters$`Cluster ID`[match(rownames(mat.m0.agg.3hr), clusters$Symbol)]
  
  clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
  mat.m0.agg.3hr$clusters = clusters$manual_fix[match(rownames(mat.m0.agg.3hr), clusters$gene)]
  mat.m0.agg.3hr$clusters = gsub("mB","mC", mat.m0.agg.3hr$clusters)
  mat.m0.agg.3hr$clusters = ifelse(mat.m0.agg.3hr$clusters=="mA","AP1",
                                   ifelse(mat.m0.agg.3hr$clusters=="mC","NFkB",
                                          ifelse(mat.m0.agg.3hr$clusters=="mD","NFkB&p38",
                                                 ifelse(mat.m0.agg.3hr$clusters=="mE","NFkB|IRF","IRF")))) 
  mat.m0.agg.3hr$clusters = factor(mat.m0.agg.3hr$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))                                  
  
  if(0){
  mat.m0.agg.3hr$label = ifelse(grepl("Nfkbiz|Zc3h12a|Ifi205|Trim21|Cmpk2|Ccl5|Cxcl10|Socs3|Saa3|Tnf$|Il1a|Il1b|Il6|Ifit3$|Tnfaip3", rownames(mat.m0.agg.3hr)), rownames(mat.m0.agg.3hr),"")
  mat.m0.agg.3hr$label = ifelse(mat.m0.agg.3hr$cc>0.4, rownames(mat.m0.agg.3hr),
                                ifelse(mat.m0.agg.3hr$dev2mean>0.75, rownames(mat.m0.agg.3hr),""))
  p1=ggplot(mat.m0.agg.3hr, aes((fano),(cc)))+ #aes(log2(dev2mean),log2(fano)))+
    geom_point(aes( fill=dev2mean), size =4,shape=21)+
    scale_fill_gradient(low="lightblue", high="red")+ggtitle(i)+
    # facet_wrap(~time, ncol=4)+
    geom_text_repel(      mapping = aes(label = label),                          
                          # data = mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr)) |(log2(mat.m0.agg.3hr$fano)<(-2)) , ] ,
                          # aes(label = c(rownames(mat.m0.agg.3hr[grepl("Tnf$|Il1a|Il1b|Il6|Ifit3$|Mx1|Mx2|Oasl1|Il10|Tnfaip3|Gsr|Bcl2a1a|Bcl2ald|Bcl2l11", rownames(mat.m0.agg.3hr) ),] ), rownames(mat.m0.agg.3hr[log2(mat.m0.agg.3hr$fano)<(-2),]) )),
                          box.padding = .5 ,size = 5 ) +
    geom_smooth(method = "lm", se=F)+stat_cor(label.y = 1.5,label.x =1, method = "pearson", p.digits = 0)+theme_bw(base_size = 14)+
    ylab("channel capacity")+xlab("avg. Fano factor")+ylim(0,1.0)
  p2=ggplot(mat.m0.agg.3hr, aes((dev2mean),(cc)))+ #aes(log2(dev2mean),log2(fano)))+
    geom_point(aes( fill=fano), size =4,shape=21)+
    scale_fill_gradient(low="lightblue", high="red")+ggtitle(i)+
    # facet_wrap(~time, ncol=4)+
    geom_text_repel(      mapping = aes(label = label),
      
      # data = mat.m0.agg.3hr[grepl("Ifi205|Trim21|Cmpk2|Ccl5|Cxcl10|Socs3|Saa3|Tnf$|Il1a|Il1b|Il6|Ifit3$|Tnfaip3", rownames(mat.m0.agg.3hr)) , ] ,
      # mapping=aes(label = c(rownames(mat.m0.agg.3hr[grepl("Ifi205|Trim21|Cmpk2|Ccl5|Cxcl10|Socs3|Saa3|Tnf$|Il1a|Il1b|Il6|Ifit3$|Tnfaip3", rownames(mat.m0.agg.3hr) ),] ) )),
                    box.padding = .5 ,size = 5 ) +
    geom_smooth(method = "lm", se=F)+stat_cor(label.y = 1.5,label.x =0, method = "pearson", p.digits = 0)+theme_bw(base_size = 14)+
    ylab("channel capacity")+xlab("mean squared deviation (MSD)")+ylim(0,1.0)
  print(p1/p2)
  }
  
  mat.m0.agg.3hr$label = ifelse(abs(mat.m0.agg.3hr$cc)>0.5, rownames(mat.m0.agg.3hr),"")
  
  pos = position_jitter(width = 0.2, seed = 1)
  p3=ggplot(mat.m0.agg.3hr[grepl("NFkB",mat.m0.agg.3hr$clusters),], aes(clusters, cc))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(aes(fill = diff), size = 4,shape=21,position = pos)+
    scale_fill_gradient(low="dodgerblue2",  high="red",)+
    # scale_fill_gradient2_tableau(palette = "Orange-Blue-White Diverging", limits = c(-6, 6))+
    geom_text_repel(mapping = aes(label = label),size = 5, position = pos )+
    theme_bw(base_size = 18)+
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+xlab("")+
    ggtitle(i)+labs(fill = "Diff")+ylab("max MI (bits)")+ylim(0,1)
  print(p3)
  #########
  # p = p+p3

}
p


# ggplot(mat.m0.agg[!grepl("^0.25$|1|8",mat.m0.agg$time),], aes(dev2mean,cc))+geom_point(aes(color = as.factor(time)))+theme_bw(base_size = 14)+
#   facet_wrap(~time, ncol=4)+geom_smooth(method = "lm", se=F)+stat_cor(method = "pearson", label.y= 1.7, p.digits = 0)
# ggplot(mat.m0.agg[!grepl("^0.25$|1|8",mat.m0.agg$time),], aes(mean/sqrt(var),cc))+geom_point(aes(color = as.factor(time)))+theme_bw(base_size = 14)+
#   facet_wrap(~time, ncol=4)+geom_smooth(method = "lm", se=F)+stat_cor(method = "pearson", p.digits = 0)+xlim(0,5)
# ggplot(mat.m0.agg[!grepl("^0.25$|1|8",mat.m0.agg$time),], aes(fano,cc))+geom_point(aes(color = as.factor(geneset)))+theme_bw(base_size = 14)+
#   facet_wrap(~time, ncol=4)+geom_smooth(method = "lm", se=F)+stat_cor(method = "pearson", p.digits = 0)+ylim(0,1.5)


ggplot(mat.m0.agg[!grepl("^0.25$|1|8",mat.m0.agg$time),], aes(fano,skew))+geom_point(aes(color = as.factor(time)))+theme_bw(base_size = 14)+
  facet_wrap(~time, ncol=4)+geom_smooth(method = "lm", se=F)+stat_cor(method = "pearson", p.digits = 0)+ylim(0,1.5)
fit <- lm(cc ~ mean+fano, 
          data = mat.m0.agg[!grepl("0$",mat.m0.agg$time),]) 
summary(fit)

table_meanff= data.frame(cor.mean =cor.mean, cor.ff=cor.ff, cor.cv, stim_pair = "all")
stim_list = c("CpG", "IFNb", "LPS", "P3CSK", "PIC", "TNF")
for (i in 1:length(stim_list)){
  for (j in (i+1):length(stim_list) ){
    print(paste(stim_list[i], stim_list[j]))
    mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T)
    mat.m0=mat.m0[grepl(paste0(stim_list[i],"|",stim_list[j]),mat.m0$stimulus),]
    mat.m0$gene.time = paste0(mat.m0$gene,"_", mat.m0$time,"hr")
    mat.m0.agg=aggregate(mat.m0[grepl("",mat.m0$stimulus),], by = list(mat.m0$gene.time), FUN = function(x){(mean(x))}) #
    mat.m0.agg$cc = collect_all$cc[match(mat.m0.agg$Group.1, paste0(collect_all$gene, "_",collect_all$time))]
    mat.m0$dev2mean = mat.m0.agg$mean[match(paste0(mat.m0$gene,"_",mat.m0$time,"hr"), mat.m0.agg$Group.1)]
    mat.m0$dev2mean = abs(mat.m0$mean-mat.m0$dev2mean)^2
    mat.m0.agg=aggregate(mat.m0[grepl("",mat.m0$time),], by = list(paste0(mat.m0$gene,"_", mat.m0$time,"hr")), FUN = function(x){(mean(x))}) #
    mat.m0.agg$cc = collect_all$cc[match(mat.m0.agg$Group.1, paste0(collect_all$gene, "_",collect_all$time))]
    cor.mean = cor.test(mat.m0.agg$cc[!grepl("^0.25$",mat.m0.agg$time)],
                        mat.m0.agg$dev2mean[!grepl("^0.25$",mat.m0.agg$time)], method = "pearson")$estimate
    cor.ff = cor.test(mat.m0.agg$cc[!grepl("^0.25$",mat.m0.agg$time)],
                      mat.m0.agg$fano[!grepl("^0.25$",mat.m0.agg$time)], method = "pearson")$estimate
    cor.cv = cor.test(mat.m0.agg$cc[!grepl("^0.25$",mat.m0.agg$time)],
                      (sqrt(mat.m0.agg$var)/mat.m0.agg$mean)[!grepl("^0.25$",mat.m0.agg$time)], method = "pearson")$estimate
    tmp = data.frame(cor.mean, cor.ff, cor.cv, stim_pair = paste0(stim_list[i],"_",stim_list[j]))
    table_meanff = rbind(table_meanff, tmp)
  }
}
library(ggrepel)
ggplot(na.omit(table_meanff), aes(abs(cor.cv), abs(cor.mean), label =stim_pair))+
  geom_point(size = 3, shape = 21, fill = "blue")+geom_text_repel()+theme_bw(base_size = 14)

# graph PMs-----
ggplot(dcast,aes(PM_B6.LFD-PM_B6.old, PM_B6.LFD-PM_B6.HFD))+geom_point(aes(color = geneset), size = 2)+
  # geom_text_repel(subset(dcast, PM_B6.old > 0.4),mapping = aes(label = gene), size = 4)+
  geom_abline(slope = 1, intercept = 0)+theme_bw(base_size = 18)

# graph WT vs MM----
dcast$cc.diff = dcast$BMDM2_WT - dcast$BMDM2_MM
ggplot(dcast,aes(BMDM2_WT, BMDM2_MM))+geom_point()+geom_abline(slope = 1, intercept = 0)+theme_bw()
#compare to 10x
MI_WT.10x = read.delim("./infotheo/SLEMI_singlegene_collectall_allBMDM10x_WT_stimulated.txt")
MI_MM.10x = read.delim("./infotheo/SLEMI_singlegene_collectall_allBMDM10x_MM_stimulated.txt")
dcast$WT_10x = MI_WT.10x$cc[match(dcast$gene, MI_WT.10x$gene)]
dcast$MM_10x = MI_MM.10x$cc[match(dcast$gene, MI_MM.10x$gene)]
dcast$cc.diff_10x =dcast$WT_10x-dcast$MM_10x
ggplot(dcast,aes(WT_10x, MM_10x))+geom_point()+geom_abline(slope = 1, intercept = 0)+theme_bw()
ggplot(dcast, aes(BMDM2_WT, WT_10x))+geom_point()+geom_abline(slope = 1, intercept = 0)
ggplot(dcast, aes(BMDM2_MM, MM_10x))+geom_point()+geom_abline(slope = 1, intercept = 0)
ggplot(dcast, aes((cc.diff), (cc.diff_10x)))+geom_point()+geom_smooth(method = "lm")
  
# other plotting----
collect_all$time.val = as.numeric(gsub("hr", "", collect_all$time))
ggplot(collect_all[!grepl("0.5hr|^5hr", collect_all$time)&grepl("Tnf$|Cxcl10|Ccl5|Il6", collect_all$gene),], aes(time, cc))+
  geom_line(aes(group = gene, color = gene), size = 1)+
  # geom_smooth(aes(group = gene, color = gene), method = "loess", se=F)+
  ylab("channel capacity")+
  geom_point()+theme_bw(base_size = 14)+ylim(c(0,1.2))

#plot M0, M1, M2 together----
collect_all.M0 = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all.M1 = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt")
collect_all.M2 = read.delim("./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt")

collect_all.M0$type = "M0"
collect_all.M1$type = "M1"
collect_all.M2$type = "M2"

collect_all = rbind(collect_all.M0, collect_all.M1, collect_all.M2)
collect_all$time.val = as.numeric(gsub("hr", "", collect_all$time))
# write.table(collect_all, "./infotheo/SLEMI_singlegene_M0M1M2_500genes.txt",sep = "\t", quote = F, row.names = F)
# write.table(collect_all, "./infotheo/SLEMI_singlegene_M0M1M2_ISnorm_500genes.txt",sep = "\t", quote = F, row.names = F)

ggplot(collect_all[!grepl("24hr", collect_all$time)&grepl("Cxcl10$|Ccl5|Tnf$|Il6$", collect_all$gene),], aes(time.val, cc))+
  # geom_line(aes(group = type, color = type), size = 1)+
  geom_smooth(aes(group = type, color = type), method = "loess", se=F)+
  ylab("max MI (bits)")+facet_grid(~gene)+
  geom_point()+theme_calc(base_size = 14)+xlab("Time (hours)")
ggplot(collect_all[grepl("3hr", collect_all$time)&grepl("Cxcl10$", collect_all$gene),], aes(type, cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = type), size = 1)+
  # geom_smooth(aes(group = gene, color = gene), method = "loess", se=F)+
  ylab("channel capacity")+theme_bw(base_size = 14)+xlab("type")


collect_all = read.delim("./infotheo/SLEMI_singlegene_M0M1M2_ISnorm_500genes.txt")
dcast0 = dcast(collect_all[!grepl("0.25|0.5hr|^5hr", collect_all$time),], gene~type+time, value.var = "cc")
dcast0$M1_1hr.diff = dcast0$M0_1hr-dcast0$M1_1hr
dcast0$M1_3hr.diff = dcast0$M0_3hr-dcast0$M1_3hr
dcast0$M1_8hr.diff = dcast0$M0_8hr-dcast0$M1_8hr
dcast0$M2_1hr.diff = dcast0$M0_1hr-dcast0$M2_1hr
dcast0$M2_3hr.diff = dcast0$M0_3hr-dcast0$M2_3hr
dcast0$M2_8hr.diff = dcast0$M0_8hr-dcast0$M2_8hr
dcast0$diff.sum_1hr = dcast0$M1_1hr.diff+dcast0$M2_1hr.diff
dcast0$diff.sum_3hr = dcast0$M1_3hr.diff+dcast0$M2_3hr.diff
dcast0$diff.sum_8hr = dcast0$M1_8hr.diff+dcast0$M2_8hr.diff
wgcna = read.delim("./wgcna/WGCNA_5modules_20210121.txt")
wgcna = read.delim("./wgcna/WGCNA_3modules_20210410.txt")
clusters = readxl::read_excel("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/1-s2.0-S2405471217300121-mmc4.xlsx")
motif.categorize = read.delim("F://scRNAseq_macro/bulk_rnaseq/motif_count_genecategorization.txt")
clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
dcast0$geneset = as.factor(wgcna$net.colors[match(dcast0$gene, wgcna$gene)])
dcast0$geneset = as.factor(motif.categorize$type[match(dcast0$gene, wgcna$gene)])
dcast0$geneset = dcast0$clusters = clusters$`Cluster ID`[match(dcast0$gene, clusters$Symbol)]
dcast0$geneset = dcast0$clusters = clusters$manual_fix[match(dcast0$gene, clusters$gene)]
dcast0$geneset = gsub("mB","mC", dcast0$geneset)
dcast0$geneset = ifelse(dcast0$geneset=="mA","AP1",
                                 ifelse(dcast0$geneset=="mC","NFkB",
                                        ifelse(dcast0$geneset=="mD","NFkB&p38",
                                               ifelse(dcast0$geneset=="mE","NFkB|IRF","IRF")))) 
dcast0$geneset = factor(dcast0$geneset, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))                                  



ggplot(na.omit(dcast0), aes(geneset, M0_1hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+
  ggplot(na.omit(dcast0), aes(geneset, M0_3hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+
  ggplot(na.omit(dcast0), aes(geneset, M0_8hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)
  
ggplot(na.omit(dcast0), aes(geneset, M0_3hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+ylim(0,1.55)+
  ggplot(na.omit(dcast0), aes(geneset, M1_3hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+ylim(0,1.55)+
  ggplot(na.omit(dcast0), aes(geneset, M2_3hr))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+ylim(0,1.55)

ggplot(dcast0, aes(M1_1hr.diff, M2_1hr.diff))+ geom_point(size = 1)+theme_minimal(base_size = 14)+geom_abline(slope = 1, intercept = 0)
  
ggplot(na.omit(dcast0[order(dcast0$geneset),]), aes(M1_3hr.diff, M2_3hr.diff))+ geom_point(aes(color = geneset),alpha = 0.8,size = 2)+
  xlim(c(-0.5,0.75))+ylim(c(-0.5,.75))+
  theme_bw(base_size = 16)+geom_abline(slope = 1, intercept = 0, size=1,linetype="dashed")+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
ggplot(na.omit(dcast0[grepl("",dcast0$geneset),]), aes(geneset, M1_3hr.diff))+ theme_bw(base_size = 14)+theme(legend.position = "None")+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)+ #ylim(-.7,1.1)+
ggplot(na.omit(dcast0[grepl("",dcast0$geneset),]), aes(geneset, M2_3hr.diff))+ theme_bw(base_size = 14)+#ylim(-.7,1.1)+
  geom_boxplot(outlier.shape = NA)+geom_point(aes(color = geneset),position = "jitter", size = 2)+geom_hline(yintercept = 0)


genes.to.label = c("Cxcl10","Tgtp1","Rsad2","Irf7","Trim21","Ifi205","Gbp7",
                   "Ccl5", "Tnf", "Tnfaip3","Gclm", "Il6",
                   "Tnfsf9","Bcl2a1a","Bcl2a1d","Socs3")
rownames(dcast0) = dcast0$gene
p1=ggplot(dcast0,  aes(M1_3hr.diff, M2_3hr.diff))+geom_point(aes(color = M0_3hr), size =2)+
  geom_point(data = subset(dcast0, subset = gene %in% genes.to.label),size = 4, aes(color = M0_3hr)) + 
  geom_point(data = subset(dcast0, subset = gene %in% genes.to.label),size = 4, color = "red",shape = 21) + 
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+theme_bw(base_size = 14)
p1 <- LabelPoints(plot = p1, points = genes.to.label, color = "red", size = I(5),repel = T, xnudge=0,ynudge=0)
p1

ggplot(dcast0, aes(M1_3hr.diff, M2_3hr.diff))+ geom_point(size = 1)+theme_minimal(base_size = 14)+geom_abline(slope = 1, intercept = 0)


ggplot(dcast0, aes(rank(M1_3hr.diff), M1_3hr.diff))+ geom_bar(stat = "identity", aes(fill = geneset))+
  scale_fill_manual(values=c("gray","#F8766D","#00BA38", "#619CFF"))+geom_hline(yintercept = 0)+
  theme_classic(base_size = 14)
ggplot(dcast0, aes(rank(M2_3hr.diff), M2_3hr.diff))+ geom_bar(stat = "identity", aes(fill = geneset))+
  scale_fill_manual(values=c("gray","#F8766D","#00BA38", "#619CFF"))+geom_hline(yintercept = 0)+
  theme_classic(base_size = 14)
# write.table(dcast0, "./infotheo/SLEMI_singlegene_M0M1M2_500genes_dcast.txt",sep = "\t", quote = F, row.names = F)
ggplot(dcast0, aes(rank(M2_8hr.diff), geneset)) + geom_point(aes(color = geneset), position = "jitter") + theme_classic(base_size = 11) + 
  scale_color_manual(values=c("gray","#F8766D","#00BA38", "#619CFF"))+
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1.2), legend.position = "none") 

dcast = dcast(collect_all[!grepl("^5hr|24hr", collect_all$time),], gene~time+type, value.var = "cc")
rownames(dcast) = dcast$gene
colnames(dcast)[-1] = paste0("X", colnames(dcast)[-1])
ggplot(dcast, aes(X1hr_M0, X1hr_M1))+ geom_point(size = 1)+theme_bw(base_size = 14)+geom_abline(slope = 1, intercept = 0)
ggplot(dcast, aes(X8hr_M0, X8hr_M1))+ geom_point(size = 1)+theme_bw(base_size = 14)+geom_abline(slope = 1, intercept = 0)
ggplot(dcast, aes(X3hr_M1, X3hr_M2))+ geom_point(size = 1)+theme_bw(base_size = 14)+geom_abline(slope = 1, intercept = 0)

pheatmap(na.omit(dcast[rowSums(dcast[,-1])>1.5,-1]), scale = "row", 
         cluster_cols = F, cluster_rows = T,clustering_method = "ward.D2", show_rownames = F,
         # gaps_col = c(4,9),
         gaps_col = c(3,6,9)
         )
pheatmap(na.omit(dcast[rowSums(dcast[,-1])>4.0,-1]), scale = "none", 
         cluster_cols = F, cluster_rows = T,clustering_method = "ward.D2", show_rownames = T,
         # gaps_col = c(4,9),
         gaps_col = c(3,6,9))
#keep if any col has cc>0.75
pheatmap(na.omit(dcast[,-1][rowSums(dcast[,-1] > .75) > 0, ]), scale = "none", 
         cluster_cols = F, cluster_rows = T,clustering_method = "average", show_rownames = T,
         colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
         # gaps_col = c(4,9),
         gaps_col = c(3,6,9))

mutual = readRDS("./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
genes = unlist(mutual[41,4])
pheatmap(na.omit(dcast[genes[order(genes)],-1][,grepl("_M0",colnames(dcast))]), scale = "none", clustering_method = "ward.D2", 
         # gaps_col = c(3,6,9),
         colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
         cluster_cols = F, cluster_rows = F)


do_kmeans_clustering(na.omit(dcast[,-1]), k_clusters = 5, cluster_cols = F,
                     show_colnames = T, show_rownames = F, colseps = c(3,6,9))

#write out clusters of channel capacity-------------------
require(pheatmap)
set.seed(1)
scaled = t(scale(t(na.omit(dcast[,-1]))))
k=5
kmcluster <- kmeans(scaled, iter.max = 1000, centers = k, 
                    nstart = 1)
mat <- cbind(na.omit(dcast[,-1]), cluster = kmcluster$cluster)
mat <- mat[order(mat$cluster), ]
count <- 0
for (i in 1:k) {
  count[i] <- length(kmcluster$cluster[kmcluster$cluster == 
                                         i])
}
rowseps <- cumsum(count)
p = pheatmap(mat[, -ncol(mat)], cluster_rows = F, 
             cluster_cols = F, scale = "row", clustering_method = "ward.D2", 
             gaps_col = c(3,6,9), gaps_row = rowseps, 
             colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu"))[2:11])(103), 
             breaks = c(-4, seq(-1.5,1.5, length = 100), 4), border_color = NA)

mat = mat[, -ncol(mat)]
count <- 1;
for(j in 1:k){
  print(j)
  count[j+1] <- length(kmcluster$cluster[kmcluster$cluster == j])
  clust = mat[cumsum(count)[j]:cumsum(count)[j+1]-1, ]
  write.table(data.frame(gene = rownames(clust)), paste0("./infotheo/SLEMI_singlegene_clust",j,"_M0M1M2.txt"), sep = '\t', quote = F, row.names = F)
}

# findMotifs.pl SLEMI_singlegene_clust1_M0M1M2.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/infotheo/MotifResults1/ -start -1000 -end 100 -len 8,10 -p 4 &
# findMotifs.pl SLEMI_singlegene_clust2_M0M1M2.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/infotheo/MotifResults2/ -start -1000 -end 100 -len 8,10 -p 4 &
# findMotifs.pl SLEMI_singlegene_clust3_M0M1M2.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/infotheo/MotifResults3/ -start -1000 -end 100 -len 8,10 -p 4 &
# findMotifs.pl SLEMI_singlegene_clust4_M0M1M2.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/infotheo/MotifResults4/ -start -1000 -end 100 -len 8,10 -p 4 &
# findMotifs.pl SLEMI_singlegene_clust5_M0M1M2.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/infotheo/MotifResults5/ -start -1000 -end 100 -len 8,10 -p 4 &


ggplot(collect_all.M1[!grepl("24hr", collect_all.M1$time),], aes(time, cc))+geom_violin()+
  geom_point(aes(color = time),position = "jitter")+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+ylab("channel capacity")+
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_text_repel(subset(collect_all.M1[!grepl("24hr", collect_all.M1$time),], cc > 0.9), 
                  mapping = aes(label = gene), size = 4)+ylim(0,1.5)

ggplot(collect_all.M2[!grepl("24hr", collect_all.M2$time),], aes(time, cc))+geom_violin()+
  geom_point(aes(color = time),position = "jitter")+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+ylab("channel capacity")+
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_text_repel(subset(collect_all.M2[!grepl("24hr", collect_all.M2$time),], cc > 0.9), 
                  mapping = aes(label = gene), size = 4)+ylim(0,1.5)

#########################################################################################
# collect by forward selection on individual gene list for each timepoint-----
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0all_gt80_500genes.txt")

#testing normalization comparison----
collect_all.0 = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_SCT_500genes.txt")
collect_all.0 = collect_all.0[order(collect_all.0$cc, decreasing = T), ]
genesets[[11]] = unique(unlist(c(collect_dimensionbest$gene[grepl("9", collect_dimensionbest$dim) ])))
collect_all = data.frame(gene = genesets[[11]])
collect_all$cc = collect_all.0$cc[match(collect_all$gene, collect_all.0$gene)]

collect = data.frame()
# for (i in c("0.25hr","1hr","3hr","8hr")){
for (i in c("3hr")){
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i, "_DBEC.rds"))
  macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  macro = subset(macro, subset= timept ==i)
  
  for (m in seq(1:23)){
    print(m)
    # genes = c(collect_all$gene[grepl(i, collect_all$time)][(1:m)])
    genes = collect_all$gene[1:m]
    print(genes)
    
    if(0){ #include unstim
      # macro = readRDS(paste0("./output/macrophage_tutorial3_subsample_8hr_sctransform.rds"))
      macro.0hr = readRDS(paste0("./output/macrophage_tutorial3_subsample_0hr_sctransform.rds"))
      macro = merge(macro.0hr, macro)
    }
    
    data = macro[["ISnorm"]]@data
    data = data.frame(data)
    meta = macro@meta.data
    colnames(data) = meta$stimulus
    data = data[rowSums((data==0))<ncol(data),] #remove genes that are all 0
    # colnames(data) = paste0(meta$stimulus, "_",meta$timept)
    # rownames(data) = gsub("-ENSMUST..*","", macro[["SCT"]]@counts@Dimnames[[1]])
    # rownames(data) = gsub("..NM..*","", rownames(data))
    my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
    my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
    rownames(my.dataframe) = seq(1:nrow(my.dataframe))
    
    my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
    str(my.dataframe)
    
    
    output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                            # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i,"_",m, "_rmUnstim"), 
                                            testing = F, boot_num = 50,testing_cores = 4)
    
    tmp = data.frame(type = i, genes = m, cc = output_capacity$cc)
    collect = rbind(collect, tmp)
    
  }
}
# write.table(collect, "./infotheo/SLEMI_M0all_gt80_forwardselection_singlegeneMI_top100.txt", quote = F, sep = "\t",row.names = F)
collect = read.delim( "./infotheo/SLEMI_M0all_gt80_forwardselection_singlegeneMI_top100.txt")
for(i in c(1:600)){
    collect$difference[i] = collect$cc[i+1]-collect$cc[i]
}
collect.tmp = collect[!grepl("0.5hr|^5hr", collect$type) &!grepl("100", collect$genes), ]
ggplot(collect[!grepl("0.5hr|^5hr", collect$type) &!grepl("100", collect$genes), ], 
       aes(fill=type, y=difference, x= genes )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("cc difference")+ 
  geom_line(aes(group = type, color = type), size = 1)+theme_bw(base_size = 20)

ggplot(collect[grepl("3hr", collect$type), ], aes(fill=type, y=cc, x= genes )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = type, color = type), size = 2)+ylim(0, 2.7)+theme_bw(base_size = 20)

ggplot(collect[!grepl("0.5hr|^5hr", collect$type), ], aes(fill=type, y=cc, x= genes )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = type, color = type), size = 2)+ylim(0, 2.7)+theme_bw(base_size = 20)

ggplot(collect[(!grepl("0.5hr|^5hr", collect$type)& grepl("14|13|^5$|^10$|20$", collect$genes)), ], 
       aes(fill=type, y=cc, x= type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = as.factor(genes), color = as.factor(genes)), size = 2)+ylim(0, 2.7)+theme_bw(base_size = 20)
#########################################################################################
# collect by keep top 20 forward selection on individual gene list for each timepoint-----
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt")
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt")
collect_all = collect_all[order(collect_all$cc, decreasing = T), ]
collect_all_best = aggregate(collect_all$cc, by = list(collect_all$time), max)
collect_all_best$dim = 0
collect_all_best$gene = c("Socs3", "Tnf", "Cmpk2", "Irf7") #M0
collect_all_best$gene = c("Tnf", "Tnf", "Ifit3","Tnf", "Ccl5", "Ccl5") #M1
collect_all_best$gene = c("Tnf", "Tnf", "Cmpk","Ifit3", "Ifit3") #M2
colnames(collect_all_best) = c("time", "cc", "dim", "gene")
collect_all_best = collect_all_best[, c("time", "dim","cc",  "gene")]

collect = data.frame()
collect_dimensionbest = data.frame()

for (i in c("0.25hr","1hr","3hr","8hr")){
  print(i)
  macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  macro = subset(macro, subset= timept ==i)
  
  collect_dimension = (collect_all[grepl(i, collect_all$time),][(1:20),]) #start 1D
    
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  
  for (d in seq(1:5)){
    print(paste0("dimension: ",d))
    
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    genesets = collect_dimension$gene[c(1:20)]
    print(genesets)
    
    collect_dimension = data.frame() #start over once got the top20
    for (g in 1:length(genesets)){
      genes = genesets[[g]]
      print(genes)
      
      other_genes = c(colnames(my.dataframe)[!colnames(my.dataframe) %in% genes])
      other_genes
      for (a in 1:length(other_genes)){
        added_gene = other_genes[a+1]
        added_gene
        my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("label", genes, added_gene)]
        str(my.dataframe.subset)
        
        output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", response = colnames(my.dataframe.subset)[-1],
                                                output_path = NULL, testing = F)
        
        tmp = data.frame(time = i,  dim = d, cc = output_capacity$cc)
        tmp$gene = list(c(genes, added_gene))
        collect_dimension = rbind(collect_dimension, tmp)
      }
      
    }
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    collect_dimensionbest = rbind(collect_dimensionbest,
                                  data.frame(collect_dimension[1,]))
    
  }
}

collect_dimensionbest = readRDS("./infotheo/collect_dimensionbest_ISnorm2.rds")
collect_dimensionbest = readRDS("./infotheo/collect_dimensionbest_M1_ISnorm2.rds")
collect_dimensionbest = readRDS("./infotheo/collect_dimensionbest_M2_ISnorm2.rds")
collect_dimensionbest = rbind(collect_dimensionbest, collect_all_best[-c(1,4,6),])
collect_dimensionbest = collect_dimensionbest[order(collect_dimensionbest$time, collect_dimensionbest$dim),]
collect_dimensionbest$dim = collect_dimensionbest$dim+1
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_M1_ISnorm_Feb2021.rds")
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_M2_ISnorm_Feb2021.rds")
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_ISnorm2_Apr2021.rds")
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_M1_ISnorm2_Apr2021.rds")
# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_M2_ISnorm2_Apr2021.rds")
collect_dimensionbest = readRDS("./infotheo/collect_dimensionbest_ISnorm2_Apr2021.rds")

ggplot(collect_dimensionbest, aes(dim, cc))+ylim(1, 2.7)+xlim(0, 20)+
  geom_point(aes(color = time, group = time), size =3)+
  geom_line(aes(color = time, group = time),size =1)+theme_bw(base_size = 20)+
  geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  xlab("dimension (#genes)")+ylab("channel capacity")
ggplot(collect_dimensionbest[grepl("3hr", collect_dimensionbest$time),], aes(dim, cc))+ylim(1, 2.7)+xlim(0, 20)+
  geom_point(aes(color = time, group = time), size =3)+
  geom_line(aes(color = time, group = time),size =1)+theme_bw(base_size = 14)+
  geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_vline(xintercept = c(2,5,10,15), linetype="dotted", size = 1)+
  xlab("dimension (#genes)")+ylab("channel capacity")
#############################################################################################
# channel capacity for select subsets of genes-----
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
genesets = list()
genesets[[1]] =  colnames(my.dataframe)[grepl("Mmp|Timp1$|Nos2|Tnf$|Il6$", colnames(my.dataframe))] #local inflammation
genesets[[2]] =  colnames(my.dataframe)[grepl("Cxcl1$|Cxcl2$|Cxcl3$|Cxcl5$|Cxcl6$|Cxcl8$", colnames(my.dataframe))] #neutrophil
genesets[[3]] =  colnames(my.dataframe)[grepl("Ccl1$|Ccl19$|Ccl20$|Cxcl9|Cxcl10|Cxcl11|Cxcl13", colnames(my.dataframe))] #lymphocyte
genesets[[4]] =  colnames(my.dataframe)[grepl("Ccl2$|Ccl7$|Ccl8$|Ccl13$|Cxcl12|Il1b", colnames(my.dataframe))] #monocyte
genesets[[5]] =  colnames(my.dataframe)[grepl("Il10$|Il4|Il11|Il13", colnames(my.dataframe))]#anti-inflammatory
genesets[[6]] =  colnames(my.dataframe)[grepl("Sod2$|Bcl2|Bcor$|Ptgs2$|Ptges$", colnames(my.dataframe))] #bcl cell death
genesets[[7]] =  colnames(my.dataframe)[grepl("Slc", colnames(my.dataframe))]#transporter
genesets[[8]] = colnames(my.dataframe)[grepl("Rhoh$|Mx1$|Rnd3$|Mx2$|Rnd1$|Arhgap26$|Rhou$|Tagap$|Iigp1$|Rhof$|Gem$|Tgtp1$", colnames(my.dataframe))] #GTPases
genesets[[9]] = unique(c(genesets[[1]],genesets[[2]],genesets[[3]],genesets[[4]],genesets[[5]],
                         genesets[[6]],genesets[[7]],genesets[[8]]))
genesets[[10]] =  colnames(my.dataframe)[grepl("Alas1$|Cmpk2$|Sod2$|Gsr$", colnames(my.dataframe))] #mitochondrial
genesets[[11]] =  colnames(my.dataframe)[grepl("Mmp13$|Tnf$|Il6$|Il27$", colnames(my.dataframe))] #local inflammation
genesets[[12]] =  colnames(my.dataframe)[grepl("Cxcl1$|Cxcl2$|Cxcl3$", colnames(my.dataframe))] #neutrophil
genesets[[13]] =  colnames(my.dataframe)[grepl("Ccl1$|Cxcl9$|Cxcl10|Cxcl11$", colnames(my.dataframe))] #lymphocyte
genesets[[14]] =  colnames(my.dataframe)[grepl("Ccl2$|Ccl7$|Ccl12$|Ccl8$", colnames(my.dataframe))] #monocyte
genesets[[15]] =  colnames(my.dataframe)[grepl("Il10$|Il4$|Il4ra$|Il13$", colnames(my.dataframe))]#anti-inflammatory
genesets[[16]] =  colnames(my.dataframe)[grepl("Bcl2|Bcor$", colnames(my.dataframe))] #bcl cell death
genesets

# top of each timepoint
collect_all.tmp1 = collect_all[grepl("1hr", collect_all$time),]
collect_all.tmp3 = collect_all[grepl("3hr", collect_all$time),]
collect_all.tmp8 = collect_all[grepl("8hr", collect_all$time),]

collect_all.tmp1$rank = rank(-collect_all.tmp1$cc)
collect_all.tmp3$rank = rank(-collect_all.tmp3$cc)
collect_all.tmp8$rank = rank(-collect_all.tmp8$cc)

collect_all.tmp1$rank3hr = collect_all.tmp3$rank[match(collect_all.tmp1$gene, collect_all.tmp3$gene)]
collect_all.tmp1$rank8hr = collect_all.tmp8$rank[match(collect_all.tmp1$gene, collect_all.tmp8$gene)]
collect_all.tmp1$sumrank = rowSums(collect_all.tmp1[, c(4:6)])/3
collect_all.tmp = collect_all.tmp1[order(collect_all.tmp1$sumrank),]
collect_all.tmp$totalrank = rank(collect_all.tmp$sumrank, na.last = T)

# genesets[[1]] = unique(collect_all.tmp$gene[c(1:5)])
# genesets[[2]] = unique(collect_all.tmp$gene[c(1:20)])
# genesets[[3]] = unique(collect_all.tmp$gene[c(1:15)])
# genesets[[4]] = unique(collect_all.tmp$gene[c(1:20)])
# genesets[[5]] = unique(collect_all.tmp$gene[c(1:30)])

#top over all timepoints
genesets[[1]] = unique(collect_all$gene)[c(1:10)]
genesets[[2]] = unique(collect_all$gene)[c(1:20)]
genesets[[3]] = unique(collect_all$gene)[c(1:30)]
genesets[[4]] = unique(collect_all$gene)[c(1:40)]
genesets[[5]] = unique(collect_all$gene)[c(1:50)]
genesets[[6]] = unique(c(collect_all.tmp1$gene[c(1:30)], 
                         collect_all.tmp3$gene[c(1:25)],
                         collect_all.tmp8$gene[c(1:25)]))

genesets[[7]] = unique(c(collect_all.tmp$gene[c(1:20)]))
genesets[[8]] = unique(c(collect_all.tmp$gene[c(1:30)]))
genesets[[9]] = unique(c(collect_all.tmp$gene[c(1:50)]))
genesets[[10]] = unique(c(intersect_all(genesets[[6]], genesets[[9]])))
# top 30 by better forward net selection----
genesets[[11]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim) ])))
genesets[[12]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("1hr|3hr", collect_dimensionbest$time)])))
genesets[[13]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("3hr|8hr", collect_dimensionbest$time) ])))
genesets[[14]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("1hr|8hr", collect_dimensionbest$time) ])))

#for M1 and M2 diff----
genesets[[15]] = dcast0$gene[order(dcast0$diff.sum_3hr, decreasing = T)][1:15]
genesets[[16]] = dcast0$gene[order(dcast0$diff.sum_8hr, decreasing = T)][1:15]
genesets[[17]] = union(dcast0$gene[order(dcast0$diff.sum_3hr, decreasing = T)][1:15], dcast0$gene[order(dcast0$diff.sum_8hr, decreasing = T)][1:15])
genesets

# collect dimension best M0, M1, M2----
collect_dimensionbest = readRDS("./infotheo/collect_dimensionbest_ISnorm2_Apr2021.rds")
collect_dimensionbest.M1 = readRDS("./infotheo/collect_dimensionbest_M1_ISnorm2_Apr2021.rds")
collect_dimensionbest.M2 = readRDS("./infotheo/collect_dimensionbest_M2_ISnorm2_Apr2021.rds")
collect_dimensionbest_all = rbind(data.frame(collect_dimensionbest, type="M0"),
                                  data.frame(collect_dimensionbest.M1, type="M1"),
                                  data.frame(collect_dimensionbest.M2, type="M2"))
collect_dimensionbest_all$type.time = paste(collect_dimensionbest_all$type,collect_dimensionbest_all$time)
ggplot(collect_dimensionbest_all[grepl("3hr", collect_dimensionbest_all$time),], aes(dim, cc))+
  xlim(0, 20)+ geom_point(aes(color = type, group = type.time), size =3)+
  geom_line(aes(color = type, group = type.time, linetype=time),size =1)+theme_bw(base_size = 20)+
  geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_vline(xintercept = c(2,5,10,15), linetype="dotted", size = 1)+
  xlab("dimension (#genes)")+ylab("channel capacity")
ggplot(collect_dimensionbest_all[grepl("3hr", collect_dimensionbest_all$time),], aes(dim, cc))+
  xlim(5, 15)+ylim(1.97, 2.3)+geom_point(aes(color = type, group = type.time), size =2)+
  geom_line(aes(color = type, group = type.time),size =1)+theme_classic(base_size = 16)+
  xlab("dimension (#genes)")+ylab("channel capacity")#+facet_grid(~time)

#plot dimension best M0, M1, M2 all timepts-----
library(plotly);library(mgcv)
collect_dimensionbest_all$timeval = as.numeric(gsub("hr","",collect_dimensionbest_all$time))
data.M0 = collect_dimensionbest_all[grepl("M0",collect_dimensionbest_all$type)&(collect_dimensionbest_all$dim<16),]
data.M1 = collect_dimensionbest_all[grepl("M1",collect_dimensionbest_all$type)&(collect_dimensionbest_all$dim<16),]
data.M2 = collect_dimensionbest_all[grepl("M2",collect_dimensionbest_all$type)&(collect_dimensionbest_all$dim<16),]

modM0 <- loess(cc ~ (timeval)*(dim), data=data.M0) #gam(cc ~ (timeval) + ti(dim), data=data.M0)
modM1 <- loess(cc ~ (timeval)*(dim), data=data.M1) #gam(cc ~ (timeval) + ti(dim), data=data.M1)
modM2 <- loess(cc ~ (timeval)*(dim), data=data.M2) #gam(cc ~ (timeval) + ti(dim), data=data.M2)

timeval.seq <- seq(0, max(collect_dimensionbest_all$timeval, na.rm=TRUE), length=25)
dim.seq <- seq(0, max(collect_dimensionbest_all$dim, na.rm=TRUE), length=25)
predfunM0 <- function(x,y){    newdat <- data.frame(timeval = x, dim=y);  predict(modM0, newdata=newdat)}
predfunM1 <- function(x,y){    newdat <- data.frame(timeval = x, dim=y);  predict(modM1, newdata=newdat)}
predfunM2 <- function(x,y){    newdat <- data.frame(timeval = x, dim=y);  predict(modM2, newdata=newdat)}

fit0 <- outer(timeval.seq, dim.seq, Vectorize(predfunM0))
fit1 <- outer(timeval.seq, dim.seq, Vectorize(predfunM1))
fit2 <- outer(timeval.seq, dim.seq, Vectorize(predfunM2))

color0 <- rep(0, length(t(fit0)));dim(color0) <- dim(t(fit0))
color1 <- rep(1, length(t(fit0)));dim(color1) <- dim(t(fit1))
color2 <- rep(2, length(t(fit0)));dim(color2) <- dim(t(fit2))

fig <- plot_ly(colors = c('red2','green4', 'blue'))  %>% 
  add_markers(x = ~data.M0$timeval, y=data.M0$dim, z=data.M0$cc, color = I("red2")) %>% 
  add_surface(x = ~timeval.seq, y = ~dim.seq, z = t(fit0), opacity = 0.5, 
              surfacecolor=color0,cauto=F,cmax=2,cmin=0)
fig <- fig %>% 
  add_markers(x = ~data.M1$timeval, y=data.M1$dim, z=data.M1$cc, color = I("green4")) %>% 
  add_surface(x = ~timeval.seq, y = ~dim.seq, z = t(fit1), opacity = 0.5, 
              surfacecolor=color1,cauto=F,cmax=2,cmin=0)
fig <- fig %>% 
  add_markers(x = ~data.M2$timeval, y=data.M2$dim, z=data.M2$cc, color = I("blue")) %>% 
  add_surface(x = ~timeval.seq, y = ~dim.seq, z = t(fit2), opacity = 0.5, 
              surfacecolor=color2,cauto=F,cmax=2,cmin=0)
fig

#overlap of best genes at dim10
library(Vennerable)
dimension = 10; time = 8
venn= Venn(list(M0=c(unlist(data.M0[(data.M0$dim==dimension &data.M0$timeval==time),]$gene)), 
     M2 = c(unlist(data.M2[data.M2$dim==dimension &data.M2$timeval==time,]$gene)),
     M1 = c(unlist(data.M1[data.M1$dim==dimension &data.M1$timeval==time,]$gene)) ))
Vennerable::plot(venn, doWeights = T )

# mtplot = ggplot(collect_dimensionbest_all, aes(x=timeval,y=dim)) + 
#   # geom_point(aes(x=timeval,y=dim,fill=cc)) +
#   geom_raster(aes(x=timeval,y=dim, fill=cc), interpolate = T) +
#   # scale_x_continuous(expand=c(0,0)) +
#   # scale_y_continuous(expand=c(0,0)) +  
#   # scale_fill_viridis_c(option = "A")+
#   facet_wrap(~type)+
#   # scale_color_continuous(limits=c(0,2.5))+
#   scale_fill_gradientn("max MI",colours = terrain.colors(10), limits=c(0,2.5))
#   # scale_fill_brewer("max MI",palette = "Spectral")
# plot_gg(mtplot, preview = TRUE)
# plot_gg(mtplot, width=3.5, multicore = TRUE, windowsize = c(1400,866), sunangle=225,
#         zoom = 0.60, phi = 30, theta = -30)
# render_snapshot(clear = TRUE)

# collect_dimensionbest_all <- apply(collect_dimensionbest_all[grepl("3hr", collect_dimensionbest_all$time),],2,as.character)
# write.table(collect_dimensionbest_all,  "F://scRNAseq_macro/SuppTables/TableS5_maxMI_dimensionbest_M0M1M2.txt",sep = "\t",row.names = F,quote=F)



genesets = list()
genesets[[1]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim) ])))
genesets[[2]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("1hr", collect_dimensionbest$time)])))
genesets[[3]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("3hr", collect_dimensionbest$time)])))
genesets[[4]] = unique(unlist(c(collect_dimensionbest$gene[grepl("14", collect_dimensionbest$dim)& grepl("8hr", collect_dimensionbest$time)])))
genesets[[5]] = unique(unlist(c(collect_dimensionbest.M1$gene[grepl("14", collect_dimensionbest.M1$dim) ])))
genesets[[6]] = unique(unlist(c(collect_dimensionbest.M1$gene[grepl("14", collect_dimensionbest.M1$dim)& grepl("1hr", collect_dimensionbest.M1$time)])))
genesets[[7]] = unique(unlist(c(collect_dimensionbest.M1$gene[grepl("14", collect_dimensionbest.M1$dim)& grepl("3hr", collect_dimensionbest.M1$time)])))
genesets[[8]] = unique(unlist(c(collect_dimensionbest.M1$gene[grepl("14", collect_dimensionbest.M1$dim)& grepl("8hr", collect_dimensionbest.M1$time)])))
genesets[[9]] = unique(unlist(c(collect_dimensionbest.M2$gene[grepl("14", collect_dimensionbest.M2$dim) ])))
genesets[[10]] = unique(unlist(c(collect_dimensionbest.M2$gene[grepl("14", collect_dimensionbest.M2$dim)& grepl("1hr", collect_dimensionbest.M2$time)])))
genesets[[11]] = unique(unlist(c(collect_dimensionbest.M2$gene[grepl("14", collect_dimensionbest.M2$dim)& grepl("3hr", collect_dimensionbest.M2$time)])))
genesets[[12]] = unique(unlist(c(collect_dimensionbest.M2$gene[grepl("14", collect_dimensionbest.M2$dim)& grepl("8hr", collect_dimensionbest.M2$time)])))
genesets

library(readxl)
biol_cat = read_excel("F://scRNAseq_macro/scRNAseq_macro/biological_categories.xlsx") 
genesets = list()
genesets[[1]] =  c(unlist(biol_cat[,1]))
genesets[[2]] =  c(unlist(biol_cat[,2]))
genesets[[3]] =  c(unlist(biol_cat[,3]))
genesets[[4]] =  c(unlist(biol_cat[,4]))
genesets[[5]] =  c(unlist(biol_cat[,5]))
genesets[[6]] =  c(unlist(biol_cat[,6]))
genesets[[7]] =  c(unlist(biol_cat[,7]))
genesets[[8]] =  c(unlist(biol_cat[,8]))
genesets[[9]] =  c(unlist(biol_cat[,9]))
genesets

assign_grs = read_excel("F://scRNAseq_macro/SuppTables/TableS3_genes2GRS.xlsx") 
genesets = list()
genesets[[1]] =  c(assign_grs$gene[assign_grs$clusters=="AP1"])
genesets[[2]] =  c(assign_grs$gene[assign_grs$clusters=="NFkB"])
genesets[[3]] =  c(assign_grs$gene[assign_grs$clusters=="NFkB&p38"])
genesets[[4]] =  c(assign_grs$gene[assign_grs$clusters=="NFkB|IRF"])
genesets[[5]] =  c(assign_grs$gene[assign_grs$clusters=="IRF"])
genesets

#test genesets--------
collect = data.frame()
# for (g in 10:16){
# for (g in 11:14){
# for (g in 1:9){
for (g in 1:5){
  print(g)
for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
# for (i in c("0.25hr","1hr","3hr","8hr")){
# for (i in c("PM_B6.LFD", "PM_B6.HFD", "PM_B6.old")){
    
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds"))
  # macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
  macro = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
  # macro = readRDS("output/macrophage_BMDM2_WT_MM_500genes_DBEC.rds")
  # macro = subset(macro, subset= timept ==i & type == "BMDM2_MM")
  
  macro = subset(macro, subset= timept ==i)
  # macro = subset(macro, subset= stimulus!="IFNb")
  
  # macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds")
  # macro = subset(macro, subset= type==i)
  
  if(0){ #include unstim
    # macro = readRDS(paste0("./output/macrophage_tutorial3_subsample_8hr_sctransform.rds"))
    macro.0hr = readRDS(paste0("./output/macrophage_tutorial3_subsample_0hr_sctransform.rds"))
    macro = merge(macro.0hr, macro)
  }
  
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  # rownames(data) = gsub("-ENSMUST..*","", macro[["SCT"]]@counts@Dimnames[[1]])
  # rownames(data) = gsub("..NM..*","", rownames(data))
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  # str(my.dataframe)
  
  genes = genesets[[g]]
  print(genes)
  my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
  str(my.dataframe.subset)
  
  output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", response = colnames(my.dataframe.subset)[-1],
                                          testing = T, testing_cores = 4, TestingSeed = 1234,boot_prob = 0.5, boot_num = 10,traintest_num = 10,
                                          output_path = NULL)
  sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
  tmp = data.frame(type = i, geneset = g, cc = output_capacity$cc, sd = sd)
  # tmp = data.frame(type = i, geneset = g, cc = output_capacity$cc)
  collect = rbind(collect, tmp)


}
}

collect = read.delim("./infotheo/SLEMI_PMexpts_Feb2021_rmUnstim_biological_subsets2.txt")
collect$geneset.state = paste0(collect$geneset, collect$type)
collect$type = factor(collect$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
ggplot(collect[grepl("10|11|13|15", collect$geneset),], aes(fill=type, y=cc, x= type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ ylim(c(0,1.25))+
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+ylim(0, 2.7)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  # geom_bar(stat="identity", position = "dodge")+
  geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_point(size =2)+theme_bw(base_size = 16)


collect = read.delim("./infotheo/SLEMI_M0_rep2only_biological_subsets3.txt")
collect0 = read.delim("./infotheo/SLEMI_M0_rep2only_biological_subsets3.txt")
collect1 = read.delim("./infotheo/SLEMI_M1_IFNg_biological_subsets3.txt")
collect2 = read.delim("./infotheo/SLEMI_M2_IL4_gt80_biological_subsets3.txt")
collect0$state = "M0"
collect1$state = "M1"
collect2$state = "M2"
collect = rbind(collect1, collect2)
collect$geneset.state = paste0(collect$geneset, collect$state)
collect$timeval = as.numeric(gsub("hr", "",collect$type))

collect0$geneset.state = paste0(collect0$geneset, collect0$state)
collect0$timeval = as.numeric(gsub("hr", "",collect0$type))
collect = rbind(collect0, collect)

#plot M0M1M2 biocategory3----
ggplot(collect[grepl("1|3|9", collect$geneset),], aes(timeval, cc))+geom_point(aes(color = state))+geom_smooth(aes(color = state,group = state))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd))+facet_wrap(~geneset, nrow = 1)+theme_calc(base_size = 14)

collect = collect[order(collect$type),]
collect$geneset.size = rep(lengths(genesets), 4)
ggplot(collect[grepl("1|2|3|4|5|6|8|9", collect$geneset)&grepl("3hr", collect$type),], aes(fill=geneset.size, y=cc, x= as.factor(geneset) )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  # geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 2.7)+theme_bw(base_size = 20)

ggplot(collect[grepl("10|11|13|15", collect$geneset),], aes(fill=timeval, y=cc, x= timeval )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ 
  # geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+ylim(0, 2.7)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  geom_smooth(aes(group = as.factor(geneset.state), color = as.factor(geneset), linetype = as.factor(state)), size = 2)+
  # geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_point(size =2)+
  theme_bw(base_size = 16)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  geom_smooth(data = collect0[grepl("10|11|13|15", collect0$geneset),], aes(group = as.factor(geneset.state), color = as.factor(geneset)), linetype="dotted",size = 1)+
  # geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_point(size =2)+
  theme_bw(base_size = 16)

p1=ggplot(collect[grepl("1|2|8", collect$geneset),], aes(fill=type, y=cc, x= timeval )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ #geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 1.5)+  theme_bw(base_size = 20)
p2=ggplot(collect[grepl("3|9", collect$geneset),], aes(fill=type, y=cc, x= timeval )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ #geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 1.5)+  theme_bw(base_size = 20)
p3=ggplot(collect[grepl("4|5|6", collect$geneset),], aes(fill=type, y=cc, x= timeval )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ #geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_line(aes(group = as.factor(geneset), color = as.factor(geneset)), size = 2)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 1.5)+  theme_bw(base_size = 20)
p1/p2/p3

#pairwise biofunctions USE THIS----
p1=ggplot(collect[grepl("1|2|8", collect$geneset)&grepl("3hr", collect$type),], aes( y=cc, x= as.factor(geneset) )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ #geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity",position="dodge", size = 2, aes(fill = as.factor(geneset)))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 1.5)+  theme_bw(base_size = 20)
p2=ggplot(collect[grepl("3|9", collect$geneset)&grepl("3hr", collect$type),], aes( y=cc, x= as.factor(geneset) )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ #geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity",position="dodge", size = 2, aes(fill = as.factor(geneset)))+
  scale_fill_manual(values = c("3"="violet", "9"="#00BFC4"))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.05))+
  ylim(0, 1.5)+  theme_bw(base_size = 20)
p1/p2

ggplot(collect[grepl("13", collect$geneset),], aes(fill=as.factor(geneset), y=cc, x= type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,  position=position_dodge(0.9))+
  ylim(0, 2.7)+theme_bw(base_size = 20)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_collectdimensionbest_subsets1.12.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_collectdimensionbest_subsets1.12.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_collectdimensionbest_subsets1.12.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_PMexpts_Feb2021_rmUnstim_collectdimensionbest_subsets1.12.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_geneset13_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_geneset13_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_geneset13_subsets.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_geneset15.16.17_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_geneset15.16.17_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_geneset15.16.17_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_PMexpts_Feb2021_geneset15.16.17_subsets.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_biological_subsets2.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_biological_subsets2.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_biological_subsets2.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_PMexpts_Feb2021_rmUnstim_biological_subsets2.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M0_rep2only_biological_subsets3.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_biological_subsets3.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_biological_subsets3.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_biological_subsets4_GRS.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_biological_subsets4_GRS.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_biological_subsets4_GRS.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M0_rep2only_5stim_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M0_2019samples.3hrXtra_5stim_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_5stim_biological_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_5stim_biological_subsets.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_top3hr_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M0_rep2only_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_BMDM2_WT_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_BMDM2_MM_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_PMexpts_Feb2021_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_consensusMI_subsets_bootstrap.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_consensusMI_subsets_bootstrap.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_consensusMI_subsets_bootstrap.txt", quote = F, sep = "\t",row.names = F)

# write.table(collect, "./infotheo/SLEMI_M0_rep2only_5stim_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M0_2019samples.3hrXtra_5stim_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M1_IFNg_5stim_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)
# write.table(collect, "./infotheo/SLEMI_M2_IL4_gt80_5stim_consensusMI_subsets.txt", quote = F, sep = "\t",row.names = F)


collect.M0= read.delim("./infotheo/SLEMI_M0_rep2only_5stim_consensusMI_subsets.txt");colnames(collect.M0)[1]="time"
collect.M0r= read.delim("./infotheo/SLEMI_M0_2019samples.3hrXtra_5stim_consensusMI_subsets.txt");colnames(collect.M0r)[1]="time"
collect.M1= read.delim("./infotheo/SLEMI_M1_IFNg_5stim_consensusMI_subsets.txt");colnames(collect.M1)[1]="time"
collect.M2= read.delim("./infotheo/SLEMI_M2_IL4_gt80_5stim_consensusMI_subsets.txt");colnames(collect.M2)[1]="time"

collect.M0= read.delim("./infotheo/SLEMI_M0_rep2only_consensusMI_subsets_bootstrap.txt");colnames(collect.M0)[1]="time"
collect.M1= read.delim("./infotheo/SLEMI_M1_IFNg_consensusMI_subsets_bootstrap.txt");colnames(collect.M1)[1]="time"
collect.M2= read.delim("./infotheo/SLEMI_M2_IL4_gt80_consensusMI_subsets_bootstrap.txt");colnames(collect.M2)[1]="time"

collect.M0= read.delim("./infotheo/SLEMI_M0_rep2only_collectdimensionbest_subsets1.12.txt");colnames(collect.M0)[1]="time"
collect.M1= read.delim("./infotheo/SLEMI_M1_IFNg_collectdimensionbest_subsets1.12.txt");colnames(collect.M1)[1]="time"
collect.M2= read.delim("./infotheo/SLEMI_M2_IL4_gt80_collectdimensionbest_subsets1.12.txt");colnames(collect.M2)[1]="time"
collect.PM= read.delim("./infotheo/SLEMI_PMexpts_Feb2021_rmUnstim_collectdimensionbest_subsets1.12.txt")

collect.WT = read.delim("./infotheo/SLEMI_BMDM2_WT_consensusMI_subsets.txt");colnames(collect.WT)[1]="time"
collect.MM = read.delim("./infotheo/SLEMI_BMDM2_MM_consensusMI_subsets.txt");colnames(collect.MM)[1]="time"
collect.PM = read.delim("./infotheo/SLEMI_PMexpts_Feb2021_consensusMI_subsets.txt")

collect.M0= read.delim("./infotheo/SLEMI_M0_rep2only_geneset15.16.17_subsets.txt");colnames(collect.M0)[1]="time"
collect.M0r= read.delim("./infotheo/SLEMI_M0_2019samples.3hrXtra_geneset15.16.17_subsets.txt");colnames(collect.M0r)[1]="time"
collect.M1= read.delim("./infotheo/SLEMI_M1_IFNg_geneset15.16.17_subsets.txt");colnames(collect.M1)[1]="time"
collect.M2= read.delim("./infotheo/SLEMI_M2_IL4_gt80_geneset15.16.17_subsets.txt");colnames(collect.M2)[1]="time"
collect.PM= read.delim("./infotheo/SLEMI_PMexpts_Feb2021_geneset15.16.17_subsets.txt");colnames(collect.M2)[1]="time"

collect.M0$type = "M0"
collect.M0r$type = "M0r"
collect.M1$type = "M1"
collect.M2$type = "M2"
collect.WT$type = "WT"
collect.MM$type = "MM"
collect_all = rbind(collect.M0, collect.M0r, collect.M1, collect.M2)
collect_all = rbind(collect.M0, collect.M1, collect.M2)
collect_all = rbind(collect.WT, collect.MM)

collect_all$time.val = as.numeric(gsub("hr", "", collect_all$time))
ggplot(collect_all[grepl("", collect_all$geneset)&!grepl("2|1|8|^5hr", collect_all$time)&grepl("", collect_all$type),], 
       aes( y=cc, x= geneset, fill=type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.9))+
  ylim(0, 2.7)+theme_bw(base_size = 20)#+facet_grid(~geneset)

#PM disease channel capacity----
collect_all = collect.PM
collect_all$type = factor(collect_all$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
ggplot(collect_all[grepl("15", collect_all$geneset)&grepl("", collect_all$type),], 
       aes( y=cc, x= type, fill=type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.9))+
  ylim(0, 2.7)+theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(collect_all[grepl("3", collect_all$geneset)&grepl("", collect_all$type),], 
       aes( y=cc, x= type, fill=type )) + #facet_wrap(~cluster, scales = "free", nrow = 2)+
  ylab("channel capacity")+ geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.2,position=position_dodge(0.9))+
  ylim(0, 2.7)+theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# heatmap of geneset ----
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = readRDS("./output/macrophage_M1_IFNg_500genes_DBEC.rds")
macro = readRDS("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")

data.subset = as.data.frame(macro[["ISnorm"]]@data)
meta = macro@meta.data
meta = meta[order(meta$timept), ]
meta = meta[!grepl("24hr", meta$timept), ]
meta$stimulus <- factor(meta$stimulus, levels = c("Unstim", "LPS", "PIC", "IFNb", "P3CSK", "CpG", "TNF"))
meta = meta[order(meta$stimulus), ]
col_order = rownames(meta)

pheatmap(data.subset[genesets[[13]],col_order ], scale = "row", cluster_cols = F, 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-2,2,length=100),4),
         annotation_col = data.frame(macro@meta.data[,c(6,7)]),
         clustering_method = "ward.D2", show_colnames = F)

pheatmap(data.subset[genesets[[17]],col_order ], scale = "row", cluster_cols = F, cluster_rows = F, 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-2,2,length=100),4),
         annotation_col = data.frame(macro@meta.data[,c(6,7)]),
         clustering_method = "ward.D2", show_colnames = F)

##########################################################################
# channel capacity for each pairwise stimuli at each timepoint----
collect = data.frame()
list = c("CpG","P3CSK","LPS", "TNF", "IFNb",  "PIC") #
# list = c("P3CSK","LPS", "TNF", "IFNb",  "PIC") #


# for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
# for (i in c("0.25hr","1hr","3hr","8hr")){

macro = readRDS("output/macrophage_BMDM2_WT_MM_500genes_DBEC.rds")
macro = subset(macro, subset = type=="BMDM2_MM")
# for (i in c("BMDM2_MM")){

macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
macro = subset(macro, subset = timept=="8hr")

macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
macro = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
macro = subset(macro, subset = timept=="3hr")


for (i in c("3hr")){
  print(i)
  
    
  for (j in seq(1:(length(list)-1) )){
    for (k in seq(j+1,length(list)) ){
      
      for (gene in c(rownames(macro))){
      # for (gene in c("Cxcl10", "Ccl5","Tnf","Il6", "Irf7", "Cmpk2", "Socs3")){
      # for (g in (1:9)){
        
        skip_to_next <- FALSE
        print(gene)
        
      tryCatch(
          {
            
        
      print(list[j])
      print(list[k])
      # print(gene)
      macro.subset = subset(x = macro, subset = stimulus==list[j]|stimulus == list[k])
      
      data = macro.subset[["ISnorm"]]@data
      data = data.frame(data)
      meta = macro.subset@meta.data
      colnames(data) = paste0(meta$stimulus, "_",meta$timept)
      # rownames(data) = gsub("-ENSMUST..*","", macro.subset[["SCT"]]@counts@Dimnames[[1]])
      # rownames(data) = gsub("..NM..*","", rownames(data))
      my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
      my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
      rownames(my.dataframe) = seq(1:nrow(my.dataframe))
      # str(my.dataframe)

      
      # genes = genesets[[g]]
      # print(genes)
      # my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
      my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("label", gene)]
      str(my.dataframe.subset)
      
      output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", response = colnames(my.dataframe.subset)[-1],testing=F, boot_prob = 0.5, boot_num = 3, testing_cores = 4)
      sd = 0; #sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
      
      #--------------------------------mi using SLEMI -----------------------
      
      #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
      # output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                              # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i,"_",list[j],".",list[k]), 
                                              # testing = F)
      
      #just a few genes
      # output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = gene,
      #                                         # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i,"_",list[j],".",list[k]),
      #                                         testing = F)
      
      tmp = data.frame(time = i, stim1 = list[j], stim2 = list[k], cc = output_capacity$cc, sd = sd, gene = gene)
      # tmp = data.frame(time = i, stim1 = list[j], stim2 = list[k], cc = output_capacity$cc,  geneset=g)
      collect = rbind(collect, tmp)

      # output_mi  <- mi_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
      #                              paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/mi_", i),
      #                              pinput=rep(1/6,6))
          
      }, error = function(e) { skip_to_next <<- TRUE})
        
      if(skip_to_next) { next }  
        
      
      }
    }
  }
}
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_500genes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_selectgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M1_IFNg_selectgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M2_IL4_gt80_selectgenes.txt",quote=F, sep="\t",row.names = F)

# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_selectgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_selectgenes2.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_selectgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_selectgenes2.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_selectgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_selectgenes2.txt",quote=F, sep="\t",row.names = F)

# write.table(collect, "./infotheo/channel_capacity_pairwise_BMDM2_MM_ISnorm_allgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_BMDM2_WT_ISnorm_allgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes_8hr.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_allgenes.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_allgenes.txt",quote=F, sep="\t",row.names = F)


# write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_biological_categories3.txt",quote=F, sep="\t",row.names = F)

#pairwise plot M1M2M0---------------------------
collect.M0 = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_selectgenes.txt")
collect.M1 = read.delim("./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_selectgenes.txt")
collect.M2 = read.delim("./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_selectgenes.txt")

collect.M0$type = "M0"
collect.M1$type = "M1"
collect.M2$type = "M2"

collect_all = rbind(collect.M0, collect.M1, collect.M2)
collect_all$time.val = as.numeric(gsub("hr", "", collect_all$time))
collect_all$pair = paste0(collect_all$stim1,"_", collect_all$stim2)
ggplot(collect_all[!grepl("24hr", collect_all$time)&grepl("Cxcl10$|Ccl5|Tnf$|Il6$", collect_all$gene)&grepl("LPS_TNF|P3CSK_LPS|TNF_PIC|LPS_PIC", collect_all$pair),], aes(time.val, cc))+
  # geom_line(aes(group = type, color = type), size = 1)+
  geom_smooth(aes(group = type, color = type), method = "loess", se=F)+
  ylab("max MI (bits)")+facet_grid(pair~gene)+
  geom_point()+theme_calc(base_size = 14)+xlab("Time (hours)")


# for the 3hr timept----------------------------
collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt")
collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes_8hr.txt")
collect = read.delim("./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_allgenes.txt")
collect = read.delim("./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_allgenes.txt")
collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_biological_categories3.txt")

#pairwise plot 3hrs-----
collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt")
collect$pair = paste0(collect$stim1,"_", collect$stim2)
collect$pair = factor(collect$pair, levels= c("CpG_LPS","CpG_P3CSK", "P3CSK_LPS",
                                              "CpG_TNF","CpG_IFNb","LPS_TNF","LPS_IFNb", "P3CSK_TNF","P3CSK_IFNb",
                                              "CpG_PIC","LPS_PIC","P3CSK_PIC",
                                              "TNF_PIC", "IFNb_PIC",
                                              "TNF_IFNb"))
collect$pair = factor(collect$pair, levels= c("CpG_IFNb","CpG_LPS","CpG_P3CSK","CpG_PIC","CpG_TNF", 
                                              "LPS_IFNb", "P3CSK_IFNb","IFNb_PIC","TNF_IFNb",
                                              "P3CSK_LPS","LPS_PIC","LPS_TNF",
                                              "P3CSK_PIC","P3CSK_TNF",
                                              "TNF_PIC"))
collect$color = as.factor(ifelse(grepl("CpG",collect$pair), "CpG", 
                       ifelse(grepl("IFNb",collect$pair), "IFNb", 
                              ifelse(grepl("LPS",collect$pair), "LPS", 
                                     ifelse(grepl("P3CSK",collect$pair), "P3CSK", 
                                            ifelse(grepl("PIC",collect$pair), "PIC",
                                                   ifelse(grepl("TNF",collect$pair), "TNF",0)
                                     )))) ))
colors_list = c("Unstim" = "gray", "CpG"="#F8766D", "IFNb"="#B79F00","LPS"= "#00BA38","P3CSK"= "#00BFC4","PIC"= "#619CFF","TNF"= "#F564E3")
ggplot(collect[!grepl("IFNb", collect$pair),], aes((cc)))+geom_density(aes(group = pair, fill = color), alpha = 0.3)+
  scale_fill_manual(values = colors_list)+
  facet_wrap(~pair, ncol = 4, scales = "free_y")+theme_bw()

p2=ggplot(collect[grepl("Tnf$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p3=ggplot(collect[grepl("Cxcl10$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p1=ggplot(collect[grepl("Ccl5$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p1|p2|p3

p2=ggplot(collect[grepl("Cmpk2$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p3=ggplot(collect[grepl("Nfkbiz$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p1=ggplot(collect[grepl("Ifit3$", collect$gene) &grepl("3hr", collect$time),], aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(0,1)
p2|p1|p3

collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_500genes.txt")
collect$pair = paste0(collect$stim1,"_", collect$stim2)
collect = collect[grepl("3hr", collect$time),]
ggplot(collect, aes(pair, cc))+geom_bar(stat = "identity",position="dodge",(aes(fill = stim1)))+
  theme_bw(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#all timepoints
collect$geneset = as.factor(collect$geneset)
collect$timeval = as.numeric(gsub("hr", "",collect$time))

ggplot(collect, aes(time, cc, group= pair))+geom_point(aes(color = pair))+geom_line(aes(color = pair))+
  theme_bw(base_size = 16)
ggplot(collect[grepl("LPS_LPSlo|LPS_TNF|P3CSK_PIC|LPSlo_TNF|P3CSK_TNF", collect$pair),], aes(time, cc, group= pair))+
  geom_point(aes(color = pair),size=3)+geom_line(aes(color = pair),size = 1)+
  theme_bw(base_size = 16)+ylab("channel_capacity")
p1=ggplot(collect[grepl("LPS_TNF", collect$pair),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p2=ggplot(collect[grepl("CpG_LPS", collect$pair),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p3=ggplot(collect[grepl("CpG_PIC", collect$pair),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p1/p2/p3 

p1=ggplot(collect[grepl("LPS_TNF", collect$pair)&grepl("3|9", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p2=ggplot(collect[grepl("LPS_PIC", collect$pair)&grepl("3|9", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p3=ggplot(collect[grepl("IFNb_PIC", collect$pair)&grepl("3|9", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p4=ggplot(collect[grepl("CpG_P3CSK", collect$pair)&grepl("3|9", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p1/p2/p3/p4 

p1=ggplot(collect[grepl("LPS_TNF", collect$pair)&grepl("1|2|8", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p2=ggplot(collect[grepl("LPS_PIC", collect$pair)&grepl("1|2|8", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p3=ggplot(collect[grepl("IFNb_PIC", collect$pair)&grepl("1|2|8", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p4=ggplot(collect[grepl("CpG_P3CSK", collect$pair)&grepl("1|2|8", collect$geneset),], aes(timeval, cc, group= geneset))+
  geom_point(size=2)+geom_line(aes(color = geneset),size = 2)+
  theme_bw(base_size = 16)+ylab("channel capacity")+ylim(c(0,1.0))
p1/p2/p3/p4 

#USE THIS 3hrs------
collect = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_biological_categories3.txt")
collect$pair = paste0(collect$stim1,"_", collect$stim2)
collect = collect[grepl("3hr", collect$time),]
collect$geneset = as.character(collect$geneset)
color_list = c("3"="violet", "9"="#00BFC4")
p1=ggplot(collect[grepl("LPS_TNF", collect$pair)&grepl("3|9", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  scale_fill_manual(values = color_list)+
  theme_bw(base_size = 16)+ylab("LPS_TNF")+ylim(c(0,1.0))+theme(legend.position = "None")
p2=ggplot(collect[grepl("LPS_PIC", collect$pair)&grepl("3|9", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  scale_fill_manual(values = color_list)+
  theme_bw(base_size = 16)+ylab("LPS_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p3=ggplot(collect[grepl("IFNb_PIC", collect$pair)&grepl("3|9", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes( fill = as.factor(geneset)),size = 2)+
  scale_fill_manual(values = color_list)+
  theme_bw(base_size = 16)+ylab("IFNb_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p4=ggplot(collect[grepl("TNF_PIC", collect$pair)&grepl("3|9", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  scale_fill_manual(values = color_list)+
  theme_bw(base_size = 16)+ylab("TNF_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p5=ggplot(collect[grepl("CpG_P3CSK", collect$pair)&grepl("3|9", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  scale_fill_manual(values = color_list)+
  theme_bw(base_size = 16)+ylab("CpG_P3CSK")+ylim(c(0,1.0))+theme(legend.position = "None")
p1|p2|p3|p4|p5

p1=ggplot(collect[grepl("LPS_TNF", collect$pair)&grepl("1|2|8", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  theme_bw(base_size = 16)+ylab("LPS_TNF")+ylim(c(0,1.0))+theme(legend.position = "None")
p2=ggplot(collect[grepl("LPS_PIC", collect$pair)&grepl("1|2|8", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  theme_bw(base_size = 16)+ylab("LPS_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p3=ggplot(collect[grepl("IFNb_PIC", collect$pair)&grepl("1|2|8", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  theme_bw(base_size = 16)+ylab("IFNb_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p4=ggplot(collect[grepl("TNF_PIC", collect$pair)&grepl("1|2|8", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  theme_bw(base_size = 16)+ylab("TNF_PIC")+ylim(c(0,1.0))+theme(legend.position = "None")
p5=ggplot(collect[grepl("CpG_P3CSK", collect$pair)&grepl("1|2|8", collect$geneset),], aes(as.factor(geneset), cc))+
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(geneset)),size = 2)+
  theme_bw(base_size = 16)+ylab("CpG_P3CSK")+ylim(c(0,1.0))+theme(legend.position = "None")
p1|p2|p3|p5 

tmp = readRDS("./infotheo/without0hr/cc_8hr_LPS.TNFoutput.rds")
wts = as.data.frame(tmp$model$wts)
rownames(wts) = c("Intercept", tmp$model$coefnames)
sumMN   <- summary(tmp$model)
coefMN <- coef(sumMN)

#################################################################################
# old extra---------------
devtools::install_github("skinnider/dismay")
library(dismay)

mat = matrix(rnorm(100), ncol = 10, dimnames = list(paste0("cell_", 1:10), 
                                                    paste0("gene_", letters[1:10])))
mat[mat < 0] = 0
tc = dismay(mat, 'jaccard')
cos = dismay(mat, 'cosine')
mi = dismay(mat, 'MI')


###
# information

library(infotheo)
library(philentropy)


# Jensen-Shannon Divergence between P and Q
P <- 1:10/sum(1:10)
Q <- 20:29/sum(20:29)
x <- rbind(P,Q)

# Jensen-Shannon Divergence between P and Q using different log bases
JSD(x, unit = "log2") # Default
JSD(x, unit = "log")
JSD(x, unit = "log10")