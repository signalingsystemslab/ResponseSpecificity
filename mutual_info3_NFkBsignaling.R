# mutual_info3_NFkBsignaling
# last mod. 5/20/21, ksheu

x.1 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_LPS100ng_756_smoothed.mat")
x.2 = readMat("F:/enhancer_dynamics/nfkb_trajectories_08142019/smoothed/nfkb_dynamics_TNF10ng_754_smoothed.mat")

if(1){
  mat.1 = data.frame(x.1$data); 
  colnames(mat.1) = paste0("X",seq(0, (ncol(mat.1)-1)*5, 5)) #timepoints every 5mins
  order = rownames(mat.1)[order(apply(mat.1[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
  mat.1 = mat.1[match(order, rownames(mat.1)), ]
  mat.1 = cbind(label = "LPS", mat.1)
  pheatmap(mat.1[,-1], cluster_cols = F, cluster_rows = F, show_rownames = F,
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
           breaks=c(0,seq(0.1,6.9,length=100),7))
  
  mat.2 = data.frame(x.2$data);
  colnames( mat.2) = paste0("X",seq(0, (ncol( mat.2)-1)*5, 5)) #timepoints every 5mins
  order = rownames( mat.2)[order(apply( mat.2[,c(1:20)], 1, max), decreasing = T)] #max within first 100 mins
  mat.2 =  mat.2[match(order, rownames( mat.2)), ]
  mat.2 = cbind(label = "TNF",  mat.2)
  pheatmap(mat.2[,-1], cluster_cols = F, cluster_rows = F, show_rownames = F,
           colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
           breaks=c(0,seq(0.1,4.9,length=100),5))
}
collect = data.frame()
for (i in 3:98){
  
  my.dataframe = rbind(mat.1, mat.2)
  my.dataframe = na.omit(my.dataframe[,c(1,i)]) #for timepoint, for timeseries use [,c(1,2:i)]
  output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                          # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                          testing = F, boot_num = 10, boot_prob = .8 ,testing_cores = 4) #without0hr
  
  tmp = data.frame(time = (i-1)*5, pair = "LPS_TNF", cc = output_capacity$cc)
  collect = rbind(collect, tmp)
}
view(collect)