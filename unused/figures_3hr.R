####functions---------------------------------------------------------------------
plot_pca=function (file, info.name, info.type, title = "", labels = TRUE, 
          PCx = "PC1", PCy = "PC2", pt.size=1, ellipse = F, conf = 0.95, density = F, 
          fliph = F, flipv = F) 
{
  require(ggplot2)
  require(ggpubr)
  require(vegan)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)]
  if (grepl("scores_VARIMAX.txt", file)) {
    PCx = gsub("PC", "V", PCx)
    PCy = gsub("PC", "V", PCy)
    if (fliph == T) {
      table[, PCx] = table[, PCx] * -1
    }
    if (flipv == T) {
      table[, PCy] = table[, PCy] * -1
    }
  }
  sdev_name = paste0(gsub("scores.txt", "", file), "sdev.txt")
  if (!grepl("scores_VARIMAX.txt", file)) {
    sdev = read.delim(paste0(gsub("scores.txt", "", file), 
                             "sdev.txt"))
    sdev$var = unlist(sdev^2)
    sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, 
                            digits = 2))
    rownames(sdev) = paste0("PC", seq(1, nrow(sdev)))
    pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + 
      geom_point(size = I(pt.size), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect()
            # axis.text.x = element_text(margin = margin(b = -2)), 
            # axis.text.y = element_text(margin = margin(l = -14))
            ) + 
      guides(color = guide_legend(title = "Type")) + labs(title = title, 
                                                          x = paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], 
                                                                     "%)"), y = paste0(PCy, " (", sdev$pve[match(PCy, 
                                                                                                                 rownames(sdev))], "%)")) + theme_bw(base_size = 18) + 
      if (labels == TRUE) {
        geom_text(data = table, mapping = aes(label = Score), 
                  check_overlap = TRUE, size = 3)
      }
  }
  else if (grepl("scores_VARIMAX.txt", file)) {
    PCx = gsub("PC", "V", PCx)
    PCy = gsub("PC", "V", PCy)
    pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + 
      geom_point(size = I(pt.size), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect() 
            # axis.text.x = element_text(margin = margin(b = -2)), 
            # axis.text.y = element_text(margin = margin(l = -14))
            ) + 
      guides(color = guide_legend(title = "Type")) + labs(title = title) + 
      theme_bw(base_size = 18) + if (labels == TRUE) {
        geom_text(data = table, mapping = aes(label = Score), 
                  check_overlap = TRUE, size = 3)
      }
  }
  else {
    pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + 
      geom_point(size = I(pt.size), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            # axis.text.x = element_text(margin = margin(b = -2)), 
            # axis.text.y = element_text(margin = margin(l = -14))
            ) + 
      guides(color = guide_legend(title = "Type")) + labs(title = title) + 
      theme_bw(base_size = 18) + if (labels == TRUE) {
        geom_text(data = table, mapping = aes(label = Score), 
                  check_overlap = TRUE, size = 3)
      }
  }
  if (ellipse == TRUE) {
    plot(table[, c(PCx, PCy)], main = title)
    ord = ordiellipse(table[, c(PCx, PCy)], table$type, 
                      kind = "sd", conf = conf)
    cov_ellipse <- function(cov, center = c(0, 0), scale = 1, 
                            npoints = 100) {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for (g in (droplevels(table$type))) {
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(table[table$type == 
                                                               g, ], cov_ellipse(ord[[g]]$cov, ord[[g]]$center, 
                                                                                 ord[[g]]$scale))), type = g))
    }
    pcx.y2 = pcx.y + geom_path(data = df_ell, aes(x = df_ell[, 
                                                             PCx], y = df_ell[, PCy], colour = type), size = 1, 
                               linetype = 1)
    print(pcx.y2)
  }
  else {
    print(pcx.y)
  }
  if (density == TRUE) {
    xplot <- ggdensity(table, PCx, fill = "type") + clean_theme()
    yplot <- ggdensity(table, PCy, fill = "type") + rotate() + 
      clean_theme()
    (ggarrange(xplot, NULL, pcx.y, yplot, ncol = 2, nrow = 2, 
               align = "hv", widths = c(2, 1), heights = c(1, 2), 
               common.legend = TRUE))
  }
  else {
    print(pcx.y)
  }
}


###
plot_pca_projection=function (file, rotated.file, info.name, info.type, info.name2, 
                              info.type2, title = "Projection", labels = F, PCx = "PC1", 
                              PCy = "PC2", pt.size=1, ellipse = F, conf = 0.95, fliph = F, flipv = F, 
                              save = F, savename = "") 
{
  require(ggplot2)
  require(vegan)
  pc.scores = read.table(file, header = TRUE, row.names = 1)
  pc.scores.reduced = pc.scores
  pc.scores.reduced$type = info.type[match(rownames(pc.scores.reduced), 
                                           info.name)]
  if (fliph == T) {
    pc.scores.reduced[, PCx] = pc.scores.reduced[, PCx] * 
      -1
  }
  if (flipv == T) {
    pc.scores.reduced[, PCy] = pc.scores.reduced[, PCy] * 
      -1
  }
  projected_data = read.table(rotated.file, header = T)
  projected_data.reduced = projected_data
  projected_data.reduced$type = info.type2[match((projected_data.reduced[, 
                                                                         1]), info.name2)]
  if (fliph == T) {
    projected_data.reduced[, PCx] = projected_data.reduced[, 
                                                           PCx] * -1
  }
  if (flipv == T) {
    projected_data.reduced[, PCy] = projected_data.reduced[, 
                                                           PCy] * -1
  }
  pcx.y <- ggplot(projected_data.reduced, aes_string(x = PCx, 
                                                     y = PCy)) + geom_point(size = I(pt.size), aes(color = factor(type))) + 
    theme(legend.position = "right", plot.title = element_text(size = 30), 
          legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
          axis.title = element_text(size = 30), legend.background = element_rect() 
          # axis.text.x = element_text(margin = margin(b = -2)), 
          # axis.text.y = element_text(margin = margin(l = -14))
          ) + 
    guides(color = guide_legend(title = "Type")) + labs(title = title) + 
    theme_bw(base_size = 18) + if (labels == TRUE) {
      geom_text(mapping = aes(label = rownames(projected_data.reduced)), 
                check_overlap = TRUE, size = 3)
    }
  pcx.y <- pcx.y + geom_point(data = pc.scores.reduced, aes_string(x = PCx, 
                                                                   y = PCy)) + geom_point(size = I(pt.size), aes(color = factor(type))) + 
    theme(legend.position = "right", plot.title = element_text(size = 30), 
          legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
          axis.title = element_text(size = 30), legend.background = element_rect() 
          # axis.text.x = element_text(margin = margin(b = -2)), 
          # axis.text.y = element_text(margin = margin(l = -14))
          ) + 
    labs(title = title) + theme_bw(base_size = 18)
  if (ellipse == TRUE) {
    plot(projected_data.reduced[, c(PCx, PCy)], main = title)
    ord = ordiellipse(projected_data.reduced[, c(PCx, PCy)], 
                      projected_data.reduced$type, kind = "sd", conf = conf)
    cov_ellipse <- function(cov, center = c(0, 0), scale = 1, 
                            npoints = 100) {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for (g in (droplevels(projected_data.reduced$type))) {
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(projected_data.reduced[projected_data.reduced$type == 
                                                                                g, ], cov_ellipse(ord[[g]]$cov, ord[[g]]$center, 
                                                                                                  ord[[g]]$scale))), type = g))
    }
    pcx.y2 = pcx.y + geom_path(data = df_ell, aes(x = df_ell[, 
                                                             PCx], y = df_ell[, PCy], colour = type), size = 1, 
                               linetype = 1)
    pcx.y2
    if (save == T) {
      png(paste0(savename, ".png"), width = 8, height = 8, 
          units = "in", res = 300)
      plot(pcx.y2)
      dev.off()
      (pcx.y2)
    }
    else {
      pcx.y2
    }
  }
  else {
    if (save == T) {
      png(paste0(savename, ".png"), width = 8, height = 8, 
          units = "in", res = 300)
      plot(pcx.y)
      dev.off()
    }
    pcx.y
  }
}

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

######### 'figures ---------------------------------------------------------------library(edgeR);library(mclust);library(pheatmap);library(e1071);library(ksheu.library1);library(ggplot2);library(svglite)
library(RColorBrewer);library(ConsensusClusterPlus);library(Matrix);library(matrixStats);library(openxlsx);
library(rTensor);library(reshape2);library(rgl);library(ggrepel);library(plot3D);library(plot3Drgl);library(gridExtra)
library(patchwork);library(plot3D);library(plot3Drgl);library(pheatmap);library(grid);library(umap);library(ggthemes)
library(ggpubr);library(ggridges);library(reshape2);library(Seurat);library(dplyr);library(ksheu.library1)
library(SLEMI);library(ggplot2);library(ggrepel);library(ksheu.library1);library(Seurat);library(pheatmap)
setwd("F://scRNAseq_macro/scRNAseq_macro/")
######## fig1
mat = read.delim("F://scRNAseq_macro/bulk_rnaseq/Cheng2017_induced.txt") #1502 induced genes
rownames(mat) = mat$gene
mat$gene = NULL

#reorder stimuli by pattern
mat.order = mat[,c(2:5, 
                   1,12:14, 1,30:32, 1,9:11, 1,45:47, 1,6:8, #PAM
                   1,18:20, 1,21:23, 1,24:26, 1,27:29, #TNF
                   1,15:17, #FLA
                   1,33:35, 1, 36:38, 1, 39:41, 1, 42:44)]
#plot only 6 stimuli
mat.order = mat[,c(2:5, 
                   1,24:26,
                   1,30:32, #IFNb
                   1,12:14, #LPS
                   1,6:8, #PAM
                   1,9:11, #PIC
                   1,27:29 )] #TNF
mat.order = cbind(BMDM_UT=rowMeans(mat.order[,1:4]), mat.order[, grepl("_3", colnames(mat.order))]) #only 3hrs
set.seed(1)
sc.inducedgenes = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
sc.inducedgenes = rownames(sc.inducedgenes@assays$RNA)
mat.scgenes = mat.order[rownames(mat.order) %in% sc.inducedgenes, ]
scaled = t(scale(t(mat.scgenes)))
k=3
kmcluster <- kmeans(na.omit(scaled), iter.max = 1000, centers = k, nstart = 1)
tmp = data.frame(kmcluster$cluster);tmp$gene = rownames(tmp)
annot_df = data.frame(gene = names(kmcluster$cluster), cluster = as.character(kmcluster$cluster), row.names = 1)
bulk.order = rownames(mat.scgenes)[order(kmcluster$cluster)]
count <- 0
for (i in c(1:k)) {
  count[i] <- length(kmcluster$cluster[kmcluster$cluster == 
                                         i])
}
rowseps <- cumsum(count)
pheatmap(mat.scgenes[bulk.order,], scale = "row", cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1.5,1.5,length=100),4),
         gaps_row = rowseps,
         show_colnames = T, show_rownames = F)

#pheatmap all--------------------------------------
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = readRDS("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds")
macro = readRDS("./output/macrophage_M1_IFNg_500genes_DBEC.rds")
macro = readRDS("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
Idents(macro) = macro$stimulus
macro = subset(macro, subset=timept=="3hr"|timept=="0.0hr")
data = GetAssayData(object = macro, assay = "ISnorm", slot = "data")
data <- as.data.frame( as.matrix(data))
meta = macro@meta.data
meta = meta[order(meta$timept), ]
meta$stimulus <- factor(meta$stimulus, levels = c("Unstim","IFNb", "LPS", "PIC", "TNF", "CpG", "P3CSK" ))
meta = meta[order(meta$stimulus), ]
count =0
for (i in c(names(table(meta$stimulus)))) {
  count[i] <- table(meta$stimulus)[names(table(meta$stimulus))==i]
}
colseps <- cumsum(count[-1])

col_order = rownames(meta)

sc.order = na.omit(rownames(data)[match(bulk.order, rownames(data))])
# p = pheatmap(data[sc.order,col_order], scale = "none", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
#              gaps_row = rowseps[-5],
#              colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
#              breaks=c(-4,seq(-1,1,length=100),4),
#              annotation_col = data.frame(macro@meta.data[,c(6,7)]), show_colnames = F)
# p
#for ISnorm.data
colors_list = list(stimulus = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
p= pheatmap(data[sc.order,col_order], scale = "row", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
            colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
            gaps_row = rowseps[-5],gaps_col = colseps[-7],
            breaks=c(-2,seq(-1,1,length=100),2),
            annotation_col = data.frame(macro@meta.data[,c(6,7)]), 
            annotation_colors = colors_list,
            show_colnames = F, show_rownames = F)
ggsave(p, filename = "heatmap_singlecell.ISnorm_M0_rep2only_3hr.png")
ggsave(p, filename = "heatmap_singlecell.ISnorm_M0_2019samps_3hr.png")
ggsave(p, filename = "heatmap_singlecell.ISnorm_M1_IFNg_3hr.png")
ggsave(p, filename = "heatmap_singlecell.ISnorm_M2_IL4_gt80_3hr.png")

# plot multiple Ridgeplot for ISnorm counts-------------------------------------------------------------
# plotted Tnf, Cxcl10, Il12b
macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
colors_list = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")
p1=RidgePlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"), sort = "increasing",
             features = c("Socs3"), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  scale_fill_manual(values= colors_list)+
  stat_summary(fun.y = median, geom='line',  size = 1, colour = "black") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p2=RidgePlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"),  sort = "increasing",
             features = c("Tnf"), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  scale_fill_manual(values= colors_list)+
  stat_summary(fun.y = median, geom='line', size = 1, colour = "black") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p3=RidgePlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"),  sort = "increasing",
             features = c("Cxcl10"), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
  scale_fill_manual(values= colors_list)+
  stat_summary(fun.y = median, geom='line', size = 1, colour = "black") +theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
p1/p2/p3


#plot fano factor a few genes-------------------------------
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
colors_list = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")

data = macro[["ISnorm"]]@data
data = data.frame(data)
meta = macro@meta.data
colnames(data) = paste0(meta$orig.ident, "_",meta$type, "_",meta$stimulus, meta$timept)
colnames(data) = paste0(meta$type, "_",meta$stimulus, meta$timept)
data.t = data.frame(gene = colnames(data), t(data))

###
# data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = mean)
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(var(x))}) #var
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(sd(x)/mean(x))^2}) #cv^2
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(var(x))/mean(x)}) #fano
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(mean(x))/var(x)}) #snr

mat = t(data.aggregate[,-1])
colnames(mat) = data.aggregate$Group.1
# mat = data.frame(mat[rownames(mat)%in% bulk$gene,])
mat = data.frame(mat)
mat[is.na(mat)] <- 0
mat.m = melt(data.frame(gene = rownames(mat), mat))

#plot var vs mean scatterplot
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = mean)
mat = t(data.aggregate[,-1])
colnames(mat) = data.aggregate$Group.1
# mat = data.frame(mat[rownames(mat)%in% bulk$gene,])
mat = data.frame(mat)
mat[is.na(mat)] <- 0
mat.m$mean = melt(data.frame(gene = rownames(mat), mat))$value
mat.select.m = mat.m[grepl("Ccl5|Cxcl10|Tnf$", mat.m$gene),]
num_genes = 3
mat.select.m$stimulus = c(rep("CpG",4*num_genes), rep("IFNb",4*num_genes),rep("LPS",4*num_genes),
                          rep("P3CSK",4*num_genes), rep("PIC",4*num_genes),rep("TNF",4*num_genes), rep("Unstim",num_genes))
mat.select.m$time = c(rep(c(rep(0.25,num_genes),rep(1,num_genes),rep(3,num_genes),rep(8,num_genes)), 6), rep(0, num_genes))
mat.select.m$gene = factor(mat.select.m$gene, levels = c("Ccl5", "Tnf", "Cxcl10"))
ggplot(mat.select.m[grepl("3hr", mat.select.m$variable),], aes(mean, value, color = stimulus))+
  geom_point(size = 5,alpha = 0.75)+scale_fill_manual(values= colors_list)+
  geom_abline(slope = 1, intercept = 0, linetype="dotted")+xlim(0,10)+ylim(0,9)+ylab("variance")+
  facet_wrap(~gene, scales = "fixed", ncol = 1)+theme_classic(base_size = 14)

#plot fano factor
# mat.m = read.delim("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt")
mat.select.m = mat.m[grepl("Ccl5|Cxcl10|Tnf$", mat.m$gene),]
num_genes = 3
mat.select.m$stimulus = c(rep("CpG",4*num_genes), rep("IFNb",4*num_genes),rep("LPS",4*num_genes),
                          rep("P3CSK",4*num_genes), rep("PIC",4*num_genes),rep("TNF",4*num_genes), rep("Unstim",num_genes))
mat.select.m$time = c(rep(c(rep(0.25,num_genes),rep(1,num_genes),rep(3,num_genes),rep(8,num_genes)), 6), rep(0, num_genes))
mat.select.m$gene = factor(mat.select.m$gene, levels = c("Ccl5", "Tnf", "Cxcl10"))
ggplot(mat.select.m[grepl("3hr", mat.select.m$variable),], aes(stimulus, value))+
  geom_bar(stat= "identity", aes(fill = stimulus))+
  geom_hline(yintercept = 1,linetype="dotted" )+
  scale_fill_manual(values= colors_list)+
  facet_wrap(~gene, scales = "free_y", ncol = 1)+theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  
#stim-specificity of the Fano factor----
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
colors_list = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")

data = macro[["ISnorm"]]@data
data = data.frame(data)
meta = macro@meta.data
colnames(data) = paste0(meta$orig.ident, "_",meta$type, "_",meta$stimulus, meta$timept)
colnames(data) = paste0(meta$type, "_",meta$stimulus, meta$timept)
data.t = data.frame(gene = colnames(data), t(data))
data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(var(x)/mean(x))}) #fano
mat = t(data.aggregate[,-1])
colnames(mat) = data.aggregate$Group.1
mat = data.frame(mat)
mat = na.omit(mat)
# mat = mat[rowSums((mat[, -1] == 0)) < ncol(mat[-1]), ]
pheatmap((mat[ ,grepl("3hr|0hr", colnames(mat))] ), scale = "none", clustering_method = "ward.D2", cutree_rows = 4,
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103), show_rownames = F,
         breaks=c(0,seq(0.5,1.5,length=100),5),
         na_col="white", main = "Fano Factor")
do_kmeans_clustering((mat[ ,grepl("3hr", colnames(mat))]), cluster_cols = F,k_clusters = 5,show_colnames = T, show_rownames = F)
# -----------------------write out the clusters for homer-----
#do kmeans
mat.test = mat[ ,grepl("3hr", colnames(mat))]
set.seed(1)
scaled = t(scale(t(mat.test)))
k = 3
kmcluster <- kmeans(scaled, iter.max = 1000, centers = k, 
                    nstart = 1)
mat.test <- cbind(mat.test, cluster = kmcluster$cluster)
mat.test <- mat.test[order(mat.test$cluster), ]
count <- 0
for (i in 1:k) {
  count[i] <- length(kmcluster$cluster[kmcluster$cluster == i])
}
rowseps <- cumsum(count)
pheatmap((mat.test[,-ncol(mat.test)]), scale = "row", clustering_method = "ward.D2",
         gaps_row = rowseps,
         # gaps_col = c(4,8,14),
         # colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu"))[2:11])(103),
         # breaks=c(-4,seq(-2,2,length=100),4),
         cluster_rows = F, cluster_cols = T, show_colnames = T, show_rownames = F)

#write out the clusters
induced <- rownames(mat.test)
count <- 0;
for(j in 1:k){
  print(j)
  count[j+1] <- length(kmcluster$cluster[kmcluster$cluster == j])
  clust = induced[(cumsum(count)[j]+1):cumsum(count)[j+1]  ]
  # write.table(clust, paste0("./output/HOMER/FanoFactor_",k,"clust",j,"_M0.txt"), sep = '\t', quote = F, row.names = F)
}



#compare to stim-specificity of means by reordering
pheatmap(mat[bulk.order , grepl("3hr", colnames(mat))], scale = "row", cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         breaks=c(-4,seq(-1.5,1.5,length=100),4),
         gaps_row = rowseps,
         show_colnames = T, show_rownames = F)
##############fig 2--------------------
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
pca.macro = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt")

macro = readRDS("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020.rds")
pca.macro = read.delim("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020_3hr_prcomp_scores.txt")

rownames(pca.macro)=pca.macro$Score
pca.macro=pca.macro[,-1]

collect_distances = data.frame()
nPCs = 20
for (m in c("3hr")){
  print(m)
  type = m
  library(fpc);library(philentropy)
  # pca.macro = as.data.frame(Embeddings(macro[["pca"]]))
  
  
  
  for (i in (c("LPS", "TNF","CpG", "P3CSK", "PIC", "IFNb"))){
    # for (i in (c("LPS", "TNF", "P3CSK", "PIC", "IFNb"))){
    for (j in (c("LPS", "TNF" ,"CpG","P3CSK", "PIC", "IFNb"))){
      
      print(i) 
      wanted.1 = rownames(macro@meta.data)[macro@meta.data$timept == m & macro@meta.data$stimulus == i ]
      wanted.2 = rownames(macro@meta.data)[macro@meta.data$timept == m & macro@meta.data$stimulus == j ]
      
      
      mu1 = apply(pca.macro[wanted.1, 1:nPCs], 2, mean)
      mu2 = apply(pca.macro[wanted.2, 1:nPCs], 2, mean)
      cov1 = cov(pca.macro[wanted.1, 1:nPCs])
      cov2 = cov(pca.macro[wanted.2, 1:nPCs])
      bd = bhattacharyya.dist(mu1, mu2, Sigma1 = cov1, Sigma2 = cov2)
      
      # dist = data.frame(condition = paste0(i, "vsNot",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("LPSvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("TNFvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("IFNbvs",i), bd = bd, type = type)
      dist = data.frame(condition = paste0(i, "vs",j),stim1=i, stim2=j, bd = bd, type = type)
      
      
      if (length(collect_distances)==0){
        collect_distances = dist
      }else{
        collect_distances = rbind(collect_distances, dist ) 
      }
      
    }
  }
}

collect_distances.mat = dcast(collect_distances, stim1~stim2, value.var="bd")
rownames(collect_distances.mat) = collect_distances.mat$stim1
collect_distances.mat = collect_distances.mat[,-1]
collect_distances.mat[upper.tri(collect_distances.mat)] <- NA
pheatmap(collect_distances.mat, 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),na_col="white",
         cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", border_color = "white")

collect_distances.m = melt(collect_distances)
ggplot(collect_distances.m[grepl("",collect_distances.m$condition),], aes(fill=type, y=value, x=condition)) + 
  geom_bar(position="dodge", stat="identity")+#facet_wrap(~cluster, scales = "free", nrow = 1)+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
collect_distances.m.agg = aggregate(collect_distances.m, by = list(collect_distances.m$type), FUN = mean)
ggplot(collect_distances.m.agg, aes(fill=Group.1, y=value, x=Group.1)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw(base_size = 14)

#PCA on M0 3hr--------------------------
PCA_from_file("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020_3hr.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020_3hr_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =F)


PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =T)
intersect_doPCA_from_file_and_project_second_dataset("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt",
                                                     "./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.txt", train_string = "2M03hrs",
                                                     center = T, scale = F)
plot_pca_projection("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt", 
                    "./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC_2M03hrs_prcomp_rotated.txt",
                    samptag.all$Cell_Index, samptag.all$stimulus, samptag.all$Cell_Index, samptag.all$type)

sdev = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_sdev.txt")
sdev$var = unlist(sdev^2)
sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
plot(sdev$pve)
sum(sdev$pve[1:20])

#UMAP of unscaled PCA scores------------
pc.scores = read.delim("./output/macrophage_M0_2019samples.3hrXtra_500genes_Dec2020_3hr_prcomp_scores.txt")
pc.scores = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt")
data = (pc.scores)
rownames(data)= data$Score
set.seed(123)
umap.projected = umap(data[,c(2:21)], n_neighbors = 30L,
                      min_dist = 0.3, seed.use = 42)
umap.projected.layout = data.frame(umap.projected$layout)
umap.projected.layout$stimulus = samptag.all$stimulus[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
umap.projected.layout$type = samptag.all$type[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
colors_list = c("Unstim" = "gray", "CpG"="#F8766D", "IFNb"="#B79F00","LPS"= "#00BA38","P3CSK"= "#00BFC4","PIC"= "#619CFF","TNF"= "#F564E3")
p1=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = stimulus),size=0.01)+theme_bw(base_size = 18)+xlab("UMAP1")+ylab("UMAP2")+
  scale_color_manual(values = colors_list)
p2=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = type),size=0.01)+theme_bw(base_size = 16)
p1|p2


# differential genes for 3hrs ----------------
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = subset(macro, subset = timept =="3hr")
Idents(object = macro) <- macro@meta.data$stimulus
all.markers <- FindAllMarkers(object = macro, test.use = "t", assay = "ISnorm")
top20 <- all.markers %>% group_by(cluster) %>% top_n(20, (avg_log2FC))
DoHeatmap(macro,  group.by = "stimulus", assay = "ISnorm", slot = "data", #cells = wanted,
          top20$gene, size = 4, angle = 90)
mat = macro[["ISnorm"]]@data[top20$gene,]
mat = mat[!duplicated(mat),]
dim(mat)
colors_list = list(stimulus = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
heat = pheatmap(mat, scale = "row", show_colnames = F,
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),
         breaks=c(-2,seq(-1,1,length=100),2),
         annotation_col = data.frame(macro@meta.data[,c(6,7)]), 
         annotation_colors = colors_list,
         cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2")

# heat refers to the original heatmap produced from the pheatmap() function
# kept.labels should be a vector of labels you wish to show
# repel.degree is a number in the range [0, 1], controlling how much the
# labels are spread out from one another
sigGenes_v = c("Mx1","Ifit3","Il6","Tnf","Socs3","Cxcl10", "Egr2","Il1a","Il1b", "Ccl5","Tnfaip3", "Dusp1", "Nfkbiz","Cmpk2")
add.flag(heat,
         kept.labels = sigGenes_v,
         repel.degree = 0)

# Figure 3 heatmap of top cc genes---------------------------------------
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
ggplot(collect_all[grepl("3hr", collect_all$time),], aes( cc))+ geom_density()+ #graph it USE THIS
  # geom_point(aes(color = time),position = position_jitter(seed = 1))+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+
  geom_vline(xintercept = 0.7, linetype="dotted", size = 1)+
  geom_vline(xintercept = (log10(6)/log10(2)), linetype="dashed", size = 1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0,2.7)) + scale_y_continuous(expand = c(0, 0))

#####USE THIS plot green cc heatmap=======-----------------
dcast = dcast(collect_all, gene~time, value.var = "cc")
rownames(dcast) = dcast$gene
p=pheatmap(na.omit(dcast[(dcast[,4] > .7),4, drop=F]), scale = "none", clustering_method = "ward.D2", 
           colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
           cluster_cols = F, cluster_rows = T, show_rownames = T)
row_order = as.character(na.omit(dcast$gene[(dcast[,4] > .7)])[p$tree_row$order])


# gene expression heatmap in same order---------------------
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
pheatmap(data[row_order,col_order], scale = "row", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
         gaps_col = colseps[-6],
         breaks=c(-2,seq(-1,1,length=100),2),
         annotation_col = data.frame(macro@meta.data[,c(7,6)]), 
         annotation_colors = colors_list,
         show_colnames = F, show_rownames = T)

# Figure 3 PCA for mean vs variance-------------------------------------
tmp = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt")
# tmp = tmp[grepl("Tnf$|Il6|Il1a|Il1b", tmp$gene),]
# write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_proinfl.txt", quote=F, sep="\t",row.names = F)
# tmp = tmp[grepl("Ifit3$|Mx1$|Mx2$|Oasl1", tmp$gene),]
# write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral.txt", quote=F, sep="\t",row.names = F)
# tmp = tmp[grepl("Bcl2a1a$|Bcl2a1d$|Bcl2l11$|Tnfaip3", tmp$gene),]
# write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_apoptosis.txt", quote=F, sep="\t",row.names = F)
# tmp = tmp[grepl("Gsr$|Sod2$|Gclm$|Gsta3$", tmp$gene),]
# write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antioxidant.txt", quote=F, sep="\t",row.names = F)

PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_proinfl.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_proinfl_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antioxidant.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antioxidant_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =F)

#subsets
tmp = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt")
cells.wanted = samptag.all$Cell_Index[grepl("LPS|TNF|PIC", samptag.all$stimulus)]
tmp = tmp[grepl("Tnf$|Il6$|Il1a$|Il1b$", tmp$gene), colnames(tmp) %in% c("gene", cells.wanted)]
write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF.txt", quote=F, sep="\t",row.names = F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =F)

tmp = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt")
cells.wanted = samptag.all$Cell_Index[grepl("LPS|PIC|IFNb|TNF", samptag.all$stimulus)]
tmp = tmp[grepl("Ifit3$|Mx1$|Mx2$|Oasl1$", tmp$gene), colnames(tmp) %in% c("gene", cells.wanted)]
write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC.txt", quote=F, sep="\t",row.names = F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=1,labels = F, ellipse =F)

#apoptosis and antiox
tmp = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt")
cells.wanted = samptag.all$Cell_Index[grepl("LPS|TNF", samptag.all$stimulus)]
tmp = tmp[grepl("Bcl2a1a$|Bcl2a1d$|Bcl2l11$|Tnfaip3$|Gsr$|Sod2$|Gclm$|Gsta3$", tmp$gene), colnames(tmp) %in% c("gene", cells.wanted)]
write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF.txt", quote=F, sep="\t",row.names = F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_LPS.TNF_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =T)

tmp = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt")
cells.wanted = samptag.all$Cell_Index[grepl("TNF|PIC|LPS", samptag.all$stimulus)]
tmp = tmp[grepl("Bcl2a1a$|Bcl2a1d$|Bcl2l11$|Tnfaip3$|Gsr$|Sod2$|Gclm$|Gsta3$", tmp$gene), colnames(tmp) %in% c("gene", cells.wanted)]
write.table(tmp, "./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC.txt", quote=F, sep="\t",row.names = F)
PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_antiviral_TNF.PIC_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.5,labels = F, ellipse =F)

# Fig 3 mean vs variance----
dist = read.delim("./output/macrophage_M2_IL4_gt80_500genes_calculatedistfeatures.txt")

############################################################################
# Figure 4-----------------------------
sdev = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_sdev.txt")
sdev$var = unlist(sdev^2)
sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
sum(sdev$pve[1:20])

#UMAP of unscaled PCA scores M0M1M2------------
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

# PCA_from_file("./output/macrophage_M0M1M2_combined_500genes_DBEC_1hr.txt")
PCA_from_file("./output/macrophage_M0M1M2_combined_500genes_DBEC_8hr.txt")

pc.scores = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt")
data = (pc.scores)
rownames(data)= data$Score
set.seed(1)
umap.projected = umap(data[,c(2:21)], n_neighbors = 30L,
                      min_dist = 0.3, seed.use = 42)
umap.projected.layout = data.frame(umap.projected$layout)
umap.projected.layout$stimulus = samptag.all$stimulus[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
umap.projected.layout$type = samptag.all$type[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
p1=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = stimulus),size=0.01)+theme_bw(base_size = 18)+xlab("UMAP1")+ylab("UMAP2")
p2=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = type),size=0.01)+theme_bw(base_size = 18)+xlab("UMAP1")+ylab("UMAP2")
p1|p2

# distance metric M0, M1, M2
macro = readRDS("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr.rds")
DimPlot(macro, reduction = "tsne", group.by = "stimulus")
DimPlot(macro, reduction = "tsne", group.by = "type")

pca.macro = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt")
rownames(pca.macro)=pca.macro$Score
pca.macro=pca.macro[,-1]

collect_distances = data.frame()
nPCs = 20
for (m in c("M0","M1_IFNg","M2_IL4")){
  print(m)
  type = m
  library(fpc);library(philentropy)
  # pca.macro = as.data.frame(Embeddings(macro[["pca"]]))
  
  
  
  for (i in (c("LPS", "TNF","CpG", "P3CSK", "PIC", "IFNb"))){
    # for (i in (c("LPS", "TNF", "P3CSK", "PIC", "IFNb"))){
    for (j in (c("LPS", "TNF" ,"CpG","P3CSK", "PIC", "IFNb"))){
      
      print(i) 
      wanted.1 = rownames(macro@meta.data)[macro@meta.data$type == m & macro@meta.data$stimulus == i ]
      wanted.2 = rownames(macro@meta.data)[macro@meta.data$type == m & macro@meta.data$stimulus == j ]
      
      
      mu1 = apply(pca.macro[wanted.1, 1:nPCs], 2, mean)
      mu2 = apply(pca.macro[wanted.2, 1:nPCs], 2, mean)
      cov1 = cov(pca.macro[wanted.1, 1:nPCs])
      cov2 = cov(pca.macro[wanted.2, 1:nPCs])
      bd = bhattacharyya.dist(mu1, mu2, Sigma1 = cov1, Sigma2 = cov2)
      
      # dist = data.frame(condition = paste0(i, "vsNot",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("LPSvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("TNFvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("IFNbvs",i), bd = bd, type = type)
      dist = data.frame(condition = paste0(i, "vs",j),stim1=i, stim2=j, bd = bd, type = type)
      
      
      if (length(collect_distances)==0){
        collect_distances = dist
      }else{
        collect_distances = rbind(collect_distances, dist ) 
      }
      
    }
  }
}
collect_distances.m = melt(collect_distances)
# collect_distances.m$type = factor(collect_distances.m$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
ggplot(collect_distances.m[grepl("",collect_distances.m$condition),], aes(fill=type, y=value, x=condition)) +
  geom_bar(position="dodge", stat="identity")+#facet_wrap(~cluster, scales = "free", nrow = 1)+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


collect_distances.mat.M0 = dcast(collect_distances[grepl("M0",collect_distances$type),], stim1~stim2, value.var="bd")
collect_distances.mat.M0 = collect_distances.mat.M0[,-1]
collect_distances.mat = dcast(collect_distances[grepl("M1",collect_distances$type),], stim1~stim2, value.var="bd")
rownames(collect_distances.mat) = collect_distances.mat$stim1
collect_distances.mat = collect_distances.mat[,-1]
collect_distances.mat = ((collect_distances.mat - collect_distances.mat.M0)/collect_distances.mat.M0)*100
collect_distances.mat[upper.tri(collect_distances.mat)] <- NA
pheatmap(collect_distances.mat, 
         colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(103),na_col="white",
         # breaks=c(0,seq(.1,20,length=100),28),
         # breaks=c(-75,seq(-60,60,length=100),75),
         show_rownames = F,
         cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", border_color = "white")


# figure 4 pairwise 3genes, M1/M2 difference-----------------------------------------
collect.M0 = read.delim("./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_selectgenes2.txt")
collect.M1 = read.delim("./infotheo/channel_capacity_pairwise_M1_IFNg_ISnorm_selectgenes2.txt")
collect.M2 = read.delim("./infotheo/channel_capacity_pairwise_M2_IL4_gt80_ISnorm_selectgenes2.txt")
collect = cbind(collect.M0, M1=collect.M1, M2=collect.M2)
collect$M1diff = collect$M1.cc-collect$cc
collect$M2diff = collect$M2.cc-collect$cc
collect$M1diff.sd = sqrt((collect$M1.sd)^2 + (collect$sd)^2)
collect$M2diff.sd = sqrt((collect$M2.sd)^2 + (collect$sd)^2)

#pairwise plot 3hrs
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
p2=ggplot(collect[grepl("Tnf$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M1diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M1diff-M1diff.sd, ymax=M1diff+M1diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p3=ggplot(collect[grepl("Cxcl10$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M1diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M1diff-M1diff.sd, ymax=M1diff+M1diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p1=ggplot(collect[grepl("Ccl5$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M1diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M1diff-M1diff.sd, ymax=M1diff+M1diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p1|p2|p3

#m2 diff
p2=ggplot(collect[grepl("Tnf$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M2diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M2diff-M2diff.sd, ymax=M2diff+M2diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p3=ggplot(collect[grepl("Cxcl10$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M2diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M2diff-M2diff.sd, ymax=M2diff+M2diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p1=ggplot(collect[grepl("Ccl5$", collect$gene) &grepl("3hr", collect$time),], aes(pair, M2diff))+geom_bar(stat = "identity",position="dodge",(aes(fill = color)))+
  scale_fill_manual(values = colors_list)+
  geom_errorbar(aes(ymin=M2diff-M2diff.sd, ymax=M2diff+M2diff.sd), width=.5,position=position_dodge(0.05))+
  theme_classic(base_size = 16) +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "None")+ylim(-0.75,0.75)
p1|p2|p3


# figure 4 mutual info density plot-------
collect_all.M0 = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
collect_all.M1 = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt")
collect_all.M2 = read.delim("./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt")

collect_all.M0$type = "M0"
collect_all.M1$type = "M1"
collect_all.M2$type = "M2"

collect_all = rbind(collect_all.M0, collect_all.M1, collect_all.M2)

#graph it USE THIS
ggplot(collect_all[grepl("3hr", collect_all$time),], aes( (cc)))+ 
  geom_density(aes(color = type,fill=type),alpha=0.2, size = 1)+ theme_bw(base_size = 16)+
  # geom_point(aes(color = time),position = position_jitter(seed = 1))+
  ggtitle("Single gene output")+
  facet_wrap(~type, ncol=3, scales = "fixed")+
  # geom_vline(xintercept = (log10(6)/log10(2)), linetype="dashed", size = 1)+
  geom_vline(xintercept = (collect_all$cc[grepl("M0", collect_all$type)&grepl("Cxcl10", collect_all$gene)&grepl("3hr", collect_all$time)]), color= "red", linetype="dashed", size = 1)+
  geom_vline(xintercept = (collect_all$cc[grepl("M1", collect_all$type)&grepl("Cxcl10", collect_all$gene)&grepl("3hr", collect_all$time)]), color= "green3", linetype="dashed", size = 1)+
  geom_vline(xintercept = (collect_all$cc[grepl("M2", collect_all$type)&grepl("Cxcl10", collect_all$gene)&grepl("3hr", collect_all$time)]), color= "blue", linetype="dashed", size = 1)

# biological functions---
library(clusterProfiler)
library(enrichplot)
# BiocManager::install("enrichTF")
# BiocManager::install("RcisTarget")
# To support paralell execution:
# BiocManager::install(c("doMC", "doRNG"))
# For the examples in the follow-up section of the tutorial:
# BiocManager::install(c("DT", "visNetwork"))
library(enrichTF);library(RcisTarget)
library(org.Mm.eg.db)
library(DOSE)

collect_all = read.delim("./infotheo/SLEMI_singlegene_M0M1M2_ISnorm_500genes.txt")
dcast0 = dcast(collect_all[!grepl("0.25|0.5hr|^5hr|1hr|8h|24hr", collect_all$time),], gene~type+time, value.var = "cc")
dcast0$M1_3hr.diff = dcast0$M0_3hr-dcast0$M1_3hr
dcast0$M2_3hr.diff = dcast0$M0_3hr-dcast0$M2_3hr
dcast0 = na.omit(dcast0)
dcast0$clusterM1 = ifelse(dcast0$M1_3hr.diff >0.1, "DOWN", #used 0.2 for BP
                          ifelse(dcast0$M1_3hr.diff <(-0.1), "UP",NA))
dcast0$clusterM2 = ifelse(dcast0$M2_3hr.diff >0.1, "DOWN",  #used 0.2 for BP
                          ifelse(dcast0$M2_3hr.diff <(-0.1), "UP",NA))
dcast0$groups = ifelse((dcast0$clusterM1=="DOWN" & dcast0$clusterM2=="DOWN"), "bothDOWN",
                       ifelse(dcast0$clusterM1=="DOWN", "M1onlyDOWN",
                       ifelse(dcast0$clusterM2=="DOWN" , "M2onlyDOWN",
                               "bothUP"
                              )))

genes.to.label = c("Cxcl10","Tgtp1","Rsad2","Irf7","Trim21","Ifi205","Gbp7", 
                   "Ccl5", "Tnf", "Tnfaip3","Gclm", "Il6",
                   "Tnfsf9","Bcl2a1a","Bcl2a1d","Socs3")
genes.to.label = c("Cxcl10","Rsad2","Irf7","Trim21","Icosl", "Cmpk2","Ifit3","Nfkbiz","Il4ra","Cxcl2",
                   "Ccl5", "Tnf", "Il6", 
                   "Tnfsf9","Socs3", "Icam1")
rownames(dcast0) = dcast0$gene
p1=ggplot(dcast0,  aes(M1_3hr.diff, M2_3hr.diff))+geom_point(aes(color = groups), size =2, alpha = 0.75)+
  # geom_point(data = subset(dcast0, subset = gene %in% genes.to.label),size = 2, aes(fill = M0_3hr)) +
  geom_point(data = subset(dcast0, subset = gene %in% genes.to.label),size = 2, color = "red",shape = 21) +
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+theme_bw(base_size = 14)
p1 <- LabelPoints(plot = p1, points = genes.to.label, color = "red", size = I(4),repel = T, xnudge=0.15,ynudge=0)
p1

#my genelist
gene <- dcast0$gene
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)

gene.df$cluster = dcast0$clusterM1[match(gene.df$SYMBOL, dcast0$gene)]
# compared heatmap clusters
ego =  clusterProfiler::compareCluster(geneClusters = ENTREZID~cluster, 
                                       data=gene.df,
                                       fun = "enrichGO", 
                                       OrgDb = org.Mm.eg.db, 
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
p1 <- dotplot(ego2, showCategory=3)
p1

gene.df$cluster = dcast0$clusterM2[match(gene.df$SYMBOL, dcast0$gene)]
# compared heatmap clusters
ego.M2 =  clusterProfiler::compareCluster(geneClusters = ENTREZID~cluster, 
                                       data=gene.df,
                                       fun = "enrichGO", 
                                       OrgDb = org.Mm.eg.db, 
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
ego2.M2 <- clusterProfiler::simplify(ego.M2, cutoff=0.7, by="p.adjust", select_fun=min)
p2 <- dotplot(ego2.M2, showCategory=3) 
p2
p1|p2

# write.table(data.frame(ego2@compareClusterResult),"F://scRNAseq_macro/SuppTables/TableS4_M1_GO.txt",sep = "\t",quote=F,row.names = F)
# write.table(data.frame(ego2.M2@compareClusterResult),"F://scRNAseq_macro/SuppTables/TableS4_M2_GO.txt",sep = "\t",quote=F,row.names = F)

#TF enrichment-----
# Load gene sets to analyze. e.g.:
gene.df$cluster = dcast0$clusterM1[match(gene.df$SYMBOL, dcast0$gene)]
geneList1 <- na.omit(gene.df$SYMBOL[gene.df$cluster=="DOWN"])
geneList2 <- na.omit(gene.df$SYMBOL[gene.df$cluster=="UP"])
gene.df$cluster = dcast0$clusterM2[match(gene.df$SYMBOL, dcast0$gene)]
geneList3 <- na.omit(gene.df$SYMBOL[gene.df$cluster=="DOWN"])
geneList4 <- na.omit(gene.df$SYMBOL[gene.df$cluster=="UP"])
geneLists <- list(M1_DOWN = geneList1,M1_UP = geneList2,M2_DOWN = geneList3, M2_UP = geneList4 )
head(geneLists$M1_DOWN)

for (i in seq(1:4)){
  print(i)
  write.table(geneLists[[i]], paste0("./output/channelcapacity_M1M2_difference_", names(geneLists[i]),".txt"), quote=F,sep="\t",row.names=F) 
}
# cd /mnt/f/scRNAseq_macro/scRNAseq_macro/output/
# findMotifs.pl channelcapacity_M1M2_difference_M1_DOWN.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/output/motif_cc/MotifResults1_M1_DOWN/ -start -1000 -end 100 -p 4
# findMotifs.pl channelcapacity_M1M2_difference_M1_UP.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/output/motif_cc/MotifResults2_M1_UP/ -start -1000 -end 100 -p 4
# findMotifs.pl channelcapacity_M1M2_difference_M2_DOWN.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/output/motif_cc/MotifResults3_M2_DOWN/ -start -1000 -end 100 -p 4
# findMotifs.pl channelcapacity_M1M2_difference_M2_UP.txt mouse /mnt/f/scRNAseq_macro/scRNAseq_macro/output/motif_cc/MotifResults4_M2_UP/ -start -1000 -end 100 -p 4

#plot motif output--------------------
files = list.files("F://scRNAseq_macro/scRNAseq_macro/output/motif_cc/", "knownResults.txt", recursive = T)
for (i in seq(1:4)){
  tmp = read.delim(paste0("F://scRNAseq_macro/scRNAseq_macro/output/motif_cc/",files[i])) #files[1]
  # tmp = tmp[tmp$q.value..Benjamini.<0.05,]
  
  frame = data.frame(tmp[,c(1,4)]) #plot lnpval
  rownames(frame) = make.unique(as.character(frame$Motif.Name))
  rownames(frame) = make.unique(gsub("/..*","", rownames(frame)))
  frame.truncate = frame#[which(frame$Log.P.value < 0),]
  frame.truncate$type = ifelse(grepl("(IRF)", rownames(frame.truncate)), "IRF", 
                               ifelse(grepl("(RHD)", rownames(frame.truncate)), "RHD", 
                                      ifelse(grepl("(bZIP)", rownames(frame.truncate)), "bZIP", 
                                             ifelse(grepl("(Zf)", rownames(frame.truncate)), "Zf", 
                                                    ifelse(grepl("(ETS|Ets)", rownames(frame.truncate)), "ETS", 
                                                           ".other")))))
  frame.truncate.500 = aggregate(frame.truncate[,c(2,3)], by = list(motif = frame.truncate$type), FUN = min)
  assign(paste0("p",i), ggplot(frame.truncate.500[!grepl("other", frame.truncate.500$motif),], aes(motif, Log.P.value*-1, fill = motif))+
           ylim(c(0,26))+
    geom_bar(stat = "identity", position = position_dodge())+ geom_hline(yintercept = 3, linetype="dashed")+
      theme_bw(base_size = 16)+theme(legend.position = "None"))
}
(p1|p2)/(p3|p4)



# geneLists <- GSEABase::GeneSet(genes, setName="geneListName") # alternative
# load files---
# # featherURL <- "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather" 
# # download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# 
# # Select motif database to use (i.e. organism and distance around TSS)
# data(motifAnnotations_hgnc)
# motifRankings <- importRankings("./mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
# 
# # Search space: 10k bp around TSS - MOSUE
# motifRankings <- importRankings("m-tss-centered-10kb-7species.mc9nr.feather")
# # Load the annotation to human transcription factors
# data(motifAnnotations_hgnc)
# # Motif enrichment analysis:
# motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
#                                          motifAnnot=motifAnnotations_hgnc)

#fano factor of polarized macs----
mat.m0 = read.table("./output/macrophage_M0_rep2only_500genes_calculatedistfeatures.txt", header = T);mat.m0$type ="M0"
mat.m1 = read.table("./output/macrophage_M1_IFNg_500genes_calculatedistfeatures.txt", header = T);mat.m1$type ="M1"
mat.m2 = read.table("./output/macrophage_M2_IL4_gt80_500genes_calculatedistfeatures.txt", header = T);mat.m2$type ="M2"
mat = rbind(mat.m0, mat.m1, mat.m2)
ggplot(mat[grepl("3",mat$time),], aes(mean, fano))+geom_point(aes(color = stimulus))+theme_bw()+
  facet_wrap(stimulus~type, ncol = 6)
mat.dcast = dcast(mat[grepl("3",mat$time),], gene~type+stimulus, value.var = "fano")
rownames(mat.dcast) = mat.dcast$gene
pheatmap(mat.dcast[ ,!grepl("gene", colnames(mat.dcast))], scale = "row", clustering_method = "ward.D2", cutree_rows = 4,
         colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103), show_rownames = F,
         # breaks=c(0,seq(0.01,2,length=100),5),
         na_col="white", main = "Fano Factor")
mat.dcast = mat.dcast[rowSums((mat.dcast[, -1] == 0)) < ncol(mat.dcast[-1]), ]
do_kmeans_clustering(mat.dcast[ ,-1], cluster_cols = T, #colseps  = c(6,12),
                     show_colnames = T, show_rownames = F)


mat.dcast.diff = cbind(gene = mat.dcast$gene, (mat.dcast[,c(8:13)]-mat.dcast[,c(2:7)])*-1,(mat.dcast[,c(14:19)]-mat.dcast[,c(2:7)])*-1 )
rownames(mat.dcast.diff) = mat.dcast.diff$gene
clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
clusters$clusters = gsub("mB","mC", clusters$manual_fix)
clusters$clusters = ifelse(clusters$clusters=="mA","AP1",
                       ifelse(clusters$clusters=="mC","NFkB",
                              ifelse(clusters$clusters=="mD","NFkB&p38",
                                     ifelse(clusters$clusters=="mE","NFkB|IRF","IRF")))) 
clusters$clusters = factor(clusters$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))   

mat.dcast.diff$cluster = clusters$clusters[match(mat.dcast.diff$gene, clusters$gene)]
mat.dcast.diff = cbind(mat.dcast.diff, mat.dcast[,c(2:7)])
mat.dcast.diff$label = ifelse(grepl("Tnf$|Cxcl10$|Ccl5$|Socs3$|Ifna4a|Ccl1$",mat.dcast.diff$gene),mat.dcast.diff$gene,"" )
mat.dcast.diff$label = ifelse(apply(mat.dcast.diff[,c(2:13)]>2, 1, any) ,mat.dcast.diff$gene,"" )

p1=ggplot(mat.dcast.diff, aes(M1_CpG, M2_CpG))+geom_point(aes(color = M0_CpG),size = 2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+ #ylim(-4,4)+xlim(-3.5,3)+ 
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
p2=ggplot(mat.dcast.diff, aes(M1_IFNb, M2_IFNb))+geom_point(aes(color = M0_IFNb),size = 2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+ #ylim(-4,4)+xlim(-3.5,3)+  
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
p3=ggplot(mat.dcast.diff, aes(M1_LPS, M2_LPS))+geom_point(aes(color = M0_LPS),size = 2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+ #ylim(-4,4)+xlim(-3.5,3)+ 
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
p4=ggplot(mat.dcast.diff, aes(M1_P3CSK, M2_P3CSK))+geom_point(aes(color = M0_P3CSK),size = 2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+ #ylim(-4,4)+xlim(-3.5,3)+  
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
p5=ggplot(mat.dcast.diff, aes(M1_PIC, M2_PIC))+geom_point(aes(color = M0_PIC),size = 2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+ #ylim(-4,4)+xlim(-3.5,3)+  
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
p6=ggplot(mat.dcast.diff, aes(M1_TNF, M2_TNF))+geom_point(aes(color = M0_TNF),size =2, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")+#ylim(-4,4)+xlim(-3.5,3)+
  scale_color_gradient(low="dodgerblue2",  high="red")#+geom_text_repel(mapping = aes(label = label),max.overlaps=100, size = 3,box.padding = 0.5)
(p1|p2|p3)/(p4|p5|p6)
ggplot(mat.dcast.diff, aes(M1_TNF, M2_TNF))+geom_point(aes(color = M0_TNF),size =1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_color_gradient(low="dodgerblue2",  high="red")
# write.table(mat.dcast.diff, "./output/macrophage_M0M1M2_FanoFactor_difference.txt",sep = "\t",quote = F, row.names = F)

cor.test(mat.dcast.diff$M1_CpG, mat.dcast.diff$M2_CpG)
cor.test(mat.dcast.diff$M1_IFNb, mat.dcast.diff$M2_IFNb)
cor.test(mat.dcast.diff$M1_LPS, mat.dcast.diff$M2_LPS)
cor.test(mat.dcast.diff$M1_P3CSK, mat.dcast.diff$M2_P3CSK)
cor.test(mat.dcast.diff$M1_PIC, mat.dcast.diff$M2_PIC)
cor.test(mat.dcast.diff$M1_TNF, mat.dcast.diff$M2_TNF)


#by cluster GRN
p1=ggplot(mat.dcast.diff, aes(M1_CpG, M2_CpG))+geom_point(aes(color =cluster),size = 1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
p2=ggplot(mat.dcast.diff, aes(M1_IFNb, M2_IFNb))+geom_point(aes(color = cluster),size = 1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
p3=ggplot(mat.dcast.diff, aes(M1_LPS, M2_LPS))+geom_point(aes(color = cluster),size = 1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
p4=ggplot(mat.dcast.diff, aes(M1_P3CSK, M2_P3CSK))+geom_point(aes(color = cluster),size = 1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
p5=ggplot(mat.dcast.diff, aes(M1_PIC, M2_PIC))+geom_point(aes(color = cluster),size = 1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
p6=ggplot(mat.dcast.diff, aes(M1_TNF, M2_TNF))+geom_point(aes(color = cluster),size =1, alpha = 0.5)+theme_bw(base_size = 11)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+theme(legend.position = "None")
(p1|p2|p3)/(p4|p5|p6)



#heatmaps
clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
anno = data.frame(gene = mat.dcast.diff$gene)
anno$clusters = clusters$manual_fix[match(anno$gene, clusters$gene)]
anno$clusters = gsub("mB","mC", anno$clusters)
anno$clusters = ifelse(anno$clusters=="mA","AP1",
                                 ifelse(anno$clusters=="mC","NFkB",
                                        ifelse(anno$clusters=="mD","NFkB&p38",
                                               ifelse(anno$clusters=="mE","NFkB|IRF","IRF")))) 
anno$clusters = factor(anno$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))   
rownames(anno)= anno$gene

pheatmap(mat.dcast.diff[,-1], cluster_cols = F, clustering_method = "ward.D2",
         colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(103),
         gaps_col = c(6), show_rownames = F,
         breaks=c(-3,seq(-2,2,length=100),3),
         annotation_row = anno[,2,drop = F])
do_kmeans_clustering(na.omit(mat.dcast.diff[,-c(1)]), k_clusters = 5,
                     annotation_col = anno[,2,drop = F])

pheatmap(mat[,-ncol(mat)], cluster_cols = F, clustering_method = "ward.D2",scale = "row",
         colorRampPalette((brewer.pal(n = 9, name ="RdBu")))(103),
         gaps_col = c(6,12), show_rownames = F,
         # breaks=c(0,seq(.1,3,length=100),5),
         annotation_row = anno[,2,drop = F])

#Figure 5------------------------------------------------

macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds");

#projected
pca.macro = read.delim("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC_2M0M1M23hrs_prcomp_rotated.txt")
rownames(pca.macro)=pca.macro$Sample
pca.macro=pca.macro[,-1]

collect_distances = data.frame()
nPCs = 20
for (m in c("PM_B6.LFD","PM_B6.old","PM_B6.HFD")){
  print(m)
  type = m
  library(fpc);library(philentropy)
  # pca.macro = as.data.frame(Embeddings(macro[["pca"]]))
  
  
  
  # for (i in (c("LPS", "TNF","CpG", "P3CSK", "PIC", "IFNb"))){
  for (i in (c("IFNb", "LPS", "P3CSK", "PIC", "TNF"))){
    for (j in (c("IFNb", "LPS", "P3CSK", "PIC", "TNF"))){
      
      print(i) 
      wanted.1 = rownames(macro@meta.data)[macro@meta.data$type == m & macro@meta.data$stimulus == i ]
      wanted.2 = rownames(macro@meta.data)[macro@meta.data$type == m & macro@meta.data$stimulus == j ]
      
      
      mu1 = apply(pca.macro[wanted.1, 1:nPCs], 2, mean)
      mu2 = apply(pca.macro[wanted.2, 1:nPCs], 2, mean)
      cov1 = cov(pca.macro[wanted.1, 1:nPCs])
      cov2 = cov(pca.macro[wanted.2, 1:nPCs])
      bd = bhattacharyya.dist(mu1, mu2, Sigma1 = cov1, Sigma2 = cov2)
      
      # dist = data.frame(condition = paste0(i, "vsNot",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("LPSvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("TNFvs",i), bd = bd, type = type)
      # dist = data.frame(condition = paste0("IFNbvs",i), bd = bd, type = type)
      dist = data.frame(condition = paste0(i, "vs",j),stim1=i, stim2=j, bd = bd, type = type)
      
      
      if (length(collect_distances)==0){
        collect_distances = dist
      }else{
        collect_distances = rbind(collect_distances, dist ) 
      }
      
    }
  }
}
collect_distances.m = melt(collect_distances)
# collect_distances.m$type = factor(collect_distances.m$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
ggplot(collect_distances.m[grepl("",collect_distances.m$condition),], aes(fill=type, y=value, x=condition)) +
  geom_bar(position="dodge", stat="identity")+#facet_wrap(~cluster, scales = "free", nrow = 1)+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

collect_distances.mat.M0 = dcast(collect_distances[grepl("LFD",collect_distances$type),], stim1~stim2, value.var="bd")
collect_distances.mat.M0 = collect_distances.mat.M0[,-1]
collect_distances.mat = dcast(collect_distances[grepl("HFD",collect_distances$type),], stim1~stim2, value.var="bd")
rownames(collect_distances.mat) = collect_distances.mat$stim1
collect_distances.mat = collect_distances.mat[,-1]
collect_distances.mat = ((collect_distances.mat - collect_distances.mat.M0)/collect_distances.mat.M0)*100
collect_distances.mat[upper.tri(collect_distances.mat)] <- NA
pheatmap(collect_distances.mat, 
         colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(103),na_col="white",
         # breaks=c(0,seq(.1,20,length=100),28),
         breaks=c(-50,seq(-20,20,length=100),50),
         show_rownames = F,
         cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", border_color = "white")

# most differential genes disease ------
macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds");
macro = subset(macro, subset = type =="PM_B6.LFD")
Idents(object = macro) <- macro@meta.data$stimulus
all.markers <- FindAllMarkers(object = macro, test.use = "t", assay = "ISnorm")
top20 <- all.markers %>% group_by(cluster) %>% top_n(5, (avg_log2FC))
# DoHeatmap(macro,  group.by = "stimulus", assay = "ISnorm", slot = "data", #cells = wanted,
#           top20$gene, size = 4, angle = 90)
mat = macro[["ISnorm"]]@data[top20$gene,]
mat = mat[!duplicated(mat),]
dim(mat)
colors_list = list(stimulus = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"),
                   type = c(PM_B6.LFD="#F8766D",PM_B6.old="#00BA38",PM_B6.HFD="#619CFF"))
heat = pheatmap(mat, scale = "row", show_colnames = F,
                colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),
                breaks=c(-2,seq(-1,1,length=100),2),
                annotation_col = data.frame(macro@meta.data[,c(5,6)]), 
                annotation_colors = colors_list,
                cluster_rows = T, cluster_cols = T, clustering_method = "ward.D2")
row.order = rownames(mat)[heat$tree_row$order]

#plot others disease
macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds");
macro = subset(macro, subset = type =="PM_B6.old")
mat = macro[["ISnorm"]]@data[top20$gene,]
mat = mat[!duplicated(mat),]
dim(mat)
heat1 = pheatmap(mat[row.order,], scale = "row", show_colnames = F,
                colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),
                breaks=c(-2,seq(-1,1,length=100),2),
                annotation_col = data.frame(macro@meta.data[,c(5,6)]), 
                annotation_colors = colors_list,
                cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2")

macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds");
macro = subset(macro, subset = type =="PM_B6.HFD")
mat = macro[["ISnorm"]]@data[top20$gene,]
mat = mat[!duplicated(mat),]
dim(mat)
heat2 = pheatmap(mat[row.order,], scale = "row", show_colnames = F,
                 colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),
                 breaks=c(-2,seq(-1,1,length=100),2),
                 annotation_col = data.frame(macro@meta.data[,c(5,6)]), 
                 annotation_colors = colors_list,
                 cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2")

macro = readRDS("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.rds")
mat = macro[["ISnorm"]]@data[top20$gene,]
mat = mat[!duplicated(mat),]
dim(mat)
heat3 = pheatmap(mat[row.order,], scale = "row", show_colnames = F,
                 colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(103),
                 breaks=c(-2,seq(-1,1,length=100),2),
                 annotation_col = data.frame(macro@meta.data[,c(5,6)]), 
                 annotation_colors = colors_list,
                 cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2")

# plot disease channel capacity differences--------
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_PMs_Feb2021_gt120_5stim_ISnorm_500genes.txt")
dcast = dcast(collect_all, gene~type, value.var = "cc")
rownames(dcast)= dcast$gene
dcast$old.diff = dcast$PM_B6.LFD-dcast$PM_B6.old
dcast$HFD.diff = dcast$PM_B6.LFD-dcast$PM_B6.HFD

genes.to.label = c("Cxcl10","AW112010","Mx1", "Ralgds", "Acod1","Phlda1","Pde4b", "Mx2","Swap70", "Tnfsf9",
                   "Slamf8","Cav1","Il4ra","Zfp36", "Icosl", "Cmpk2","Ifit3","Nfkbiz","Icam1","Zc3h12a","Peli1",
                   "Ccl5", "Tnf",  "Il6")
p1=ggplot(dcast, aes(PM_B6.LFD-PM_B6.old, PM_B6.LFD-PM_B6.HFD))+geom_point(aes(color = PM_B6.LFD), size =2)+
  geom_point(data = subset(dcast, subset = gene %in% genes.to.label),size = 4, aes(color = PM_B6.LFD)) + 
  geom_point(data = subset(dcast, subset = gene %in% genes.to.label),size = 4, color = "red",shape = 21) + 
    geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+theme_bw(base_size = 14)
p1 <- LabelPoints(plot = p1, points = genes.to.label, color = "red", size = I(5),repel = T, xnudge=0,ynudge=0)
p1


#figure 4 WT/MM UMAP------
macro = readRDS("output/macrophage_BMDM2_WT_MM_500genes_LPT.rds")
i = "BMDM2_WT"
wanted = rownames(macro@meta.data)[macro@meta.data$type==i]
list = as.list(wanted);  names(list) = macro@meta.data$stimulus[match(list, rownames(macro@meta.data))]
list = split(unlist(list, use.names = FALSE), rep(names(list), lengths(list))) #merge under unique names
p1=DimPlot(object = macro, reduction = "umap", cells.highlight = list, 
           # cols.highlight = c(rev(hue_pal() (length( unique(names(list)) )) ) ),
           cols.highlight = c( TNF= "#F564E3", PIC= "#619CFF",LPS= "#00BA38"),
           sizes.highlight = .1,
           group.by = "stimulus", pt.size = .1)+theme(legend.position = "None")
i = "BMDM2_MM"
wanted = rownames(macro@meta.data)[macro@meta.data$type==i]
list = as.list(wanted);  names(list) = macro@meta.data$stimulus[match(list, rownames(macro@meta.data))]
list = split(unlist(list, use.names = FALSE), rep(names(list), lengths(list))) #merge under unique names
p2=DimPlot(object = macro, reduction = "umap", cells.highlight = list, 
           # cols.highlight = c(rev(hue_pal() (length( unique(names(list)) )) ) ),
           cols.highlight = c( TNF= "#F564E3", PIC= "#619CFF",LPS= "#00BA38"),
           sizes.highlight = .1,
           group.by = "stimulus", pt.size = .1)+theme(legend.position = "None")
p1|p2

#figure 4 WT/MM single genes vs signaling------
genes.mi = read.delim("F://scRNAseq_macro/scRNAseq_macro/infotheo/SLEMI_singlegene_collectall_BMDM2_WTMM_3stim_ISnorm_500genes.txt")
genes.mi = reshape2::dcast(genes.mi, gene~type, value.var = 'cc')
genes.mi$diff = genes.mi$BMDM2_WT - genes.mi$BMDM2_MM
clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
genes.mi$clusters = clusters$manual_fix[match((genes.mi$gene), clusters$gene)]
genes.mi$clusters = gsub("mB","mC", genes.mi$clusters)
genes.mi$clusters = ifelse(genes.mi$clusters=="mA","AP1",
                                 ifelse(genes.mi$clusters=="mC","NFkB",
                                        ifelse(genes.mi$clusters=="mD","NFkB&p38",
                                               ifelse(genes.mi$clusters=="mE","NFkB|IRF","IRF")))) 
genes.mi$clusters = factor(genes.mi$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))                                  
ggplot(genes.mi, aes(clusters, diff))+geom_boxplot(outlier.shape = NA)+
  theme_bw(base_size = 14)+ylab("max MI diff (WT-MM)")+ylim(-0.3,0.2)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
genes.mi$rank = rank(genes.mi$diff)
genes.mi$label = ifelse(abs(genes.mi$diff)>0.1, genes.mi$gene,"")
genes.mi$label = ifelse(grepl("Saa3|Tnfaip3$|Nfkbie|Vcam1|Tlr2|Cd40|Clec4e|Il1a|Il1b|Ptgs2|Sod2|Mt2|Slc28a2", genes.mi$gene), genes.mi$gene,"")
ggplot(genes.mi[grepl("NFkB$",genes.mi$clusters),], aes(rank, diff))+geom_bar(stat = "identity")+
  theme_bw(base_size = 14)+ylab("max MI diff (WT-MM)")+ylim(-0.3,0.2)+geom_hline(yintercept = -0.1,linetype="dashed"  )+geom_hline(yintercept = 0.1,linetype="dashed"  )+
  geom_hline(yintercept = -0.2,linetype="dashed"  )+geom_hline(yintercept = 0.2,linetype="dashed"  )+
  geom_text_repel(mapping = aes(label = label), xlim = c(50,450), box.padding = 0.3) 

  

nfkb.mi = rbind(read.delim("./infotheo/SLEMI_3stim_NFkBsignalingcodons_IkBaWT_zscore.txt"),
                read.delim("./infotheo/SLEMI_3stim_NFkBsignalingcodons_IkBaMM_zscore.txt"))
nfkb.mi = dcast(nfkb.mi, codon~time, value.var = 'cc')
nfkb.mi$diff = nfkb.mi$WT - nfkb.mi$MM
ggplot(nfkb.mi, aes(codon, diff))+geom_bar(stat="identity")+theme_bw(base_size = 14)+ylab("max MI diff (WT-MM)")


#makes supp tables-----
clusters = readxl::read_excel("F://scRNAseq_macro/SuppTables/TableS3_genes2GRS.xlsx")
clusters$clusters = gsub("mB","mC", clusters$manual_fix)
clusters$clusters = ifelse(clusters$clusters=="mA","AP1",
                           ifelse(clusters$clusters=="mC","NFkB",
                                  ifelse(clusters$clusters=="mD","NFkB&p38",
                                         ifelse(clusters$clusters=="mE","NFkB|IRF","IRF")))) 
clusters$clusters = factor(clusters$clusters, level = c("AP1","NFkB","NFkB&p38","NFkB|IRF","IRF"))   
write.table(clusters, "F://scRNAseq_macro/SuppTables/TableS3_genes2GRS.txt", quote=F,sep="\t",row.names = F)


###############################
# 3D plots of data----
macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = subset(macro, subset = timept == "3hr")
data = GetAssayData(object = macro, assay = "ISnorm", slot = "data")
data <- as.data.frame( as.matrix(data))
colnames(data)= macro$stimulus[match(colnames(data), colnames(macro))]
data.t= data.frame(t(data))
data.t=cbind(stimulus=rownames(data.t),data.t)
data.t$stimulus = gsub("\\..*","",data.t$stimulus)
ggplot(data.t, aes(Icam1, stimulus))+geom_violin(outlier.shape = NA)+
  geom_point(aes(color = stimulus), alpha = 0.2,position = "jitter")+theme_bw(base_size = 14)
ggplot(data.t, aes(Cmpk2, Nfkbiz))+geom_point(aes(color = stimulus),alpha = 0.5)+theme_bw(base_size = 14)
ggplot(data.t, aes(Cmpk2, Icam1))+geom_point(aes(color = stimulus),alpha = 0.5)+theme_bw(base_size = 14)
ggplot(data.t, aes(Nfkbiz, Icam1))+geom_point(aes(color = stimulus),alpha = 0.5)+theme_bw(base_size = 14)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")

## An empty plot with 3 axes
qplot(x=0, y=0, z=0, geom="blank") + 
  theme_void() +
  axes_3D()
## Axes can be populated with points using the function stat_3D.
ggplot(data.t, aes(x=Cmpk2, y=Nfkbiz, z=Icam1, color=stimulus)) +
# ggplot(data.t, aes(x=Ccl5, y=Cxcl10, z=Tnf, color=stimulus)) +   
  theme_void()+
  axes_3D() +
  stat_3D()+
  labs_3D(labs=c("Cmpk2", "Nfkbiz", "Icam1"),
  # labs_3D(labs=c("Ccl5", "Cxcl10", "Tnf"),
          angle=c(45,-45,90),hjust=c(0,1,2),vjust=c(2,2,2))


library(rgl)
colors_list = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")
with(data.t, scatter3D(x=Cmpk2, y=Nfkbiz, z=Icam1, colvar = as.integer(as.factor(stimulus)), 
                       col = colors_list, bty = "b2",pch = 19, cex = 0.6, alpha = 0.5,  
                       xlab = "Cmpk2", ylab = "Nfkbiz", zlab = "Icam1"))
