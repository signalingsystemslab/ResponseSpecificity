# Response Specificity code
# Contact: katherinesheu[at]ucla.edu
# last updated Oct 2021

# anxillary functions----
PCA_from_file = function (file, center = TRUE, scale = FALSE, fread = FALSE, 
          tol = sqrt(.Machine$double.eps)) {
  require(data.table)
  if (fread == T) {
    data = fread(file)
    data = data[rowSums((data[, -1, with = F] == 0)) < ncol(data[-1]), 
    ]
    t.data = t(data[, -1, with = F])
    pca <- prcomp(t.data, scale = scale, center = center, 
                  tol = tol)
    pca_scores = pca$x
    pca_scores = cbind(Score = gsub("-", ".", rownames(pca_scores)), 
                       pca_scores)
    pca_loadings = pca$rotation
    pca_loadings = cbind(Loading = data[, 1, with = F], 
                         pca_loadings)
    pca_evalues = pca$sdev
  }
  else {
    data = read.delim(file, header = T, stringsAsFactors = F)
    data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), 
    ]
    t.data = t(data[, -1])
    pca <- prcomp(t.data, scale = scale, center = center)
    pca_scores = pca$x
    pca_scores = cbind(Score = rownames(pca_scores), pca_scores)
    pca_loadings = pca$rotation
    pca_loadings = cbind(Loading = data[, 1], pca_loadings)
    pca_evalues = pca$sdev
  }
  name = sub(".txt", "", file)
  savename = paste(name, "_prcomp_scores.txt", sep = "")
  write.table(pca_scores, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  savename = paste(name, "_prcomp_loadings.txt", sep = "")
  write.table(pca_loadings, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  savename = paste(name, "_prcomp_sdev.txt", sep = "")
  write.table(pca_evalues, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  print(summary(pca))
  screeplot(pca)
}

plot_pca=function (file, info.name, info.type, title = "", labels = TRUE, 
                   PCx = "PC1", PCy = "PC2", pt.size=1, ellipse = F, conf = 0.95, density = F, 
                   fliph = F, flipv = F) {
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

intersect_doPCA_from_file_and_project_second_dataset = function (file, file2, train_string, center = TRUE, scale = FALSE, 
          fread = F) {
  if (fread == T) {
    data1 = fread(file)
    data1 = data[rowSums((data1[, -1, with = F] == 0)) < 
                   ncol(data1[-1]), ]
    data2 = fread(file2)
    data1 = data1[!duplicated(data1[, 1, with = F]), ]
    data2 = data2[!duplicated(data2[, 1, with = F]), ]
    common.genes = intersect_all(data[, 1, with = F], data2[, 
                                                            1, with = F])
    data1 = data1[data1[, 1, with = F] %in% common.genes, 
    ]
    data2 = data2[data2[, 1, with = F] %in% common.genes, 
    ]
    data1 = data1[order(data1[, 1, with = F]), ]
    data2 = data2[order(data2[, 1, with = F]), ]
    rownames(data1) = make.names(data1[, 1], unique = TRUE)
    t.data = data.frame(t(data1[, -1, with = F]))
  }
  else {
    data1 = read.delim(file, header = T, stringsAsFactors = F)
    data2 = read.delim(file2, header = T, stringsAsFactors = F)
    data1 = data1[!duplicated(data1[, 1]), ]
    data2 = data2[!duplicated(data2[, 1]), ]
    data1 = data1[rowSums((data1[, -1] == 0)) < ncol(data1[-1]), 
    ]
    common.genes <- intersect((data1[, 1]), (data2[, 1]))
    data <- data1[(data1[, 1]) %in% common.genes, ]
    data2 <- data2[(data2[, 1]) %in% common.genes, ]
    data = data[order(data[, 1]), ]
    data2 = data2[order(data2[, 1]), ]
    t.data = t(data[, -1])
  }
  pca <- prcomp(t.data, scale = scale, center = center)
  pca_scores = pca$x
  pca_scores = cbind(Score = rownames(pca_scores), pca_scores)
  pca_loadings = pca$rotation
  pca_loadings = cbind(Loading = data[, 1], pca_loadings)
  pca_evalues = pca$sdev
  pca_scale = pca$scale
  name = sub(".txt", "", file)
  savename = paste(name, "_prcomp_scores.txt", sep = "")
  write.table(pca_scores, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  savename = paste(name, "_prcomp_loadings.txt", sep = "")
  write.table(pca_loadings, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  savename = paste(name, "_prcomp_sdev.txt", sep = "")
  write.table(pca_evalues, savename, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  print(summary(pca))
  screeplot(pca)
  t.data2 = t(data2[, -1])
  rotated.data2 = scale(t.data2, pca$center, pca$scale) %*% 
    pca$rotation
  rotated.data2 = cbind(Sample = rownames(rotated.data2), 
                        rotated.data2)
  name2 = sub(".txt", "", file2)
  savename_intermed = paste(name2, train_string, sep = "_")
  savename2 = paste(savename_intermed, "_prcomp_rotated.txt", 
                    sep = "")
  write.table(rotated.data2, savename2, sep = "\t", row.names = FALSE, 
              quote = FALSE)
  rotated.data2
}

plot_pca_projection=function (file, rotated.file, info.name, info.type, info.name2, 
                              info.type2, title = "Projection", labels = F, PCx = "PC1", 
                              PCy = "PC2", pt.size=1, ellipse = F, conf = 0.95, fliph = F, flipv = F, 
                              save = F, savename = "") {
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

# load libraries----
library(edgeR);library(mclust);library(pheatmap);library(e1071);library(ksheu.library1);library(ggplot2);library(svglite)
library(RColorBrewer);library(ConsensusClusterPlus);library(Matrix);library(matrixStats);library(openxlsx);
library(rTensor);library(reshape2);library(rgl);library(ggrepel);library(plot3D);library(plot3Drgl);library(gridExtra)
library(patchwork);library(plot3D);library(plot3Drgl);library(pheatmap);library(grid);library(umap);library(ggthemes)
library(ggpubr);library(ggridges);library(reshape2);library(Seurat);library(dplyr);library(ksheu.library1)
library(SLEMI);library(ggplot2);library(ggrepel);library(ksheu.library1);library(Seurat);library(pheatmap)


###################################################### Figure1----
# Figure 1b----
files = list.files(".output/motif_results/", "knownResults.txt", recursive = T)
tmp = read.delim(paste0("F://scRNAseq_macro/bulk_rnaseq/motif_results/",files[1])) #files[1]
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
frame.truncate.500 = aggregate(frame.truncate[,c(2,3)], by = list(motif = frame.truncate$type), FUN = mean)
ggplot(frame.truncate.500, aes(motif, Log.P.value*-1, fill = motif))+geom_bar(stat = "identity", position = position_dodge())+theme_bw(base_size = 16)

tmp = read.delim(paste0("F://scRNAseq_macro/bulk_rnaseq/motif_results/",files[2])) #files[1]
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
frame.truncate.all = aggregate(frame.truncate[,c(2,3)], by = list(motif = frame.truncate$type), FUN = mean)
ggplot(frame.truncate.all, aes(motif, Log.P.value*-1, fill = motif))+geom_bar(stat = "identity", position = position_dodge())+theme_bw(base_size = 16)

frame.truncate = rbind(data.frame(frame.truncate.500, type = "selected"),
                       data.frame(frame.truncate.all, type = "induced"))
ggplot(frame.truncate[!grepl("other", frame.truncate$motif),], aes(type.1, log(Log.P.value*-1), fill = motif))+
  geom_bar(stat = "identity", position = position_dodge())+theme_bw(base_size = 16)+
  ylab("ln(-ln(p-value))")

# Figure 1c----
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(DOSE)

#my genelist
gene <- rownames(mat.scgenes)
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
gene.df$cluster = tmp$kmcluster.cluster[match(gene.df$SYMBOL, tmp$gene)]

induced = read.delim("output/Cheng2017_induced.txt")
induced.df <- bitr(induced$gene, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db)
induced.df$selected = ifelse(induced.df$SYMBOL %in% sc.inducedgenes, 0, 1)

# compare induced vs selected
ego =  clusterProfiler::compareCluster(geneClusters = ENTREZID~selected, 
                                       data=induced.df,
                                       fun = "enrichGO", 
                                       OrgDb = org.Mm.eg.db, 
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
p1 <- dotplot(ego2, showCategory=3) + ggtitle("")
p1

# Figure 1d----
mat = read.delim("./output/Cheng2017_induced.txt") #1502 induced genes
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

# Figure 1d
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
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
colors_list = list(stimulus = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
p= pheatmap(data[sc.order,col_order], scale = "row", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D2", 
            colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103),
            gaps_row = rowseps[-5],gaps_col = colseps[-7],
            breaks=c(-2,seq(-1,1,length=100),2),
            annotation_col = data.frame(macro@meta.data[,c(6,7)]), 
            annotation_colors = colors_list,
            show_colnames = F, show_rownames = F)
ggsave(p, filename = "heatmap_singlecell.ISnorm_M0_rep2only_3hr.png")

# Figure 1e----
# plotted Ccl5, Tnf, Cxcl10 Ridgeplots
macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
colors_list = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")
p1=RidgePlot(object = subset(macro, subset= timept=="3hr"|timept=="0.0hr"), sort = "increasing",
             features = c("Ccl5"), group.by = "stimulus", assay = "ISnorm" , y.max = ymax) +
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

# Figure 1f----
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
colors_list = c(Unstim = "gray", CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")

data = macro[["ISnorm"]]@data
data = data.frame(data)
meta = macro@meta.data
colnames(data) = paste0(meta$orig.ident, "_",meta$type, "_",meta$stimulus, meta$timept)
colnames(data) = paste0(meta$type, "_",meta$stimulus, meta$timept)
data.t = data.frame(gene = colnames(data), t(data))

data.aggregate = aggregate(data.t[,-1], by = list((data.t$gene)), FUN = function(x){(var(x))}) #var
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


########################## Figure 2
# Figure 2a----
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

PCA_from_file("./output/macrophage_M0_rep2only_500genes_DBEC_3hr.txt", center =T, scale = F)
plot_pca("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt", 
         samptag.all$Cell_Index, as.factor(samptag.all$stimulus), PCx="PC1", PCy="PC2",pt.size=0.1,labels = F, ellipse =T)

# Figure 2b----
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
pca.macro = read.delim("./output/macrophage_M0_rep2only_500genes_DBEC_3hr_prcomp_scores.txt")
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


#umap
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


# Figure 2cd----
for (i in c("3hr")){
  
  print(i)
  macro = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
  # macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
  
  filename = "M0all.Dec2020"
  genesetname = "500genes" 
  macro = subset(macro, subset= timept==i)
  # macro = subset(macro, subset= stimulus!="CpG")
  
  # data = macro[["RNA"]]@data
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  # my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
  
  #--------------------------------------ML using CARET-----------------------------------------
  library(caret)
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
  training <- my.dataframe[ inTraining,]
  testing  <- my.dataframe[-inTraining,]
  
  if (0) {#plsr --------------------------------------------------------------------------------
    fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
    fit_pls <- train(label ~ .,  data = training,
                     method = "pls",
                     ## Center and scale the predictors for the training
                     ## set and all future samples.
                     # preProc = c("center", "scale"),
                     tuneLength = 20,
                     trControl = fitControl,
                     metric = "ROC")
    fit_pls
    ggplot(fit_pls)
    plsClasses <- predict(fit_pls, newdata = testing)
    plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
    head(plsClasses)
    head(plsProbs)
    confusion = confusionMatrix(data = plsClasses, as.factor(testing$label))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T, fontsize_number = 14,main = paste0("pls_",filename,i)) 
    # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
    # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
    
    # varImp = varImp(fit_pls); varImp = varImp$importance
    
    # save the model to disk
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".rds"))
  }
  
  if (1){ #random forest------------------------------------------------------------------------
    print("running rf")
    fitControl <- trainControl(method='repeatedcv', number=10, repeats=3)
    metric <- "Accuracy" #Metric compare model is Accuracy
    set.seed(1)
    mtry <- sqrt(ncol(training) -1) #Number randomely variable selected is mtry
    tunegrid <- expand.grid(.mtry=mtry)
    fit_rf_default <- train(label~., 
                            data=training, 
                            method="rf", 
                            metric='Accuracy', 
                            tuneGrid=tunegrid, 
                            trControl=fitControl)
    
    print(fit_rf_default)
    rfClasses <- predict(fit_rf_default, newdata = testing)
    confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("rf_",filename,i))
    rf_prob=predict(fit_rf_default, newdata = testing, type = "prob")
    # varImp = varImp(fit_rf_default); 
    # varImp = data.frame(varImp$importance)
    
    
    # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
    #   geom_point() +
    #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
    #   scale_color_discrete(name="Variable Group") +
    #   ylab("IncNodePurity") +
    #   xlab("Variable Name") +
    #   coord_flip()
    # plot(varImp, top = 20)
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_rf_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(fit_rf_default, paste0("./analysis_rhapsody_500genes/MLfit_rf_",filename,"_",genesetname,"_",i, ".rds")) #on poseidon then moved back
  }
  
  if (0){ # SVM -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Linear <- train(label ~., data = training, 
                        method = "svmLinear",
                        trControl=trctrl,
                        # preProcess = c("center", "scale"),
                        tuneLength = 10)
    
    svm_Linear
    test_pred <- predict(svm_Linear, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("svmL_",filename,i)) 
    # varImp = varImp(svm_Linear); #varImp = varImp$importance
    # plot(varImp, top = 20)
    ggsave(p,filename=paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for memory problems
    
  }
  
  if (0){ # SVM radial -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Radial <- train(label ~., data = training,
                        method = "svmRadial",
                        trControl = trctrl,
                        # preProcess = c("center","scale"),
                        tuneLength = 10)
    
    # Print the best tuning parameter sigma and C that maximizes model accuracy
    svm_Radial$bestTune
    
    svm_Radial
    test_pred <- predict(svm_Radial, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("svmR_",filename,i)) 
    # varImp = varImp(svm_Radial); #varImp = varImp$importance
    # plot(varImp, top = 20)
    # saveRDS(svm_Radial, paste0("./analysis_rhapsody/MLfit_svmRadial_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for memory problems
    
  }
  
  if (0){ # kNN-------------------------------------------------------------------
    # my.dataframe = data.frame(cbind(label = colnames(data), macro@reductions$pca@cell.embeddings))
    # set.seed(1)
    # inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
    # training <- my.dataframe[ inTraining,]
    # testing  <- my.dataframe[-inTraining,]
    ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
    knnFit <- train(label ~ ., 
                    data = training, 
                    method = "knn", 
                    trControl = ctrl, 
                    # preProcess = c("center","scale"), 
                    tuneLength = 20)
    
    #Output of kNN fit
    knnFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(knnFit)
    knnPredict <- predict(knnFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(knnPredict, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("knn_",filename,i)) 
    # varImp = varImp(knnFit); varImp = varImp$importance
    # saveRDS(knnFit, paste0("./analysis_rhapsody/MLfit_knn_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(knnFit, paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for memory problems
    
  }
  
  if (0){ #naive bayes
    library(e1071)
    ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
    nbFit <- train(label ~ ., 
                   data = training, 
                   method = "nb", 
                   trControl = ctrl, 
                   # preProcess = c("center","scale"),
                   tuneLength = 20)
    
    #Output of nb fit
    nbFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(nbFit)
    nbPredict <- predict(nbFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(nbPredict, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T, fontsize_number = 14,main = paste0("nb_",filename,i)) 
    # varImp = varImp(nbFit); 
    # varImp = varImp$importance
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(nbFit, paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for memory problems
    
  }
}

# Load previously trained and saved machine learning models 
list = list.files(pattern= "MLfit_rf_M0all")
# for (i in c("knn", "pls", "rf", "svmLinear", "svmRadial") ){
for (i in c("rf") ){
  for (j in c("3hr") ){
    tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M0all.Dec2020_500genes_", j, ".rds"))
       assign(paste0(i, "_", j), tmp)
  }
}

# prediction probabilities
for (i in c("3hr") ){
  
  tmp = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))

  tmp = subset(tmp, subset=timept==i)
  tmp.meta = tmp@meta.data
  tmp2 = tmp[["ISnorm"]]@data
  # tmp2 = tmp[["RNA"]]@data
  tmp2 = data.frame(tmp2)
  colnames(tmp2) = paste0(tmp.meta$stimulus, "_",tmp.meta$timept)
  # rownames(tmp2) = gsub("-ENSMUST..*","", tmp[["SCT"]]@counts@Dimnames[[1]])
  # rownames(tmp2) = gsub("..NM..*","", rownames(tmp2))
  
  my.dataframe = cbind(label = colnames(tmp2), data.frame(t(tmp2)))
  # my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
  
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
  training <- my.dataframe[ inTraining,]
  assign(paste0("testing", "_",i), my.dataframe[-inTraining,])
  testing = my.dataframe[-inTraining,]
  
  
  fit_rf_default = get(paste0("rf_",i))
  print(fit_rf_default)
  rfClasses <- predict(fit_rf_default, newdata = testing)
  confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
  confusion.table = (confusion$table)
  confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
  pheatmap(confusion.table, cluster_rows = F, cluster_cols = F,
           # colorRampPalette(c("lightblue", "white", "red"))(50),
           # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30),
           breaks =  seq(0, 1, by = .01),
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
           display_numbers = T, fontsize_number = 16)
  
  predProbs <- extractProb(list(fit_rf_default),
                           testX = testing[,-1], testY = testing$label)
  predProbs = predProbs[order(predProbs$obs, decreasing = T),]
  colors_list = list(obs = c(Unstim = "gray", CpG_0.25hr="#F8766D", IFNb_0.25hr="#B79F00",LPS_0.25hr= "#00BA38",P3CSK_0.25hr= "#00BFC4",PIC_0.25hr= "#619CFF",TNF_0.25hr= "#F564E3",
                             CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3",
                             CpG_3hr="#F8766D", IFNb_3hr="#B79F00",LPS_3hr= "#00BA38",P3CSK_3hr= "#00BFC4",PIC_3hr= "#619CFF",TNF_3hr= "#F564E3",
                             CpG_8hr="#F8766D", IFNb_8hr="#B79F00",LPS_8hr= "#00BA38",P3CSK_8hr= "#00BFC4",PIC_8hr= "#619CFF",TNF_8hr= "#F564E3"),
                     pred = c(Unstim = "gray", CpG_0.25hr="#F8766D", IFNb_0.25hr="#B79F00",LPS_0.25hr= "#00BA38",P3CSK_0.25hr= "#00BFC4",PIC_0.25hr= "#619CFF",TNF_0.25hr= "#F564E3",
                              CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3",
                              CpG_3hr="#F8766D", IFNb_3hr="#B79F00",LPS_3hr= "#00BA38",P3CSK_3hr= "#00BFC4",PIC_3hr= "#619CFF",TNF_3hr= "#F564E3",
                              CpG_8hr="#F8766D", IFNb_8hr="#B79F00",LPS_8hr= "#00BA38",P3CSK_8hr= "#00BFC4",PIC_8hr= "#619CFF",TNF_8hr= "#F564E3"))
  pheatmap(predProbs[predProbs$dataType == "Test",c(1:6)], scale = "none",cluster_rows = F, cluster_cols = F,
           annotation_row = predProbs[predProbs$dataType == "Test",c(7:8)],
           show_rownames = F, main = "prediction probabilities",
           annotation_colors = colors_list)
  plotClassProbs(predProbs)
  plotClassProbs(predProbs[predProbs$dataType == "Test",])
  
  
}


# Figure 2e----
table = data.frame()
for (i in c("3hr") ){
  print(i)
  fit_rf_default = get(paste0("rf_",i) )
  testing = get( paste0("testing_",i) )
  
  print(fit_rf_default)
  rfClasses <- predict(fit_rf_default, newdata = testing)
  confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
  confusion.table = (confusion$table)
  values = as.data.frame(confusion$byClass)
  
  table = rbind(table, values)
  
}
#for testing M0
table$time = c(
  # rep(0.25, 6), 
  # rep(0.5, 5), 
  # rep(1, 6), 
  rep(3, 6) 
  # rep(5, 5), 
  # rep(8, 6)
)
table$stim = c(  rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 1))
table[is.na(table)] <-0
table$FPR = 1-table$Specificity
table$FDR = 1-table$Precision

ggplot(table[grepl("3", table$time),], aes(FDR, F1))+ geom_point(aes(color = stim), size=5)+
  theme_bw(base_size = 16)
ggplot(table[grepl("3", table$time),], aes(stim, F1,fill=stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

table.m = melt(table)
ggplot(table.m[grepl("FPR|FDR", table.m$variable),], aes(stim, value,fill=stim))+ geom_bar(position="dodge", stat="identity")+
  facet_grid(~variable, scales = "free_y")+ylim(0,0.25)+
  theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Figure 2f----
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

################################# Figure 3----
# run mutual_info3_pairwise_SLEMI.R script (pasted below)
# mutual information SLEMI----
# last mod. 5/17/21, ksheu

install_github("sysbiosig/SLEMI")
library(SLEMI)

# channel capacity for each pairwise stimuli at each timepoint----
collect = data.frame()
list = c("CpG","P3CSK","LPS", "TNF", "IFNb",  "PIC") #

for (i in c("3hr")){
  print(i)
  
  macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  macro = subset(macro, subset = timept==i)
  
  
  for (j in seq(1:(length(list)-1) )){
    for (k in seq(j+1,length(list)) ){
      
      for (gene in c(rownames(macro))){
        
        skip_to_next <- FALSE
        print(gene)
        
        tryCatch(
          {
            
            
            print(list[j])
            print(list[k])
            
            macro.subset = subset(x = macro, subset = stimulus==list[j]|stimulus == list[k])
            
            data = macro.subset[["ISnorm"]]@data
            data = data.frame(data)
            meta = macro.subset@meta.data
            colnames(data) = paste0(meta$stimulus, "_",meta$timept)
            my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
            my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
            rownames(my.dataframe) = seq(1:nrow(my.dataframe))
            my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("label", gene)]
            str(my.dataframe.subset)
            
            #--------------------------------mi using SLEMI -----------------------
            
            output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "label", 
                                                    response = colnames(my.dataframe.subset)[-1],
                                                    testing=T, boot_prob = 0.5, boot_num = 30, testing_cores = 4)
            sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
            
            
            tmp = data.frame(time = i, stim1 = list[j], stim2 = list[k], cc = output_capacity$cc, sd = sd, gene = gene)
            collect = rbind(collect, tmp)
            
          }, error = function(e) { skip_to_next <<- TRUE})
        
        if(skip_to_next) { next }  
        
        
      }
    }
  }
}
view(collect)
write.table(collect, "./infotheo/channel_capacity_pairwise_M0rep2only_ISnorm_allgenes.txt",row.names = F, sep = "\t",quote = F)

# Figure 3b ----
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

# Figure 3c----
collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
ggplot(collect_all[grepl("3hr", collect_all$time),], aes( cc))+ geom_density()+ #graph it USE THIS
  # geom_point(aes(color = time),position = position_jitter(seed = 1))+
  ggtitle("Single gene output")+theme_bw(base_size = 20)+
  geom_vline(xintercept = 0.7, linetype="dotted", size = 1)+
  geom_vline(xintercept = (log10(6)/log10(2)), linetype="dashed", size = 1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0,2.7)) + scale_y_continuous(expand = c(0, 0))


dcast = dcast(collect_all, gene~time, value.var = "cc")
rownames(dcast) = dcast$gene
p=pheatmap(na.omit(dcast[(dcast[,4] > .7),4, drop=F]), scale = "none", clustering_method = "ward.D2", 
           colorRampPalette((brewer.pal(n = 11, name = "Greens")))(103),
           cluster_cols = F, cluster_rows = T, show_rownames = T)
row_order = as.character(na.omit(dcast$gene[(dcast[,4] > .7)])[p$tree_row$order])

# Figure 3e----
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


############################### Figure 4
# Figure 4cde----
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


# Figure 4f----
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

features.wanted = readxl::read_excel("./infotheo/data_table_KS.xlsx", sheet = 3)

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
features.wanted = readxl::read_excel("./infotheo/data_table_KS.xlsx", sheet = 2) #WT zscore
features.wanted = readxl::read_excel("./infotheo/data_table_KS.xlsx", sheet = 4) #MM zscore

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
# clusters = readxl::read_excel("F://scRNAseq_macro/scRNAseq_macro/fit_mean_modv3_p38input/reassign_genes_to_grns_unweightedcost_500genes_Jan2021_ks.xlsx")
clusters = readxl::read_excel("./SuppTables/TableS3_genes2GRS.xlsx")
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


##################################Figure 5-----
# Figure 5b----
sdev = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_sdev.txt")
sdev$var = unlist(sdev^2)
sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
sum(sdev$pve[1:20])

#UMAP of unscaled PCA scores M0M1M2------------
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

PCA_from_file("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr.txt")

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

macro = readRDS("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr.rds")

if(1){
  library(scales)
  time="3hr"
  red = "pca"
  wanted = rownames(macro@meta.data)[macro@meta.data$timept==time&macro@meta.data$type=="M0"]
  list = as.list(wanted);  names(list) = macro@meta.data$stimulus[match(list, rownames(macro@meta.data))]
  list = split(unlist(list, use.names = FALSE), rep(names(list), lengths(list))) #merge under unique names
  p1=DimPlot(object = macro, cells.highlight = list, reduction = red, group.by = "stimulus", raster = F,sizes.highlight=0.1,cols.highlight = c(rev(hue_pal() (length( unique(names(list)) )) ) ),pt.size = 0.1)+theme(legend.position = "None")
  wanted = rownames(macro@meta.data)[macro@meta.data$timept==time&macro@meta.data$type=="M1_IFNg"]
  list = as.list(wanted);  names(list) = macro@meta.data$stimulus[match(list, rownames(macro@meta.data))]
  list = split(unlist(list, use.names = FALSE), rep(names(list), lengths(list))) #merge under unique names
  p2=DimPlot(object = macro, cells.highlight = list, reduction = red, group.by = "stimulus", raster = F,sizes.highlight=0.1,cols.highlight = c(rev(hue_pal() (length( unique(names(list)) )) ) ),pt.size = 0.1)+theme(legend.position = "None")
  wanted = rownames(macro@meta.data)[macro@meta.data$timept==time&macro@meta.data$type=="M2_IL4"]
  list = as.list(wanted);  names(list) = macro@meta.data$stimulus[match(list, rownames(macro@meta.data))]
  list = split(unlist(list, use.names = FALSE), rep(names(list), lengths(list))) #merge under unique names
  p3=DimPlot(object = macro, cells.highlight = list, reduction = red, group.by = "stimulus", raster = F,sizes.highlight=0.1,cols.highlight = c(rev(hue_pal() (length( unique(names(list)) )) ) ),pt.size = 0.1)+theme(legend.position = "None")
  p1|p2|p3
}

# Figure 5cd (distance metric M0, M1, M2) ---- 
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


# Figure 5e----
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

# Figure 5fg----
# Load previously trained M0 machine learning models 
list = list.files(pattern= "MLfit_rf_M0all")
# for (i in c("knn", "pls", "rf", "svmLinear", "svmRadial") ){
for (i in c("rf") ){
  for (j in c("3hr") ){
    tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M0all.Dec2020_500genes_", j, ".rds"))
    assign(paste0(i, "_", j), tmp)
  }
}
#load testing data for M1 M2, and test on M0 models----------
for (i in c("3hr") ){
  
  # tmp = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))   # for testing M1
  tmp = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")) # for testing M2
  
  
  # tmp = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
  
  tmp = subset(tmp, subset=timept==i)
  tmp.meta = tmp@meta.data
  tmp2 = tmp[["ISnorm"]]@data
  # tmp2 = tmp[["RNA"]]@data
  tmp2 = data.frame(tmp2)
  colnames(tmp2) = paste0(tmp.meta$stimulus, "_",tmp.meta$timept)
  # rownames(tmp2) = gsub("-ENSMUST..*","", tmp[["SCT"]]@counts@Dimnames[[1]])
  # rownames(tmp2) = gsub("..NM..*","", rownames(tmp2))
  
  my.dataframe = cbind(label = colnames(tmp2), data.frame(t(tmp2)))
  # my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
  
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
  training <- my.dataframe[ inTraining,]
  assign(paste0("testing", "_",i), my.dataframe[-inTraining,])
  testing = my.dataframe[-inTraining,]
  

  fit_rf_default = get(paste0("rf_",i))
  print(fit_rf_default)
  rfClasses <- predict(fit_rf_default, newdata = testing)
  confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
  confusion.table = (confusion$table)
  confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
  pheatmap(confusion.table, cluster_rows = F, cluster_cols = F,
           # colorRampPalette(c("lightblue", "white", "red"))(50),
           # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30),
           breaks =  seq(0, 1, by = .01),
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
           display_numbers = T, fontsize_number = 16)
  
  predProbs <- extractProb(list(fit_rf_default),
                           testX = testing[,-1], testY = testing$label)
  predProbs = predProbs[order(predProbs$obs, decreasing = T),]
  colors_list = list(obs = c(Unstim = "gray", CpG_0.25hr="#F8766D", IFNb_0.25hr="#B79F00",LPS_0.25hr= "#00BA38",P3CSK_0.25hr= "#00BFC4",PIC_0.25hr= "#619CFF",TNF_0.25hr= "#F564E3",
                             CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3",
                             CpG_3hr="#F8766D", IFNb_3hr="#B79F00",LPS_3hr= "#00BA38",P3CSK_3hr= "#00BFC4",PIC_3hr= "#619CFF",TNF_3hr= "#F564E3",
                             CpG_8hr="#F8766D", IFNb_8hr="#B79F00",LPS_8hr= "#00BA38",P3CSK_8hr= "#00BFC4",PIC_8hr= "#619CFF",TNF_8hr= "#F564E3"),
                     pred = c(Unstim = "gray", CpG_0.25hr="#F8766D", IFNb_0.25hr="#B79F00",LPS_0.25hr= "#00BA38",P3CSK_0.25hr= "#00BFC4",PIC_0.25hr= "#619CFF",TNF_0.25hr= "#F564E3",
                              CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3",
                              CpG_3hr="#F8766D", IFNb_3hr="#B79F00",LPS_3hr= "#00BA38",P3CSK_3hr= "#00BFC4",PIC_3hr= "#619CFF",TNF_3hr= "#F564E3",
                              CpG_8hr="#F8766D", IFNb_8hr="#B79F00",LPS_8hr= "#00BA38",P3CSK_8hr= "#00BFC4",PIC_8hr= "#619CFF",TNF_8hr= "#F564E3"))
  pheatmap(predProbs[predProbs$dataType == "Test",c(1:6)], scale = "none",cluster_rows = F, cluster_cols = F,
           annotation_row = predProbs[predProbs$dataType == "Test",c(7:8)],
           show_rownames = F, main = "prediction probabilities",
           annotation_colors = colors_list)
  plotClassProbs(predProbs)
  plotClassProbs(predProbs[predProbs$dataType == "Test",])
  
  
}


#for testing M1/M2 plot summary stats-----
table$time = c(
  rep(3, 6))
table$stim = c(rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 1))
table$FPR = 1-table$Specificity
table$FDR = 1-table$Precision
table[is.na(table)] <-0

ggplot(table[grepl("3", table$time),], aes(FPR, F1))+ geom_point(aes(color = stim), size=5)+
  theme_bw(base_size = 16)

ggplot(table[grepl("3", table$time),], aes(as.factor(time), FPR, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+ylim(c(0,.4))

ggplot(table, aes(as.factor(time), fdr, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+ylim(c(0,1))
ggplot(table, aes(as.factor(time), `Pos Pred Value`, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+ylim(c(0,1))


table.m = melt(table[grepl("3",table$time),])
ggplot(table.m[grepl("FPR|FDR", table.m$variable),], aes(stim, value, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+facet_grid(~variable)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(table.m[grepl("F1|Balanced", table.m$variable),], aes(stim, value, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+facet_grid(~variable)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Figure 5h
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

###################################Figure 6
# Figure 6a----
# collect by keep top 20 forward selection on individual gene list for each timepoint-----
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M0_rep2only_ISnorm_500genes.txt")
# collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M1_IFNg_ISnorm_500genes.txt")
# collect_all = read.delim("./infotheo/SLEMI_singlegene_collectall_M2_IL4_gt80_ISnorm_500genes.txt")
collect_all = collect_all[order(collect_all$cc, decreasing = T), ]

collect = data.frame()
collect_dimensionbest = data.frame()

for (i in c("3hr")){
  print(i)
  macro = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
  # macro = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
  macro = subset(macro, subset= timept ==i)
  
  collect_dimension = (collect_all[grepl(i, collect_all$time),][(1:20),]) #start 1D
  # collect_dimension = readRDS(paste0("./infotheo/old/collect_dimension",i,"_10.rds")) #start 10D
  # collect_dimension = readRDS(paste0("./infotheo/collect_dimension",i,"_20.rds")) #start 20D
  # collect_dimension = readRDS(paste0("./infotheo/collect_dimension",i,"_30.rds")) #start 30D
  
  data = macro[["ISnorm"]]@data
  # data = macro[["SCT"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe$label = as.numeric(factor(my.dataframe$label, levels = sort(unique(my.dataframe$label))))
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  
  for (d in c(1:25)){
    print(paste0("dimension: ",d))
    
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    collect_dimension = collect_dimension[!duplicated(collect_dimension[,c(1:3)]), ]
    genesets = collect_dimension$gene[c(1:20)]
    print(genesets)
    
    collect_dimension_seed = data.frame() #start over to collect new dim's set once got the top20
    # for (g in 1:length(genesets)){
    collect_dimension = foreach(g = 1:length(genesets), .combine = rbind, .packages=c("SLEMI")) %dopar% {
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
        collect_dimension_seed = rbind(collect_dimension_seed, tmp)
      }
      collect_dimension_seed #rbinds all to collect_dimension
    }
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    saveRDS(collect_dimension, paste0("./infotheo/collect_dimension",i,"_",d,"_ISnorm3.rds"))
    
    collect_dimensionbest = rbind(collect_dimensionbest,
                                  data.frame(collect_dimension[1,]))
    
  }
}

# saveRDS(collect_dimensionbest, "./infotheo/collect_dimensionbest_ISnorm2.rds")
# saveRDS(collect_dimensionbest, "./infotheo/collect_dimensionbest_M1_ISnorm2.rds")
# saveRDS(collect_dimensionbest, "./infotheo/collect_dimensionbest_M2_ISnorm2.rds")


ggplot(collect_dimensionbest, aes(dim+1, cc))+ylim(0, 2.7)+
  geom_line(aes(color = time, group = time),size =2)+theme_bw(base_size = 14)+
  geom_hline(yintercept = (log10(6)/log10(2)), linetype="dotted", size = 1)

# Figure 6b----
mutual = readRDS("./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
genes.top2 = unlist(mutual[28,4])
genes.top5 = unlist(mutual[31,4])
genes.top15 = unlist(mutual[41,4])
genes = genes.top15 #c(genes.top2, genes.top5, genes.top15)

genesetname = "MItop15"

for (i in c("3hr")){
  
  print(i)

    macro = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
 
  filename = "M0all.Dec2020"

    macro = subset(macro, subset= timept==i)
  # macro = subset(macro, subset= stimulus!="CpG")
  
  # data = macro[["RNA"]]@data
  data = macro[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  # my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
  
  #--------------------------------------ML using CARET-----------------------------------------
  library(caret)
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
  training <- my.dataframe[ inTraining,]
  testing  <- my.dataframe[-inTraining,]
  
  if (0) {#plsr --------------------------------------------------------------------------------
    fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
    fit_pls <- train(label ~ .,  data = training,
                     method = "pls",
                     ## Center and scale the predictors for the training
                     ## set and all future samples.
                     # preProc = c("center", "scale"),
                     tuneLength = 20,
                     trControl = fitControl,
                     metric = "ROC")
    fit_pls
    ggplot(fit_pls)
    plsClasses <- predict(fit_pls, newdata = testing)
    plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
    head(plsClasses)
    head(plsProbs)
    confusion = confusionMatrix(data = plsClasses, as.factor(testing$label))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T, fontsize_number = 14,main = paste0("pls_",filename,i)) 
    # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
    # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
    
    # varImp = varImp(fit_pls); varImp = varImp$importance
    
    # save the model to disk
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".rds"))
  }
  
  if (1){ #random forest------------------------------------------------------------------------
    print("running rf")
    fitControl <- trainControl(method='repeatedcv', number=10, repeats=3)
    metric <- "Accuracy" #Metric compare model is Accuracy
    set.seed(1)
    mtry <- sqrt(ncol(training) -1) #Number randomely variable selected is mtry
    tunegrid <- expand.grid(.mtry=mtry)
    fit_rf_default <- train(label~., 
                            data=training, 
                            method="rf", 
                            metric='Accuracy', 
                            tuneGrid=tunegrid, 
                            trControl=fitControl)
    
    print(fit_rf_default)
    rfClasses <- predict(fit_rf_default, newdata = testing)
    confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("rf_",filename,i))
    rf_prob=predict(fit_rf_default, newdata = testing, type = "prob")
    # varImp = varImp(fit_rf_default); 
    # varImp = data.frame(varImp$importance)
    
    
    # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
    #   geom_point() +
    #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
    #   scale_color_discrete(name="Variable Group") +
    #   ylab("IncNodePurity") +
    #   xlab("Variable Name") +
    #   coord_flip()
    # plot(varImp, top = 20)
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_rf_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(fit_rf_default, paste0("./analysis_rhapsody_500genes/MLfit_rf_",filename,"_",genesetname,"_",i, ".rds")) #on poseidon then moved back
  }
  
  if (0){ # SVM -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Linear <- train(label ~., data = training, 
                        method = "svmLinear",
                        trControl=trctrl,
                        # preProcess = c("center", "scale"),
                        tuneLength = 10)
    
    svm_Linear
    test_pred <- predict(svm_Linear, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("svmL_",filename,i)) 
    # varImp = varImp(svm_Linear); #varImp = varImp$importance
    # plot(varImp, top = 20)
    ggsave(p,filename=paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ # SVM radial -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Radial <- train(label ~., data = training,
                        method = "svmRadial",
                        trControl = trctrl,
                        # preProcess = c("center","scale"),
                        tuneLength = 10)
    
    # Print the best tuning parameter sigma and C that maximizes model accuracy
    svm_Radial$bestTune
    
    svm_Radial
    test_pred <- predict(svm_Radial, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("svmR_",filename,i)) 
    # varImp = varImp(svm_Radial); #varImp = varImp$importance
    # plot(varImp, top = 20)
    # saveRDS(svm_Radial, paste0("./analysis_rhapsody/MLfit_svmRadial_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ # kNN-------------------------------------------------------------------
    # my.dataframe = data.frame(cbind(label = colnames(data), macro@reductions$pca@cell.embeddings))
    # set.seed(1)
    # inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
    # training <- my.dataframe[ inTraining,]
    # testing  <- my.dataframe[-inTraining,]
    ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
    knnFit <- train(label ~ ., 
                    data = training, 
                    method = "knn", 
                    trControl = ctrl, 
                    # preProcess = c("center","scale"), 
                    tuneLength = 20)
    
    #Output of kNN fit
    knnFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(knnFit)
    knnPredict <- predict(knnFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(knnPredict, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T,fontsize_number = 14,main = paste0("knn_",filename,i)) 
    # varImp = varImp(knnFit); varImp = varImp$importance
    # saveRDS(knnFit, paste0("./analysis_rhapsody/MLfit_knn_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(knnFit, paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ #naive bayes
    library(e1071)
    ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
    nbFit <- train(label ~ ., 
                   data = training, 
                   method = "nb", 
                   trControl = ctrl, 
                   # preProcess = c("center","scale"),
                   tuneLength = 20)
    
    #Output of nb fit
    nbFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(nbFit)
    nbPredict <- predict(nbFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(nbPredict, as.factor(testing$label) )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
               # colorRampPalette(c("lightblue", "white", "red"))(50),
               # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
               breaks =  seq(0, 1, by = .01),
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
               display_numbers = T, fontsize_number = 14,main = paste0("nb_",filename,i)) 
    # varImp = varImp(nbFit); 
    # varImp = varImp$importance
    ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".png"))
    saveRDS(nbFit, paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
    
  }
}

# Figure 6d----
# test random sets of 15 genes -----
store = data.frame()
set.seed(1)
for (i in c("3hr")){
  
  print(i)

    macro = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))

  for (r in c(1:50)){
    print(r)
    genes = sample(rownames(macro[["ISnorm"]]@data), 15)
    filename = "M0all.Dec2020"
    genesetname = paste0("15genes.rnd",r) #"500genes" 
    macro = subset(macro, subset= timept==i)
    # macro = subset(macro, subset= stimulus!="CpG")
    
    # data = macro[["RNA"]]@data
    data = macro[["ISnorm"]]@data
    data = data.frame(data)
    meta = macro@meta.data
    colnames(data) = paste0(meta$stimulus, "_",meta$timept)
    my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
    my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
    
    #--------------------------------------ML using CARET-----------------------------------------
    library(caret)
    set.seed(1)
    inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
    training <- my.dataframe[ inTraining,]
    testing  <- my.dataframe[-inTraining,]
    
    if (0) {#plsr --------------------------------------------------------------------------------
      fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
      fit_pls <- train(label ~ .,  data = training,
                       method = "pls",
                       ## Center and scale the predictors for the training
                       ## set and all future samples.
                       # preProc = c("center", "scale"),
                       tuneLength = 20,
                       trControl = fitControl,
                       metric = "ROC")
      fit_pls
      ggplot(fit_pls)
      plsClasses <- predict(fit_pls, newdata = testing)
      plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
      head(plsClasses)
      head(plsProbs)
      confusion = confusionMatrix(data = plsClasses, as.factor(testing$label))
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T, fontsize_number = 14,main = paste0("pls_",filename,i)) 
      # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
      # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
      
      # varImp = varImp(fit_pls); varImp = varImp$importance
      
      # save the model to disk
      ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes/MLfit_pls_",filename,"_",genesetname,"_",i, ".rds"))
    }
    
    if (1){ #random forest------------------------------------------------------------------------
      print("running rf")
      fitControl <- trainControl(method='repeatedcv', number=10, repeats=3)
      metric <- "Accuracy" #Metric compare model is Accuracy
      set.seed(1)
      mtry <- sqrt(ncol(training) -1) #Number randomely variable selected is mtry
      tunegrid <- expand.grid(.mtry=mtry)
      fit_rf_default <- train(label~., 
                              data=training, 
                              method="rf", 
                              metric='Accuracy', 
                              tuneGrid=tunegrid, 
                              trControl=fitControl)
      
      print(fit_rf_default)
      rfClasses <- predict(fit_rf_default, newdata = testing)
      confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      values = as.data.frame(confusion$byClass)
      
      store = rbind(store, data.frame(time = i, genes = list(c(genes)), accuracy.train = fit_rf_default$results$Accuracy))
      
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T,fontsize_number = 14,main = paste0("rf_",filename,i))
      rf_prob=predict(fit_rf_default, newdata = testing, type = "prob")
      # varImp = varImp(fit_rf_default); 
      # varImp = data.frame(varImp$importance)
      
      
      # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
      #   geom_point() +
      #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
      #   scale_color_discrete(name="Variable Group") +
      #   ylab("IncNodePurity") +
      #   xlab("Variable Name") +
      #   coord_flip()
      # plot(varImp, top = 20)
      ggsave(p, filename=paste0("./analysis_rhapsody_15rndgenes/MLfit_rf_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(fit_rf_default, paste0("./analysis_rhapsody_15rndgenes/MLfit_rf_",filename,"_",genesetname,"_",i, ".rds")) #on poseidon then moved back
    }
    
    if (0){ # SVM -------------------------------------------------------------------
      
      trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
      set.seed(1)
      
      svm_Linear <- train(label ~., data = training, 
                          method = "svmLinear",
                          trControl=trctrl,
                          # preProcess = c("center", "scale"),
                          tuneLength = 10)
      
      svm_Linear
      test_pred <- predict(svm_Linear, newdata = testing)
      test_pred
      confusion = confusionMatrix(test_pred, as.factor(testing$label) )
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T,fontsize_number = 14,main = paste0("svmL_",filename,i)) 
      # varImp = varImp(svm_Linear); #varImp = varImp$importance
      # plot(varImp, top = 20)
      ggsave(p,filename=paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
      
    }
    
    if (0){ # SVM radial -------------------------------------------------------------------
      
      trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
      set.seed(1)
      
      svm_Radial <- train(label ~., data = training,
                          method = "svmRadial",
                          trControl = trctrl,
                          # preProcess = c("center","scale"),
                          tuneLength = 10)
      
      # Print the best tuning parameter sigma and C that maximizes model accuracy
      svm_Radial$bestTune
      
      svm_Radial
      test_pred <- predict(svm_Radial, newdata = testing)
      test_pred
      confusion = confusionMatrix(test_pred, as.factor(testing$label) )
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T,fontsize_number = 14,main = paste0("svmR_",filename,i)) 
      # varImp = varImp(svm_Radial); #varImp = varImp$importance
      # plot(varImp, top = 20)
      # saveRDS(svm_Radial, paste0("./analysis_rhapsody/MLfit_svmRadial_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
      ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
      
    }
    
    if (0){ # kNN-------------------------------------------------------------------
      # my.dataframe = data.frame(cbind(label = colnames(data), macro@reductions$pca@cell.embeddings))
      # set.seed(1)
      # inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
      # training <- my.dataframe[ inTraining,]
      # testing  <- my.dataframe[-inTraining,]
      ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
      knnFit <- train(label ~ ., 
                      data = training, 
                      method = "knn", 
                      trControl = ctrl, 
                      # preProcess = c("center","scale"), 
                      tuneLength = 20)
      
      #Output of kNN fit
      knnFit
      #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
      plot(knnFit)
      knnPredict <- predict(knnFit,newdata = testing )
      #Get the confusion matrix to see accuracy value and other parameter values
      confusion = confusionMatrix(knnPredict, as.factor(testing$label) )
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T,fontsize_number = 14,main = paste0("knn_",filename,i)) 
      # varImp = varImp(knnFit); varImp = varImp$importance
      # saveRDS(knnFit, paste0("./analysis_rhapsody/MLfit_knn_M0all_gt80_500genes_",i, ".rds")) #run on poseidon for mememory problems
      ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(knnFit, paste0("./analysis_rhapsody_500genes/MLfit_knn_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
      
    }
    
    if (0){ #naive bayes
      library(e1071)
      ctrl <- trainControl(method="repeatedcv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
      nbFit <- train(label ~ ., 
                     data = training, 
                     method = "nb", 
                     trControl = ctrl, 
                     # preProcess = c("center","scale"),
                     tuneLength = 20)
      
      #Output of nb fit
      nbFit
      #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
      plot(nbFit)
      nbPredict <- predict(nbFit,newdata = testing )
      #Get the confusion matrix to see accuracy value and other parameter values
      confusion = confusionMatrix(nbPredict, as.factor(testing$label) )
      confusion.table = (confusion$table)
      confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
                 # colorRampPalette(c("lightblue", "white", "red"))(50),
                 # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
                 breaks =  seq(0, 1, by = .01),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
                 display_numbers = T, fontsize_number = 14,main = paste0("nb_",filename,i)) 
      # varImp = varImp(nbFit); 
      # varImp = varImp$importance
      ggsave(p, filename=paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".png"))
      saveRDS(nbFit, paste0("./analysis_rhapsody_500genes/MLfit_nb_",filename,"_",genesetname,"_",i, ".rds")) #run on poseidon for mememory problems
      
    }
  }
}
ggplot(store[grepl("", store$time),], aes(x=accuracy.train, color=time, fill=time)) +
  geom_density(alpha=0.2)+  theme_bw(base_size = 20)+xlim(0,1)

# Figure 6e----
# plot FDR and FPR----
#collect FPR and FDR for each stimulus from rnd15
table = data.frame()
macro.all = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
set.seed(1)
for (i in c("3hr") ){
  print(i)
  for (r in c(1:50)){
    print(r)
    filename = "M0all.Dec2020"
    genesetname = paste0("15genes.rnd",r)
    
    fit_rf_default = readRDS(paste0("./analysis_rhapsody_15rndgenes/MLfit_rf_",filename,"_",genesetname,"_",i, ".rds"))
    macro = subset(macro.all, subset= timept==i)
    data = macro[["ISnorm"]]@data
    data = data.frame(data)
    meta = macro@meta.data
    colnames(data) = paste0(meta$stimulus, "_",meta$timept)
    my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
    
    library(caret)
    set.seed(1)
    inTraining <- createDataPartition(my.dataframe$label, p = .7, list = FALSE)
    training <- my.dataframe[ inTraining,]
    testing  <- my.dataframe[-inTraining,]
    
    print(fit_rf_default)
    rfClasses <- predict(fit_rf_default, newdata = testing)
    confusion = confusionMatrix(data = rfClasses, as.factor(testing$label))
    confusion.table = (confusion$table)
    values = as.data.frame(confusion$byClass)
    values$rnd = r
    table = rbind(table, values)
    
  }
}
table$fpr = 1-table$Specificity
table$fdr = 1-table$Precision
table$stim = c(
  rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 50))
table[is.na(table)] <-0

table.m = melt(table[,c(13:15)])
ggplot(table.m[grepl("fdr", table.m$variable),], aes(x=value, color=stim, fill=stim)) +
  geom_density(alpha=0.2)+  theme_bw(base_size = 20)+xlim(0,1)+
  geom_vline(xintercept = 0.1949, color = "#F8766D", linetype="dashed", size=1)+
  geom_vline(xintercept = 0.0013, color = "#B79F00", linetype="dashed", size=1)+
  geom_vline(xintercept = 0.1499, color = "#00BA38", linetype="dashed", size=1)+
  geom_vline(xintercept = 0.2513, color = "#00BFC4", linetype="dashed", size=1)+
  geom_vline(xintercept = 0.1598, color = "#619CFF", linetype="dashed", size=1)+
  geom_vline(xintercept = 0.0972, color = "#F564E3", linetype="dashed", size=1)

ggplot(table.m[grepl("fdr", table.m$variable),], aes(x= stim, y=value, color=stim, fill=stim)) +
  geom_violin(alpha=0.2)+theme_bw(base_size = 20)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(aes(color = stim), position = position_jitter(w = 0.05, h = 0))+
  geom_hline(yintercept = 0.1949, color = "#F8766D", linetype="dashed", size=1)+
  geom_hline(yintercept = 0.0013, color = "#B79F00", linetype="dashed", size=1)+
  geom_hline(yintercept = 0.1499, color = "#00BA38", linetype="dashed", size=1)+
  geom_hline(yintercept = 0.2513, color = "#00BFC4", linetype="dashed", size=1)+
  geom_hline(yintercept = 0.1598, color = "#619CFF", linetype="dashed", size=1)+
  geom_hline(yintercept = 0.0972, color = "#F564E3", linetype="dashed", size=1)
# geom_hline(yintercept = 0.24, color = "#F8766D", linetype="dashed", size=1)+
# geom_hline(yintercept = 0.01, color = "#B79F00", linetype="dashed", size=1)+
# geom_hline(yintercept = 0.18, color = "#00BA38", linetype="dashed", size=1)+
# geom_hline(yintercept = 0.22, color = "#00BFC4", linetype="dashed", size=1)+
# geom_hline(yintercept = 0.19, color = "#619CFF", linetype="dashed", size=1)+
# geom_hline(yintercept = 0.07, color = "#F564E3", linetype="dashed", size=1)


###################################Figure 7


# Figure 6f----
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

mutual = readRDS("./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
genes = unlist(mutual[41,4])
macro = readRDS(paste0("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr.rds"))

collect = data.frame()
for (i in c("M0","M1_IFNg","M2_IL4")){
  print(i)
  for (s in c(
    "LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",
    "LPS|PIC|IFNb|TNF|P3CSK", "LPS|P3CSK|PIC","LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF", "LPS|P3CSK")){
    print(s)
    meta = na.omit(samptag.all)
    wanted = (meta$Cell_Index[meta$type==i])
    
    data = macro[["ISnorm"]]@data
    data = data.frame(data)
    
    my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
    my.dataframe = my.dataframe[rownames(my.dataframe)%in%wanted, colnames(my.dataframe) %in% c("label", genes)]
    
    
    my.dataframe$label = meta$stimulus[match(rownames(my.dataframe), rownames(meta))]
     str(my.dataframe)
    
    
    my.dataframe=my.dataframe[grepl(s, my.dataframe$label),]
    #--------------------------------mi using SLEMI -----------------------
    
    #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
    output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                            # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                            testing = T, boot_num = 50, boot_prob = .5 ,testing_cores = 4) #without0hr
    sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
    tmp = data.frame(type = i, stim_cells = s, cc = output_capacity$cc, sd = sd)
    collect = rbind(collect, tmp)
    
    # tempoutput_probs  <- prob_discr_pairwise(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
    #                                          output_path=paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_pairwiseProb_", i)) 
  }
}
write.table(collect, "./infotheo/collect_polarization.top15genes_3hr_bootstrap.txt", sep="\t",quote = F, row.names = F)

collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4" ))
collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",#"LPS|P3CSK",
                                                          "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
ggplot(na.omit(collect), aes(stim_cells, cc, fill = type))+geom_bar(stat="identity", position="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))

# Figure 6g----
#calculate MI on PC scores----
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

pc.scores = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt")
colnames(pc.scores)[1]="Sample"
pc.scores$stimulus = samptag.all$stimulus[match(pc.scores$Sample, samptag.all$Cell_Index)]
pc.scores$type = samptag.all$type[match(pc.scores$Sample, samptag.all$Cell_Index)]
projected.data = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_2M03hrs_prcomp_rotated.txt")
projected.data$stimulus = samptag.all$stimulus[match(projected.data$Sample, samptag.all$Cell_Index)]
projected.data$type = samptag.all$type[match(projected.data$Sample, samptag.all$Cell_Index)]

mi.frame = rbind(pc.scores, projected.data)
# mi.frame = projected.data
rownames(mi.frame)= mi.frame$Sample
mi.frame$Sample = mi.frame$stimulus
colnames(mi.frame)[1]="label"
# mi.frame$label = factor(mi.frame$label, levels=c("PIC","P3CSK","LPS","IFNb","TNF"))


collect = data.frame()
nPCs = 3
for (i in c("M0","M1_IFNg","M2_IL4")){
  print(i)
  for (s in c(
    "LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",
    "LPS|PIC|IFNb|TNF|P3CSK", "LPS|P3CSK|PIC","LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF", "LPS|P3CSK")){
    # for (s in c("PIC|TNF")){
    print(s)
    meta = na.omit(samptag.all)
    wanted = (meta$Cell_Index[meta$type==i])
    my.dataframe= mi.frame[rownames(mi.frame)%in%wanted, c(1:(nPCs+1))]
    
    my.dataframe=my.dataframe[grepl(s, my.dataframe$label),]
    #--------------------------------mi using SLEMI -----------------------
    
    #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
    output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                            # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                            testing = T, boot_num = 50, boot_prob = .5 ,testing_cores = 4) #without0hr
    sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
    tmp = data.frame(type = i, stim_cells = s, cc = output_capacity$cc, sd = sd)
    collect = rbind(collect, tmp)
    
    # tempoutput_probs  <- prob_discr_pairwise(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
    #                                          output_path=paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_pairwiseProb_", i)) 
  }
}

collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
colors_list = (c(M0="gray",M1_IFNg="gray",M2_IL4="gray",PM_B6.LFD="#F8766D",PM_B6.old="#00BA38",PM_B6.HFD="#619CFF"))
ggplot(collect, aes(type, cc, fill = type))+geom_bar(stat="identity", position="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  facet_grid(~stim_cells)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)
# write.table(collect, "./infotheo/collect_diseasemodel_bootstrap.txt", sep="\t",quote = F, row.names = F)
# write.table(collect, "./infotheo/collect_Polarization_BMDMdiseasemodel_bootstrap.txt", sep="\t",quote = F, row.names = F)

collect = read.delim("./infotheo/collect_diseasemodel_bootstrap.txt")
collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK", # "LPS|P3CSK","LPS|P3CSK|PIC",
                                                          "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
ggplot(na.omit(collect), aes(stim_cells, cc, fill = type))+geom_bar(stat="identity", position="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  scale_fill_manual(values=colors_list)


#plot just the M0/M1/M2 bars------------------
collect = read.delim("./infotheo/collect_Polarization_BMDMdiseasemodel_bootstrap.txt")
collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4"))
collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",#"LPS|P3CSK",
                                                          "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
colors_list = (c(M0="#F8766D",M1_IFNg="#00BA38",M2_IL4="#619CFF"))
ggplot(na.omit(collect[grepl("M0|M1|M2", collect$type),]), aes(stim_cells, cc, fill = type))+
  geom_bar(stat="identity", position ="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)+
  geom_hline(yintercept = 1, linetype="dotted")+
  # geom_hline(yintercept = log(3)/log(2), linetype="dotted")+
  geom_hline(yintercept = 2, linetype="dotted")


######################################## Figure 7----
# Figure 7d ----
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


# Figure 7ef----
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
      
      # dist = data.frame(condition = paste0(i, "vsNot",i), bd = bd, type = type) # use this for Figure 7e
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

# Figure 7g----
#project disease new Feb 2021--------------------
intersect_doPCA_from_file_and_project_second_dataset("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr.txt",
                                                     "./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC.txt", train_string = "2M0M1M23hrs",
                                                     center = T, scale = F)
# varimax_from_file("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt",
#                   "./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_loadings.txt",comp = 3, normalize = F)

# debug(plot_pca)
plot_pca("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt", samptag.all$Cell_Index, samptag.all$type, pt.size = 0.1,ellipse = F, labels = F)
plot_pca_projection("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt", 
                    "./output/macrophage_PMexpts_Feb2021_500genes_DBEC_2M0M1M23hrs_prcomp_rotated.txt",
                    samptag.all$Cell_Index, samptag.all$stimulus, samptag.all$Cell_Index, samptag.all$stimulus)
plot_pca_projection_all("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt", 
                        "./output/macrophage_PMexpts_Feb2021_500genes_DBEC_2M0M1M23hrs_prcomp_rotated.txt",
                        samptag.all$Cell_Index, as.factor(samptag.all$type), 
                        samptag.all$Cell_Index, as.factor(samptag.all$stimulus), ellipse =T)

# plot projected
projected.data = read.delim("./output/macrophage_PMexpts_Feb2021_500genes_DBEC_2M0M1M23hrs_prcomp_rotated.txt")
projected.data$stimulus = samptag.all$stimulus[match(projected.data$Sample, samptag.all$Cell_Index)]
projected.data$type = samptag.all$type[match(projected.data$Sample, samptag.all$Cell_Index)]
projected.data$type = factor(projected.data$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
p1=ggplot(projected.data, aes(PC1, PC2))+geom_point(aes(color = stimulus),size=0.01)+theme_bw(base_size = 16)
p2=ggplot(projected.data, aes(PC1, PC2))+geom_point(aes(color = type),size=0.01)+theme_bw(base_size = 16)
p1|p2

data = rbind(pc.scores, projected.data)
rownames(data)= data$Sample
set.seed(123)
umap.projected = umap(data[,c(2:21)], n_neighbors = 30L,
                      min_dist = 0.3, seed.use = 42)
umap.projected.layout = data.frame(umap.projected$layout)
umap.projected.layout$stimulus = samptag.all$stimulus[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
umap.projected.layout$type = samptag.all$type[match(rownames(umap.projected.layout), samptag.all$Cell_Index)]
# umap.projected.layout$type = factor(umap.projected.layout$type, levels= c("PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
p1=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = stimulus),size=0.01)+theme_bw(base_size = 16)
p2=ggplot(umap.projected.layout, aes(X1, X2))+geom_point(aes(color = type),size=0.01)+theme_bw(base_size = 16)
p1|p2

#caluclate MI on PC scores
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all$Cell_Index = paste0("X", samptag.all$Cell_Index)

pc.scores = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_prcomp_scores.txt")
colnames(pc.scores)[1]="Sample"
pc.scores$stimulus = samptag.all$stimulus[match(pc.scores$Sample, samptag.all$Cell_Index)]
pc.scores$type = samptag.all$type[match(pc.scores$Sample, samptag.all$Cell_Index)]
projected.data = read.delim("./output/macrophage_PMexpts_Feb2021_rmUnstim_500genes_DBEC_2M0M1M23hrs_prcomp_rotated.txt")
# projected.data = read.delim("./output/macrophage_M0M1M2_combined_500genes_DBEC_3hr_2M03hrs_prcomp_rotated.txt")
projected.data$stimulus = samptag.all$stimulus[match(projected.data$Sample, samptag.all$Cell_Index)]
projected.data$type = samptag.all$type[match(projected.data$Sample, samptag.all$Cell_Index)]

mi.frame = rbind(pc.scores, projected.data)
# mi.frame = projected.data
rownames(mi.frame)= mi.frame$Sample
mi.frame$Sample = mi.frame$stimulus
colnames(mi.frame)[1]="label"
# mi.frame$label = factor(mi.frame$label, levels=c("PIC","P3CSK","LPS","IFNb","TNF"))


collect = data.frame()
nPCs = 3
for (i in c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD")){
# for (i in c("M0","M1_IFNg","M2_IL4")){
  print(i)
  for (s in c(
    "LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",
    "LPS|PIC|IFNb|TNF|P3CSK", "LPS|P3CSK|PIC","LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF", "LPS|P3CSK")){
    # for (s in c("PIC|TNF")){
    print(s)
    meta = na.omit(samptag.all)
    wanted = (meta$Cell_Index[meta$type==i])
    my.dataframe= mi.frame[rownames(mi.frame)%in%wanted, c(1:(nPCs+1))]
    
    my.dataframe=my.dataframe[grepl(s, my.dataframe$label),]
    #--------------------------------mi using SLEMI -----------------------
    
    #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
    output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                            # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                            testing = T, boot_num = 50, boot_prob = .5 ,testing_cores = 4) #without0hr
    sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
    tmp = data.frame(type = i, stim_cells = s, cc = output_capacity$cc, sd = sd)
    collect = rbind(collect, tmp)
    
    # tempoutput_probs  <- prob_discr_pairwise(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
    #                                          output_path=paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_pairwiseProb_", i)) 
  }
}

collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
colors_list = (c(M0="gray",M1_IFNg="gray",M2_IL4="gray",PM_B6.LFD="#F8766D",PM_B6.old="#00BA38",PM_B6.HFD="#619CFF"))
ggplot(collect, aes(type, cc, fill = type))+geom_bar(stat="identity", position="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  facet_grid(~stim_cells)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)
# write.table(collect, "./infotheo/collect_diseasemodel_bootstrap.txt", sep="\t",quote = F, row.names = F)

collect = read.delim("./infotheo/collect_diseasemodel_bootstrap.txt")
collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK", 
                                                          "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
ggplot(na.omit(collect), aes(stim_cells, cc, fill = type))+geom_bar(stat="identity", position="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  scale_fill_manual(values=colors_list)


#plot just the M0/M1/M2 bars------------------
collect = read.delim("./infotheo/collect_Polarization_BMDMdiseasemodel_bootstrap.txt")
collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4"))
collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK|CpG", "LPS|P3CSK|CpG",
                                                          "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
colors_list = (c(M0="#F8766D",M1_IFNg="#00BA38",M2_IL4="#619CFF"))
ggplot(na.omit(collect[grepl("M0|M1|M2", collect$type),]), aes(stim_cells, cc, fill = type))+
  geom_bar(stat="identity", position ="dodge")+theme_bw(base_size = 16)+ylab("channel capacity")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)+
  geom_hline(yintercept = 1, linetype="dotted")+
  # geom_hline(yintercept = log(3)/log(2), linetype="dotted")+
  geom_hline(yintercept = 2, linetype="dotted")

#calc RSS score----
if(0){
  collect = read.delim("./infotheo/collect_diseasemodel_bootstrap.txt")
  collect$type = factor(collect$type, levels= c("M0","M1_IFNg","M2_IL4","PM_B6.LFD","PM_B6.old","PM_B6.HFD" ))
  collect$stim_cells = factor(collect$stim_cells, levels= c("LPS|PIC|IFNb|TNF|P3CSK", # "LPS|P3CSK","LPS|P3CSK|PIC",
                                                            "LPS|PIC|TNF", "LPS|PIC|IFNb", "IFNb|TNF", "PIC|TNF"))
  tmp = na.omit(collect)
}
tmp.cast=  dcast(tmp, stim_cells~type, value.var = 'cc')
tmp.cast=  cbind(tmp.cast, sd = (dcast(tmp, stim_cells~type, value.var = 'sd')[,-1])^2)# convert to variances

tmp.cast$M0.score = (tmp.cast$M0 - tmp.cast$M0)^2
tmp.cast$M1.score = (tmp.cast$M1_IFNg - tmp.cast$M0)^2
tmp.cast$M2.score = (tmp.cast$M2_IL4 - tmp.cast$M0)^2
tmp.cast$LFD.score = (tmp.cast$PM_B6.LFD - tmp.cast$M0)^2
tmp.cast$old.score = (tmp.cast$PM_B6.old - tmp.cast$M0)^2
tmp.cast$HFD.score = (tmp.cast$PM_B6.HFD - tmp.cast$M0)^2

tmp.cast$M0.sd = (tmp.cast$sd.M0 + tmp.cast$sd.M0)
tmp.cast$M1.sd = (tmp.cast$sd.M1_IFNg + tmp.cast$sd.M0)
tmp.cast$M2.sd = (tmp.cast$sd.M2_IL4 + tmp.cast$sd.M0)
tmp.cast$LFD.sd = (tmp.cast$sd.PM_B6.LFD + tmp.cast$sd.M0)
tmp.cast$old.sd = (tmp.cast$sd.PM_B6.old + tmp.cast$sd.M0)
tmp.cast$HFD.sd = (tmp.cast$sd.PM_B6.HFD + tmp.cast$sd.M0)

score = data.frame(score = sqrt(colSums(tmp.cast[,c(14:19)])),
                   sd = sqrt(colSums(tmp.cast[,c(20:25)])))
score$type = rownames(score)
score$type = factor(score$type, levels= c("M0.score","M1.score","M2.score","LFD.score","old.score","HFD.score" ))
colors_list = (c(M0.score="gray",M1.score="gray",M2.score="gray",LFD.score="#F8766D",old.score="#00BA38",HFD.score="#619CFF"))
ggplot(score, aes(type, score))+geom_bar(stat= "identity", aes(fill = type))+
  geom_errorbar(aes(ymin=score-sd, ymax=score+sd), width=.3,position=position_dodge(1))+
  scale_fill_manual(values=colors_list)+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

