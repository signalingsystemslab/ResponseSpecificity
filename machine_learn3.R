#Ml classification

library(caret);library(dplyr)
library(xgboost);library(randomForest);library(vcdExtra)
# library(ksheu.library1);
library(pheatmap);library(RColorBrewer);library(Seurat)

#from normalized and downsampled matrix----
# macro = readRDS("./output/macrophage_tutorial3_4sets_subsample_sctransform.rds")


#training on all timepts and stim together----
if(1){
  # macro.1 = readRDS("./output/macrophage_M2_IL4_500genes_DBEC.rds")
  # name = "M2"
  # macro = readRDS("./output/macrophage_M1_IFNg_500genes_DBEC.rds")
  # name = "M1"
  # macro = readRDS("./output/macrophage_M0all_500genes_DBEC.rds")
  # name = "M0"
  data = macro[["RNA"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  # rownames(data) = gsub("-ENSMUST..*","", macro[["SCT"]]@counts@Dimnames[[1]])
  # rownames(data) = gsub("..NM..*","", rownames(data))
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  #ML using CARET
  library(caret)
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
  training <- my.dataframe[ inTraining,]
  testing  <- my.dataframe[-inTraining,]
  testing = testing[!grepl("24hr|0.25hr|LPSlo_",testing$label), colnames(testing) %in% c("label",macro.1[["SCT"]]@counts@Dimnames[[1]]) ]
  
  if (0) {#plsr --------------------------------------------------------------------------------
    fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
    fit_pls <- train(label ~ .,  data = training,
                     method = "pls",
                     ## Center and scale the predictors for the training
                     ## set and all future samples.
                     preProc = c("center", "scale"),
                     tuneLength = 20,
                     trControl = fitControl,
                     metric = "ROC")
    fit_pls
    ggplot(fit_pls)
    plsClasses <- predict(fit_pls, newdata = testing)
    plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
    head(plsClasses)
    head(plsProbs)
    confusion = confusionMatrix(data = plsClasses, testing$label)
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
    # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
    
    varImp = varImp(fit_pls); varImp = varImp$importance
    
    # save the model to disk
    saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes/MLfit_pls_all",name,".rds"))
  }
  
  if (1){ #random forest------------------------------------------------------------------------
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
    confusion = confusionMatrix(data = rfClasses,
                                factor(testing$label, levels = levels(rfClasses)))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(fit_rf_default); 
    varImp = data.frame(varImp$importance)
    
    
    # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
    #   geom_point() +
    #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
    #   scale_color_discrete(name="Variable Group") +
    #   ylab("IncNodePurity") +
    #   xlab("Variable Name") +
    #   coord_flip()
    # plot(varImp, top = 20)
    saveRDS(fit_rf_default, paste0("./analysis_rhapsody_500genes/MLfit_rf_all_",name,".rds")) #on poseidon then moved back
  }
  
  if (0){ # SVM -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Linear <- train(label ~., data = training, 
                        method = "svmLinear",
                        trControl=trctrl,
                        preProcess = c("center", "scale"),
                        tuneLength = 10)
    
    svm_Linear
    test_pred <- predict(svm_Linear, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Linear); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_all_",name,".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ # SVM radial -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Radial <- train(label ~., data = training,
                        method = "svmRadial",
                        trControl = trctrl,
                        preProcess = c("center","scale"),
                        tuneLength = 10)
    
    # Print the best tuning parameter sigma and C that maximizes model accuracy
    svm_Radial$bestTune
    
    svm_Radial
    test_pred <- predict(svm_Radial, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Radial); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes/MLfit_svmRadial_all_",name,".rds")) #run on poseidon for mememory problems
    
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
                    preProcess = c("center","scale"), 
                    tuneLength = 20)
    
    #Output of kNN fit
    knnFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(knnFit)
    knnPredict <- predict(knnFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(knnPredict, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(knnFit); varImp = varImp$importance
    saveRDS(knnFit, paste0("./analysis_rhapsody_500genes/MLfit_knn_all_",name,".rds")) #run on poseidon for mememory problems
    
  }
} #ran on poseidon wdblack


#try on timepoints individually----
mutual = readRDS("./infotheo/collect_dimensionbest_ISnorm_Feb2021.rds")
genes1 = unlist(mutual[15,4])
genes3 = unlist(mutual[41,4])
genes8 = unlist(mutual[67,4])
genes = c(genes1, genes3, genes8)
genes = c(genes3, genes8)
genes = genes[!duplicated(genes)]

for (i in c("0.25hr","1hr","3hr","8hr")){
# for (i in c("3hr")){
  
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i, "_DBEC.rds"))
  macro = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
  # macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
  
  filename = "M0all.Dec2020"
  genesetname = "500genes" #"MItop15.38hrmerge" #"MItop15.1hr"
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

####################################################################################
# test random sets of 15 genes -----
store = data.frame()
set.seed(1)
for (i in c("0.25hr","1hr","3hr","8hr")){
  # for (i in c("3hr")){
  
  print(i)
  # macro = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i, "_DBEC.rds"))
  macro = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
  # macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
  
  for (r in c(1:50)){
    print(r)
    genes = sample(rownames(macro[["ISnorm"]]@data), 15)
  filename = "M0all.Dec2020"
  genesetname = paste0("15genes.rnd",r) #"500genes" #"MItop15.38hrmerge" #"MItop15.1hr"
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

#plot FDR and FPR----
setwd("F://scRNAseq_macro/scRNAseq_macro/")
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

#######################################################################################
  
#ensemble for M0 above knn, svm radial, rf,  nb-----
# rf_fit = readRDS("./analysis_rhapsody_500genes/MLfit_rf_M0all.Dec2020_500genes_1hr.rds")
# rf_fit = readRDS("./analysis_rhapsody_500genes/MLfit_rf_M0all_MItop15.1hr_1hr.rds")
# predProbs <- extractProb(list(rf_fit),
#                          testX = testing_1hr[,-1], testY = testing_1hr$label)
# predProbs = predProbs[order(predProbs$obs),]
# colors_list = list(obs = c(Unstim = "gray", CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3"),
#                    pred = c(Unstim = "gray", CpG_1hr="#F8766D", IFNb_1hr="#B79F00",LPS_1hr= "#00BA38",P3CSK_1hr= "#00BFC4",PIC_1hr= "#619CFF",TNF_1hr= "#F564E3"))
# pheatmap(predProbs[predProbs$dataType == "Test",c(1:6)], scale = "none",cluster_rows = F, cluster_cols = F,
#          annotation_row = predProbs[predProbs$dataType == "Test",c(7:8)],
#          show_rownames = F, main = "prediction probabilities",
#          annotation_colors = colors_list)
# plotClassProbs(predProbs)
# plotClassProbs(predProbs[predProbs$dataType == "Test",])


knn_fit = readRDS("./analysis_rhapsody_500genes/MLfit_knn_M0all_MItop15.1hr_1hr.rds")
rf_fit = readRDS("./analysis_rhapsody_500genes/MLfit_rf_M0all_MItop15.1hr_1hr.rds")
svmR_fit = readRDS("./analysis_rhapsody_500genes/MLfit_svmRadial_M0all_MItop15.1hr_1hr.rds")
nb_fit = readRDS("./analysis_rhapsody_500genes/MLfit_nb_M0all_MItop15.1hr_1hr.rds")
svmL_fit = readRDS("./analysis_rhapsody_500genes/MLfit_svmLinear_M0all_MItop15.1hr_1hr.rds")

knn_fit = readRDS("./analysis_rhapsody_500genes/MLfit_knn_M0all.noCpG_MItop15.38hrmerge_3hr.rds")
rf_fit = readRDS("./analysis_rhapsody_500genes/MLfit_rf_M0all.noCpG_MItop15.38hrmerge_3hr.rds")
svmR_fit = readRDS("./analysis_rhapsody_500genes/MLfit_svmRadial_M0all.noCpG_MItop15.38hrmerge_3hr.rds")
nb_fit = readRDS("./analysis_rhapsody_500genes/MLfit_nb_M0all.noCpG_MItop15.38hrmerge_3hr.rds")
svmL_fit = readRDS("./analysis_rhapsody_500genes/MLfit_svmLinear_M0all.noCpG_MItop15.38hrmerge_3hr.rds")

knn_dataframe = predict(knn_fit,newdata = testing )
rf_dataframe = predict(rf_fit,newdata = testing )
svmR_dataframe = predict(svmR_fit,newdata = testing )
nb_dataframe = predict(nb_fit,newdata = testing )
svmL_dataframe = predict(svmL_fit,newdata = testing )


ensemble_dataframe = cbind(knn_dataframe, rf_dataframe, svmR_dataframe, nb_dataframe, svmL_dataframe)
ensemble_dataframe = cbind(ensemble_dataframe, data.frame(consensus=apply(ensemble_dataframe,1,function(x) names(which.max(table(x))))))

confusionensemblehrs=confusionMatrix(reference = as.factor(as.numeric(as.factor(testing$label))), data = as.factor(ensemble_dataframe$consensus), mode='everything')
confusionensemblehrs
confusionensemblehrs.table = (confusionensemblehrs$table)
confusionensemblehrs.table = sweep(confusionensemblehrs.table, 2, colSums(confusionensemblehrs.table), FUN = '/')
pheatmap(confusionensemblehrs.table, cluster_rows = F, cluster_cols = F, 
         # colorRampPalette(c("lightblue", "white", "red"))(50),
         # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
         breaks =  seq(0, 1, by = .01),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
         display_numbers = T, fontsize_number = 14,main = paste0("ensemble5_",filename,i)) 


# load M0all.noCPG model and test disease PMs-----
rf_M0.nCpG.3hr = readRDS("./analysis_rhapsody_500genes/MLfit_rf_M0all.noCpG_MItop15.38hrmerge_3hr.rds")

macro = readRDS("./output/macrophage_PMexpts_Feb2021_500genes_DBEC.rds")
macro = subset(macro, subset= type=="PM_B6.HFD"&stimulus!="Unstim")
data = macro[["ISnorm"]]@data
data = data.frame(data)
meta = macro@meta.data
colnames(data) = paste0(meta$stimulus, "_",meta$timept)
my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
my.dataframe = my.dataframe[, colnames(my.dataframe) %in% c("label", genes)]
testing = my.dataframe

rf_dataframe = predict(rf_M0.nCpG.3hr,newdata = my.dataframe )
confusion = confusionMatrix(data = rf_dataframe, as.factor(my.dataframe$label))
confusion.table = (confusion$table)
confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
           # colorRampPalette(c("lightblue", "white", "red"))(50),
           # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
           breaks =  seq(0, 1, by = .01),
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
           display_numbers = T,fontsize_number = 14,main = paste0("rf_",filename,i)) 


#all timepoint individually for M1--------

for (i in c("0.5hr","1hr","3hr","5hr","8hr", "24hr")){
  print(i)
  macro = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_subsample_",i, "_DBEC.rds"))
  
  data = macro[["RNA"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  
  #--------------------------------------ML using CARET-----------------------------------------
  library(caret)
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
  training <- my.dataframe[ inTraining,]
  testing  <- my.dataframe[-inTraining,]
  
  if (0) {#plsr --------------------------------------------------------------------------------
    fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
    fit_pls <- train(label ~ .,  data = training,
                     method = "pls",
                     ## Center and scale the predictors for the training
                     ## set and all future samples.
                     preProc = c("center", "scale"),
                     tuneLength = 20,
                     trControl = fitControl,
                     metric = "ROC")
    fit_pls
    ggplot(fit_pls)
    plsClasses <- predict(fit_pls, newdata = testing)
    plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
    head(plsClasses)
    head(plsProbs)
    confusion = confusionMatrix(data = plsClasses, testing$label)
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
    # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
    
    varImp = varImp(fit_pls); varImp = varImp$importance
    
    # save the model to disk
    saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes//MLfit_pls_M1_IFNg_500genes_",i, ".rds"))
  }
  
  if (1){ #random forest------------------------------------------------------------------------
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
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T, fontsize_number = 16) 
    varImp = varImp(fit_rf_default); 
    varImp = data.frame(varImp$importance)
    
    
    # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
    #   geom_point() +
    #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
    #   scale_color_discrete(name="Variable Group") +
    #   ylab("IncNodePurity") +
    #   xlab("Variable Name") +
    #   coord_flip()
    # plot(varImp, top = 20)
    saveRDS(fit_rf_default, paste0("./analysis_rhapsody_500genes//MLfit_rf_M1_IFNg_500genes_",i, ".rds")) #on poseidon then moved back
  }
  
  if (0){ # SVM -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Linear <- train(label ~., data = training, 
                        method = "svmLinear",
                        trControl=trctrl,
                        preProcess = c("center", "scale"),
                        tuneLength = 10)
    
    svm_Linear
    test_pred <- predict(svm_Linear, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Linear); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_M1_IFNg_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ # SVM radial -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Radial <- train(label ~., data = training,
                        method = "svmRadial",
                        trControl = trctrl,
                        preProcess = c("center","scale"),
                        tuneLength = 10)
    
    # Print the best tuning parameter sigma and C that maximizes model accuracy
    svm_Radial$bestTune
    
    svm_Radial
    test_pred <- predict(svm_Radial, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Radial); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes//MLfit_svmRadial_M1_IFNg_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
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
                    preProcess = c("center","scale"), 
                    tuneLength = 20)
    
    #Output of kNN fit
    knnFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(knnFit)
    knnPredict <- predict(knnFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(knnPredict, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(knnFit); varImp = varImp$importance
    saveRDS(knnFit, paste0("./analysis_rhapsody_500genes//MLfit_knn_M1_IFNg_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
  }
}

#all timept individually for M2----------------------------------
for (i in c("0.5hr","1hr","3hr","5hr","8hr")){
  print(i)
  macro = readRDS(paste0("./output/macrophage_M2_IL4_500genes_subsample_",i, "_DBEC.rds"))
  
  data = macro[["RNA"]]@data
  data = data.frame(data)
  meta = macro@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$timept)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  
  #--------------------------------------ML using CARET-----------------------------------------
  library(caret)
  set.seed(1)
  inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
  training <- my.dataframe[ inTraining,]
  testing  <- my.dataframe[-inTraining,]
  
  if (0) {#plsr --------------------------------------------------------------------------------
    fitControl <- trainControl(method = "repeatedcv", repeats = 10, classProbs = T)
    fit_pls <- train(label ~ .,  data = training,
                     method = "pls",
                     ## Center and scale the predictors for the training
                     ## set and all future samples.
                     preProc = c("center", "scale"),
                     tuneLength = 20,
                     trControl = fitControl,
                     metric = "ROC")
    fit_pls
    ggplot(fit_pls)
    plsClasses <- predict(fit_pls, newdata = testing)
    plsProbs <- predict(fit_pls, newdata = testing, type = "prob")
    head(plsClasses)
    head(plsProbs)
    confusion = confusionMatrix(data = plsClasses, testing$label)
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    # grid.text("Actual", y=-0.01, gp=gpar(fontsize=16))
    # grid.text("Predicted", x=-0.02, rot=90, gp=gpar(fontsize=16))
    
    varImp = varImp(fit_pls); varImp = varImp$importance
    
    # save the model to disk
    saveRDS(fit_pls, paste0("./analysis_rhapsody_500genes//MLfit_pls_M2_IL4_500genes_",i, ".rds"))
  }
  
  if (1){ #random forest------------------------------------------------------------------------
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
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T, fontsize_number = 16) 
    varImp = varImp(fit_rf_default); 
    varImp = data.frame(varImp$importance)
    
    
    # ggplot(varImp, aes(x=reorder(varnames, IncNodePurity), y=IncNodePurity, color=as.factor(var_categ))) + 
    #   geom_point() +
    #   geom_segment(aes(x=varnames,xend=varnames,y=0,yend=IncNodePurity)) +
    #   scale_color_discrete(name="Variable Group") +
    #   ylab("IncNodePurity") +
    #   xlab("Variable Name") +
    #   coord_flip()
    # plot(varImp, top = 20)
    saveRDS(fit_rf_default, paste0("./analysis_rhapsody_500genes//MLfit_rf_M2_IL4_500genes_",i, ".rds")) #on poseidon then moved back
  }
  
  if (0){ # SVM -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Linear <- train(label ~., data = training, 
                        method = "svmLinear",
                        trControl=trctrl,
                        preProcess = c("center", "scale"),
                        tuneLength = 10)
    
    svm_Linear
    test_pred <- predict(svm_Linear, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Linear); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Linear, paste0("./analysis_rhapsody_500genes/MLfit_svmLinear_M2_IL4_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
  }
  
  if (0){ # SVM radial -------------------------------------------------------------------
    
    trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
    set.seed(1)
    
    svm_Radial <- train(label ~., data = training,
                        method = "svmRadial",
                        trControl = trctrl,
                        preProcess = c("center","scale"),
                        tuneLength = 10)
    
    # Print the best tuning parameter sigma and C that maximizes model accuracy
    svm_Radial$bestTune
    
    svm_Radial
    test_pred <- predict(svm_Radial, newdata = testing)
    test_pred
    confusion = confusionMatrix(test_pred, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(svm_Radial); #varImp = varImp$importance
    plot(varImp, top = 20)
    saveRDS(svm_Radial, paste0("./analysis_rhapsody_500genes//MLfit_svmRadial_M2_IL4_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
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
                    preProcess = c("center","scale"), 
                    tuneLength = 20)
    
    #Output of kNN fit
    knnFit
    #Plotting yields Number of Neighbours Vs accuracy (based on repeated cross validation)
    plot(knnFit)
    knnPredict <- predict(knnFit,newdata = testing )
    #Get the confusion matrix to see accuracy value and other parameter values
    confusion = confusionMatrix(knnPredict, testing$label )
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
             # colorRampPalette(c("lightblue", "white", "red"))(50),
             # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
             breaks =  seq(0, 1, by = .01),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
             display_numbers = T) 
    varImp = varImp(knnFit); varImp = varImp$importance
    saveRDS(knnFit, paste0("./analysis_rhapsody_500genes//MLfit_knn_M2_IL4_500genes_",i, ".rds")) #run on poseidon for mememory problems
    
  }
}


# load the model for training on all, need to load the all data and get test data too, above.
# did not run rf on all due to run time problems
fit_pls <- readRDS("./analysis_rhapsody/MLfit_pls_all.rds")
print(pls_model)
confusion.table = (confusion$table)
confusion.table = collapse.table(confusion.table, Reference=c(rep("CpG",5), rep("LPS",5), rep("LPSlo",5), rep("P3CSK",5), rep("PIC",5), rep("TNF",5), "Unstim"),
                                 Prediction = c(rep("CpG",5), rep("LPS",5), rep("LPSlo",5), rep("P3CSK",5), rep("PIC",5), rep("TNF",5), "Unstim" ) )
confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
         # colorRampPalette(c("lightblue", "white", "red"))(50),
         # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
         breaks =  seq(0, 1, by = .01),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
         display_numbers = T) 


#load all the models and the testing data----
# setwd("F://scRNAseq_macro/scRNAseq_macro/analysis_rhapsody_500genes/")

# Load machine learning models, and #read in by time point data, add labels to data
list = list.files(pattern= "MLfit_rf_M0all")
# for (i in c("knn", "pls", "rf", "svmLinear", "svmRadial") ){
for (i in c("rf") ){
  # for (j in c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") ){
  for (j in c("0.25hr","1hr","3hr","8hr") ){
    # tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M0all_MItop15.3hr_", j, ".rds"))
    tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M0all.Dec2020_500genes_", j, ".rds"))
    # tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M0all_gt80_500genes_", j, ".rds"))
    assign(paste0(i, "_", j), tmp)
  }
}

#load testing data for M1 M2 and test on M0 models----------
setwd("F://scRNAseq_macro/scRNAseq_macro/")
# for (i in c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") ){
# for (i in c("0.5hr","1hr","3hr","5hr","8hr") ){
# for (i in c("0.25hr","1hr","3hr","8hr") ){
# for (i in c("1hr","3hr","8hr") ){
for (i in c("3hr") ){
  
  # i = "0.5hr"
  # tmp = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_DBEC.rds"))
  # tmp = readRDS(paste0("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds"))
  tmp = readRDS(paste0("./output/macrophage_M0all_500genes_Dec2020.rds"))
  
  # tmp = readRDS(paste0("./output/macrophage_M0all_gt80_rmLPSlo_500genes_DBEC.rds"))
  # tmp = readRDS(paste0("./output/macrophage_M0_rep2only_500genes_DBEC.rds"))
  
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
  
  if(0){
    #0.5hr predicted on M0 0.25hr
    fit_rf_default = rf_0.25hr
    testing$label = gsub("0.5hr","0.25hr", testing$label)
  }
  
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


########################################################################
#load M1 models------
list = list.files(pattern= "MLfit_rf_M1_IFNg")
# for (i in c("knn", "pls", "rf", "svmLinear", "svmRadial") ){
for (i in c("rf") ){
  for (j in c("0.5hr","1hr","3hr","5hr","8hr","24hr") ){
    tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M1_IFNg_500genes_", j, ".rds"))
    assign(paste0(i, "_", j), tmp)
  }
}
#load testing data for M0 M2, and test on M1 models----
setwd("F://scRNAseq_macro/scRNAseq_macro/")
for (i in c("0.5hr","1hr","3hr","5hr","8hr") ){
  i = "8hr"
  # tmp = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_subsample_",i,"_DBEC.rds"))
  # tmp = readRDS(paste0("./output/macrophage_M2_IL4_500genes_subsample_",i,"_DBEC.rds"))
  tmp = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i,"_DBEC.rds"))
  tmp.meta = tmp@meta.data
  tmp2 = tmp[["RNA"]]@data
  tmp2 = data.frame(tmp2)
  colnames(tmp2) = paste0(tmp.meta$stimulus, "_",tmp.meta$timept)
  # rownames(tmp2) = gsub("-ENSMUST..*","", tmp[["SCT"]]@counts@Dimnames[[1]])
  # rownames(tmp2) = gsub("..NM..*","", rownames(tmp2))
  
  my.dataframe = cbind(label = colnames(tmp2), data.frame(t(tmp2)))
  # set.seed(1)
  # inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
  # training <- my.dataframe[ inTraining,]
  # assign(paste0("testing", "_",i), my.dataframe[-inTraining,])
  testing = my.dataframe
  
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
  
}


#load M2 models ---------------------------------------------
for (i in c("rf") ){
  for (j in c("0.5hr","1hr","3hr","5hr","8hr") ){
    tmp = readRDS(paste0("./analysis_rhapsody_500genes/MLfit_", i, "_M2_IL4_500genes_", j, ".rds"))
    assign(paste0(i, "_", j), tmp)
  }
}
#test other data on M2 models -------------------------------
for (i in c("0.5hr","1hr","3hr","5hr","8hr") ){
  # i = "1hr"
  # tmp = readRDS(paste0("./output/macrophage_M2_IL4_500genes_subsample_",i,"_DBEC.rds"))
  # tmp = readRDS(paste0("./output/macrophage_M1_IFNg_500genes_subsample_",i,"_DBEC.rds"))
  tmp = readRDS(paste0("./output/macrophage_M0all_gt80_500genes_subsample_",i,"_DBEC.rds"))
  tmp.meta = tmp@meta.data
  tmp2 = tmp[["RNA"]]@data
  tmp2 = data.frame(tmp2)
  colnames(tmp2) = paste0(tmp.meta$stimulus, "_",tmp.meta$timept)
  # rownames(tmp2) = gsub("-ENSMUST..*","", tmp[["SCT"]]@counts@Dimnames[[1]])
  # rownames(tmp2) = gsub("..NM..*","", rownames(tmp2))
  
  my.dataframe = cbind(label = colnames(tmp2), data.frame(t(tmp2)))
  # set.seed(1)
  # inTraining <- createDataPartition(my.dataframe$label, p = .3, list = FALSE)
  # training <- my.dataframe[ inTraining,]
  # assign(paste0("testing", "_",i), my.dataframe[-inTraining,])
  testing = my.dataframe
  
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
  
}




########################################################################
#plot total accuracy over time for rf----
# sensitivity = data.frame(0)
# specificity = data.frame(0)
# pospredval = data.frame(0)
table = data.frame()
# for (i in c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") ){
# for (i in c("0.5hr","1hr","3hr","5hr","8hr") ){
# for (i in c("0.25hr","1hr","3hr","8hr") ){
# for (i in c("1hr","3hr","8hr") ){
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


#for testing M0----
table$time = c(
  # rep(0.25, 6), 
               # rep(0.5, 5), 
               # rep(1, 6), 
               rep(3, 6) 
               # rep(5, 5), 
               # rep(8, 6)
               )
table$stim = c(
  # rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 1),
               # rep(c("CpG", "LPS",  "P3C", "PIC", "TNF"), 1),
               # rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 2),
               # rep(c("CpG", "LPS",  "P3C", "PIC", "TNF"), 1),
               rep(c("CpG","IFNb", "LPS",  "P3C", "PIC", "TNF"), 1))
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

ggplot(table, aes(as.factor(time), fpr, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+ylim(c(0,.5))

ggplot(table, aes(time, F1, group = stim))+ geom_point(aes(color = stim), size = 3)+
  geom_line(aes(time, F1, group = stim, color = stim), size = 1)+ylim(0,1)+xlim(0,8.2)+ theme_bw(base_size = 20)
ggplot(table, aes(time, Specificity, group = stim))+ geom_point(aes(color = stim), size = 3)+
  geom_line(aes(time, Specificity, group = stim, color = stim), size = 1)+ylim(0,1)+xlim(0,8.2)+ theme_bw(base_size = 20)

table.m = melt(table[grepl("3",table$time),c(13:15)])
ggplot(table.m, aes(variable, value, fill = stim))+ geom_bar(position="dodge", stat="identity")+
  theme_bw(base_size = 20)+ylim(c(0,1))

#for testing M1/M2-----
table$time = c(
  # rep(0.5, 6), 
  # rep(1, 6),
  rep(3, 6))
               # rep(5, 6), 
               # rep(8, 6))
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

ggplot(table[!grepl("Class: IFNb_5hr", rownames(table)),], aes(time, F1, group = stim))+ geom_point(aes(color = stim), size = 2)+
  geom_line(aes(time, F1, group = stim, color = stim), size = 1)+ylim(0,1)+xlim(0,8.2)+ theme_bw(base_size = 20)
ggplot(table[!grepl("Class: IFNb_5hr", rownames(table)),], aes(time, `Balanced Accuracy`, group = stim))+ geom_point(aes(color = stim), size = 2)+
  geom_line(aes(time, `Balanced Accuracy`, group = stim, color = stim), size = 1)+ylim(0,1)+xlim(0,8.2)+ theme_bw(base_size = 20)


# cv accuarcy
# acc0.5 = rf_0.5hr$results$Accuracy
# acc1 = rf_1hr$results$Accuracy
# acc3 = rf_3hr$results$Accuracy
# acc5 = rf_5hr$results$Accuracy
# acc8 = rf_8hr$results$Accuracy




#plot changing variable importances for rf----

# for (i in c("0.25hr","0.5hr","1hr","3hr","5hr","8hr") ){
for (i in c("3hr") ){
  print(i)
  fit_rf_default = get(paste0("rf_",i) )
  # testing = get( paste0("testing_",i) )
  
  print(fit_rf_default)
  # rfClasses <- predict(fit_rf_default, newdata = testing)
  # confusion = confusionMatrix(data = rfClasses, testing$label)
  # confusion.table = (confusion$table)
  # confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
  # pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
  #          # colorRampPalette(c("lightblue", "white", "red"))(50),
  #          # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
  #          breaks =  seq(0, 1, by = .01),
  #          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
  #          display_numbers = T) 
  varImp = varImp(fit_rf_default); 
  varImp = data.frame(varImp$importance)
  assign(paste0("varImp_",i), varImp)
}
# write.table(varImp_3hr, "./../SuppTables/TableS2_ML_varImp3hrs.txt", row.names = F, quote = F, sep = "\t")

intersect.genes = intersect_all(rownames(varImp_0.25hr),
                                rownames(varImp_0.5hr), rownames(varImp_1hr),rownames(varImp_3hr),
                                rownames(varImp_5hr),rownames(varImp_8hr))
collect = cbind(varImp_0.25hr[rownames(varImp_0.25hr) %in% intersect.genes, ], 
                varImp_0.5hr[rownames(varImp_0.5hr) %in% intersect.genes, ], 
                varImp_1hr[rownames(varImp_1hr) %in% intersect.genes, ],
                varImp_3hr[rownames(varImp_3hr) %in% intersect.genes, ],
                varImp_5hr[rownames(varImp_5hr) %in% intersect.genes, ],
                varImp_8hr[rownames(varImp_8hr) %in% intersect.genes, ])
rownames(collect) = intersect.genes; colnames(collect) = c("0.25hr","0.5hr","1hr","3hr","5hr","8hr")
pheatmap(collect[rowSums(collect)>100, ], cluster_rows = T, cluster_cols = F, scale = 'none', 
         show_rownames = T, show_colnames = T, clustering_method = "ward.D2") #for M0
pheatmap(collect[rowSums(collect)>80, ], cluster_rows = T, cluster_cols = F, scale = 'none', 
         show_rownames = T, show_colnames = T, clustering_method = "ward.D2") #for M1/M2


mat.melt = melt(collect)
ggplot(mat.melt[grepl("Il1b|Irf7|Adgre1|Ccl5", mat.melt$Var1),], 
       # aes(x=rep(c(0.5, 1, 3, 5, 8), nrow(collect[grepl("^Tnf$|^Cxcl10$|Il1b|Irf7", rownames(collect)),])), y=value,
       aes(x=Var2, y=value, group = Var1))+ 
  xlab("Time (hrs)")+ ylab("Importance")+
  geom_point()+  geom_line(aes(color = Var1), size=1)+theme_bw(base_size = 18)
ggplot(mat.melt[grepl("Cxcl10|Ccl3|Tnfrsf1b|Ddx58", mat.melt$Var1),], 
       aes(x=Var2, y=value, group = Var1))+ 
  xlab("Time (hrs)")+ ylab("Importance")+
  geom_point()+  geom_line(aes(color = Var1), size=1)+theme_bw(base_size = 18)
ggplot(mat.melt[grepl("^Tnf$|^Cxcl2$", mat.melt$Var1),], 
       aes(x=Var2, y=value, group = Var1))+ 
  xlab("Time (hrs)")+ ylab("Importance")+ylim(c(0,100))+
  geom_point()+  geom_line(aes(color = Var1), size=1)+theme_bw(base_size = 18)



# EXTRA/NOT YET WORKING STUFF----
# stochastic gradient boosting------------------------------------------------------------------
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 3) ## repeated ten times
fit_gbm <- train(label ~ ., data = training, 
                 method = "gbm", 
                 trControl = fitControl,
                 verbose = T)
fit_gbm



# neural net----------------------------------------------------------------------------------
library(neuralnet)
library(boot)
library(plyr)

# train on dim reduced data----
data.reduc = data.frame(macro@reductions$pca@cell.embeddings)
my.dataframe = cbind(label = rownames(data.reduc), data.frame( (data.reduc) ))

## Scale data for neural network
# max = apply(my.dataframe[,-1] , 2 , max)
# min = apply(my.dataframe[,-1], 2 , min)
# scaled = cbind(label = my.dataframe$label, as.data.frame(scale(my.dataframe[,-1]) )) #z-score

# train and test sets
set.seed(1)
train = sample(seq_len ( nrow ( my.dataframe ) ), 0.03*nrow(my.dataframe) )
datatrain = my.dataframe[train, ]
datatest = my.dataframe[-train, ]

# fit neural network
nn=neuralnet(label~ . , data=datatrain, hidden=2, act.fct = "logistic",
             linear.output = FALSE)
plot(nn)

## Prediction using neural network
predict_testNN = neuralnet::compute(nn, datatest)
results <- data.frame(actual = datatest$label, prediction = nn.results$net.result)

#confusion matrix
roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
attach(roundedresultsdf)
table(actual,prediction)

# xgboost-----------------------------------------------------------------------------------------
set.seed(1337) 
inTrain <- createDataPartition(y = my.dataframe$label, p = 0.85, list = FALSE)
x_train = xgb.DMatrix(as.matrix(my.dataframe[inTrain, ] %>% select(-label)))
y_train = my.dataframe[inTrain, ]$label
x_test = xgb.DMatrix(as.matrix(my.dataframe[-inTrain, ] %>% select(-label)))
y_test = my.dataframe[-inTrain, ]$label


param_search <- function(xtrain, ytrain, xtest, ytest){ 
  # Cross validation init
  xgb_trcontrol = trainControl(method = "cv", number = 5, allowParallel = TRUE, 
                               verboseIter = T, returnData = FALSE)
  # Param grid
  xgbGrid <- expand.grid(nrounds = 50, #nrounds = c(10,20,30,40),
                         max_depth = 30, #max_depth = c(3, 5, 10, 15, 20, 30),
                         colsample_bytree = 0.6, #colsample_bytree = seq(0.5, 0.9, length.out = 5),
                         eta = 0.005, #eta = c(0.001, 0.0015, 0.005, 0.1),
                         gamma=0, min_child_weight = 1, subsample = 1)
  # xgbGrid <- expand.grid(nrounds = c(10,20,30,40),
  #                        max_depth = c(3, 5, 10, 15, 20, 30),
  #                        colsample_bytree = seq(0.5, 0.9, length.out = 5),
  #                        eta = c(0.001, 0.0015, 0.005, 0.1),
  #                        gamma=0, min_child_weight = 1, subsample = 1)
  # Model and parameter search
  xgb_model = train(x_train, y_train, trControl = xgb_trcontrol,
                    tuneGrid = xgbGrid, method = "xgbTree",
                    verbose=2,
                    objective="multi:softprob",
                    eval_metric="mlogloss")
  #num_class=3)
  # Evaluate new model
  xgb.pred = predict(xgb_model, x_test, reshape=T)
  xgb.pred = as.data.frame(xgb.pred, col.names=c("pred"))
  result = sum(xgb.pred$xgb.pred==y_test) / nrow(xgb.pred)
  print(paste("Final Accuracy =",sprintf("%1.2f%%", 100*result)))
  
  return(xgb_model)
}

result.table = (table(xgb.pred$xgb.pred, y_test))
result.table = sweep(result.table, 2, colSums(result.table), FUN = '/')
pheatmap(result.table, cluster_rows = F, cluster_cols = F, 
         gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30))
pheatmap(result.table, cluster_rows = F, cluster_cols = F, 
         gaps_col = c(2,4,6,8,10,12), gaps_row = c(2,4,6,8,10,12))

# pheatmap(result.table, cluster_rows = F, cluster_cols = F)


best.model <- xgboost(
  data = as.matrix(my.dataframe[inTrain, ] %>% select(-IMPORTANCE)),
  label = as.matrix(as.numeric(my.dataframe[inTrain,]$IMPORTANCE)-1),
  nrounds = xgb_model$bestTune$nrounds,
  max_depth = xgb_model$bestTune$max_depth,
  eta = xgb_model$bestTune$eta,
  gamma = xgb_model$bestTune$gamma,
  colsample_bytree = xgb_model$bestTune$colsample_bytree,
  min_child_weight = xgb_model$bestTune$min_child_weight,
  subsample = xgb_model$bestTune$subsample,
  objective = "multi:softprob", num_class=3)

xgb_feature_imp <- xgb.importance(
  colnames(donnees[inTrain, ] %>% select(-label)), 
  model = best.model
)


gg <- xgb.ggplot.importance(xgb_feature_imp, 40); gg




# LASSO to pick inportant features---------------------------------------------------------------------------

# Loaging the library
library(glmnet)
lambda_seq = as.numeric(10^seq(2, -2, by = -.1))
# data.modelmatrix = model.matrix(label~., my.dataframe)

# Splitting the data into test and train
set.seed(1337)
inTrain <- createDataPartition(y = my.dataframe$label, p = 0.35, list = FALSE)
x_train = model.matrix(label~., my.dataframe[inTrain, ])[,-1]
y_train = my.dataframe[inTrain, ]$label
x_test = model.matrix(label~., my.dataframe[-inTrain, ])[,1]
y_test = my.dataframe[-inTrain, ]$label


# cv_output <- cv.glmnet(x=(x_train), y=(y_train), alpha = 1, nfolds = 5, lambda = lambda_seq, family = "multinomial")

# identifying best lamda
best_lam <- cv_output$lambda.min
best_lam = 0.01

# Rebuilding the model with best lamda value identified
lasso_best <- glmnet(x_train, (y_train), alpha = 1, lambda = best_lam, family = "multinomial")
pred <- predict(lasso_best, s = best_lam, newx = x_test)

coef(lasso_best)


final <- data.frame(cbind(y_test, pred))
# Checking the first six obs
head(final)


actual <- final$y_test
preds <- final$X1
rss <- sum((preds - actual) ^ 2)
tss <- sum((actual - mean(actual)) ^ 2)
rsq <- 1 - rss/tss
rsq

