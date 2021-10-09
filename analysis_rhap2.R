# analysis of rhapsody scRNAseq data

# install.packages("brms");
# library(brms)
library(stringi); library(reshape2);library(tidyr);library(dplyr);library(Matrix)
library(ksheu.library1);library(M3C) #library(tsne); 
setwd("F:/scRNAseq_macro/scRNAseq_macro/")
# data_dir = "./data/rhapsody/"
# out_dir = "./analysis_rhapsody/"

data_dir = "./data/rhapsody_500genes/"
out_dir = "./analysis_rhapsody_500genes/"


# load all the data and metadata ------------------------------------------

# new set OCt 2020 - 500 genes --------------------------------------------
mrna1 = read.delim("data/rhapsody_500genes/Combined_SampleTag1-Rev-Primer1_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#1 expt
mrna2 = read.delim("data/rhapsody_500genes/Combined_SampleTag2-Rev-Primer2_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#2
mrna3 = read.delim("data/rhapsody_500genes/Combined_SampleTag3-Rev-Primer3_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#3
mrna4 = read.delim("data/rhapsody_500genes/Combined_SampleTag4-Rev-Primer4_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#4
mrna5 = read.delim("data/rhapsody_500genes/Combined_SampleTag5-Rev-Primer1_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#5
mrna6 = read.delim("data/rhapsody_500genes/Combined_SampleTag6-Rev-Primer4_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#6 expt
mrna7 = read.delim("data/rhapsody_500genes/Combined_SampleTag7-Rev-Primer1_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#7
mrna8 = read.delim("data/rhapsody_500genes/Combined_SampleTag8-Rev-Primer2_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#8
mrna9 = read.delim("data/rhapsody_500genes/Combined_SampleTag9-Rev-Primer3_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#9
mrna10 = read.delim("data/rhapsody_500genes/Combined_SampleTag10-Rev-Primer4_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#10
mrna14 = read.delim("data/rhapsody_500genes/Combined_scRhapsodyLibraryPooled_PM_20201113_mRNA3-Rev_Primer3_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate#11PMs Dec2020

# Feb 2021 - PMs----------------------------------------------------
mrna15 = read.delim("data/rhapsody_500genes/Combined_mRNAPMs1-Rev-Primer1_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #PMs 3hrsfatVSnormal Feb2021
mrna16 = read.delim("data/rhapsody_500genes/Combined_mRNAPMs2-Rev-Primer2_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #PMs 3hrs old Feb2021

#original
mrna2019.1 = read.delim("data/rhapsody_500genes/Combined_SampleTag2-2019sample-Rev-Primer2_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate2019, 8,0.5
mrna2019.2 = read.delim("data/rhapsody_500genes/Combined_SampleTag3-2019sample-Rev-Primer3_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #Plate2019, 5,3
mrna3hrXtra = read.delim("data/rhapsody_500genes/_2_Combined_SampleTag4-Rev-Primer4_DBEC_MolsPerCell.csv", comment.char = "#", sep = ",") #PlateJan2020, 3hrXtras


samptag1 = read.delim("data/rhapsody_500genes/SampleTag1-Rev-Primer1_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag2 = read.delim("data/rhapsody_500genes/SampleTag2-Rev-Primer2_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag3 = read.delim("data/rhapsody_500genes/SampleTag3-Rev-Primer3_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag4 = read.delim("data/rhapsody_500genes/SampleTag4-Rev-Primer4_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag5 = read.delim("data/rhapsody_500genes/SampleTag5-Rev-Primer1_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag6 = read.delim("data/rhapsody_500genes/SampleTag6-Rev-Primer4_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag7 = read.delim("data/rhapsody_500genes/SampleTag7-Rev-Primer1_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag8 = read.delim("data/rhapsody_500genes/SampleTag8-Rev-Primer2_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag9 = read.delim("data/rhapsody_500genes/SampleTag9-Rev-Primer3_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag10 = read.delim("data/rhapsody_500genes/SampleTag10-Rev-Primer4_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag14 = read.delim("data/rhapsody_500genes/scRhapsodyLibraryPooled_PM_20201113_mRNA3-Rev_Primer3_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")

samptag15 = read.delim("data/rhapsody_500genes/mRNAPMs1-Rev-Primer1_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag16 = read.delim("data/rhapsody_500genes/mRNAPMs2-Rev-Primer2_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")


samptag2019.1 = read.delim("data/rhapsody_500genes/SampleTag2-2019sample-Rev-Primer2_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag2019.2 = read.delim("data/rhapsody_500genes/SampleTag3-2019sample-Rev-Primer3_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")
samptag3hrXtra = read.delim("data/rhapsody_500genes/_2_SampleTag4-Rev-Primer4_Sample_Tag_Calls.csv", comment.char = "#", sep = ",")


# concatenate--------------------------------------------------------------

mrna1$Cell_Index = paste0("1_",mrna1$Cell_Index)
mrna2$Cell_Index = paste0("2_",mrna2$Cell_Index)
mrna3$Cell_Index = paste0("3_",mrna3$Cell_Index)
mrna4$Cell_Index = paste0("4_",mrna4$Cell_Index)
mrna5$Cell_Index = paste0("5_",mrna5$Cell_Index)
mrna6$Cell_Index = paste0("6_",mrna6$Cell_Index)
mrna7$Cell_Index = paste0("7_",mrna7$Cell_Index)
mrna8$Cell_Index = paste0("8_",mrna8$Cell_Index)
mrna9$Cell_Index = paste0("9_",mrna9$Cell_Index)
mrna10$Cell_Index = paste0("10_",mrna10$Cell_Index)
mrna14$Cell_Index = paste0("14_",mrna14$Cell_Index)
mrna.all = rbind(mrna1, mrna2, mrna3, mrna4, mrna5, mrna6, mrna7, mrna8, mrna9, mrna10, mrna14)

mrna15$Cell_Index = paste0("15_",mrna15$Cell_Index)
mrna16$Cell_Index = paste0("16_",mrna16$Cell_Index)
mrna.all = rbind(mrna15, mrna16)

mrna2019.1$Cell_Index = paste0("11_",mrna2019.1$Cell_Index)
mrna2019.2$Cell_Index = paste0("12_",mrna2019.2$Cell_Index)
mrna3hrXtra$Cell_Index = paste0("13_", mrna3hrXtra$Cell_Index)
mrna2019.all = rbind(mrna2019.1, mrna2019.2, mrna3hrXtra)

mrna.all = rbind(mrna.all, mrna2019.all)

samptag1$Cell_Index = paste0("1_",samptag1$Cell_Index)
samptag2$Cell_Index = paste0("2_",samptag2$Cell_Index)
samptag3$Cell_Index = paste0("3_",samptag3$Cell_Index)
samptag4$Cell_Index = paste0("4_",samptag4$Cell_Index)
samptag5$Cell_Index = paste0("5_",samptag5$Cell_Index)
samptag6$Cell_Index = paste0("6_",samptag6$Cell_Index)
samptag7$Cell_Index = paste0("7_",samptag7$Cell_Index)
samptag8$Cell_Index = paste0("8_",samptag8$Cell_Index)
samptag9$Cell_Index = paste0("9_",samptag9$Cell_Index)
samptag10$Cell_Index = paste0("10_",samptag10$Cell_Index)
samptag14$Cell_Index = paste0("14_",samptag14$Cell_Index)

samptag15$Cell_Index = paste0("15_",samptag15$Cell_Index) #PMs1
samptag16$Cell_Index = paste0("16_",samptag16$Cell_Index) #PMs2

samptag2019.1$Cell_Index = paste0("11_",samptag2019.1$Cell_Index)
samptag2019.2$Cell_Index = paste0("12_",samptag2019.2$Cell_Index)
samptag3hrXtra$Cell_Index = paste0("13_",samptag3hrXtra$Cell_Index)

samptag1$Sample_Name = paste0("1_",samptag1$Sample_Name)
samptag2$Sample_Name = paste0("2_",samptag2$Sample_Name)
samptag3$Sample_Name = paste0("3_",samptag3$Sample_Name)
samptag4$Sample_Name = paste0("4_",samptag4$Sample_Name)
samptag5$Sample_Name = paste0("5_",samptag5$Sample_Name)
samptag6$Sample_Name = paste0("6_",samptag6$Sample_Name)
samptag7$Sample_Name = paste0("7_",samptag7$Sample_Name)
samptag8$Sample_Name = paste0("8_",samptag8$Sample_Name)
samptag9$Sample_Name = paste0("9_",samptag9$Sample_Name)
samptag10$Sample_Name = paste0("10_",samptag10$Sample_Name)
samptag14$Sample_Name = paste0("14_",samptag14$Sample_Name)

samptag15$Sample_Name = paste0("15_",samptag15$Sample_Name)
samptag16$Sample_Name = paste0("16_",samptag16$Sample_Name)
samptag.all = rbind(samptag1, samptag2, samptag3, samptag4, samptag5, samptag6, samptag7, samptag8, samptag9, samptag10, samptag14)


# 2019 samples 500g genes ----------------------------------------------------------------
samptag2019.1$Sample_Name = paste0("11_",samptag2019.1$Sample_Name)
samptag2019.2$Sample_Name = paste0("12_",samptag2019.2$Sample_Name)
samptag3hrXtra$Sample_Name = paste0("13_",samptag3hrXtra$Sample_Name)
samptag2019.all = rbind(samptag2019.1, samptag2019.2, samptag3hrXtra)

samptag2019.all$timept = ifelse(grepl("11_", samptag2019.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag2019.all$Sample_Tag), "8hr",
                         ifelse(grepl("11_", samptag2019.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag2019.all$Sample_Tag), "0.5hr",
                                ifelse(grepl("12_", samptag2019.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag2019.all$Sample_Tag), "5hr",
                                       ifelse(grepl("12_", samptag2019.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag2019.all$Sample_Tag), "3hr", 
                                              ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("01|02|03|04|05", samptag2019.all$Sample_Name), "3hr", 
                                                     ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("06", samptag2019.all$Sample_Name), "0hr", NA
                                       ))))))
samptag2019.all$stimulus = ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("01|07", samptag2019.all$Sample_Tag), "LPS",
                                  ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("02|08", samptag2019.all$Sample_Tag), "TNF",
                                         ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("03|09", samptag2019.all$Sample_Tag), "P3CSK",
                                                ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("04|10", samptag2019.all$Sample_Tag), "CpG",
                                                       ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("05|11", samptag2019.all$Sample_Tag), "PIC",
                                                              ifelse(grepl("11_|12_", samptag2019.all$Cell_Index) & grepl("06|12", samptag2019.all$Sample_Tag), "LPSlo",
                                                                     ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("01|02", samptag2019.all$Sample_Name), "LPS",
                                                                            ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("03", samptag2019.all$Sample_Name), "TNF",
                                                                                   ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("04", samptag2019.all$Sample_Name), "CpG",
                                                                                          ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("05", samptag2019.all$Sample_Name), "PIC",
                                                                                                 ifelse(grepl("13_", samptag2019.all$Cell_Index) & grepl("06", samptag2019.all$Sample_Name), "Unstim", NA
                                                              )))))))))))
samptag2019.all$type = "M0"
samptag2019.all$batch = gsub("_..*", "", samptag2019.all$Cell_Index)
samptag2019.all$replicate = "rep1"                                                     
#--------------------------------------------------------------                                         

samptag.all$timept = ifelse(grepl("^1_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "5hr",
                            ifelse(grepl("^1_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "3hr",
                                   ifelse(grepl("^2_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "8hr",
                                          ifelse(grepl("^2_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "0.5hr",
                                                 ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "1hr",
                                                        ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "0hr",
                                                               ifelse(grepl("^4_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "24hr",
                                                                      ifelse(grepl("^4_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "1hr",
                                                                             ifelse(grepl("5_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "5hr", 
                                                                                    ifelse(grepl("5_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "3hr", 
                                                                                           ifelse(grepl("6_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "8hr",
                                                                                                  ifelse(grepl("6_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "0.5hr",
                                                                                                         ifelse(grepl("7_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "8hr",
                                                                                                                ifelse(grepl("7_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "0.25hr",
                                                                                                                       ifelse(grepl("8_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "3hr",
                                                                                                                              ifelse(grepl("8_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "1hr",
                                                                                                                                     ifelse(grepl("9_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "8hr",
                                                                                                                                            ifelse(grepl("9_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "8hr",
                                                                                                                                                   ifelse(grepl("10_|14_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "8hr",
                                                                                                                                                          ifelse(grepl("10_|14_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "8hr",NA
                                                                      ))))))))))))))))))))
samptag.all$stimulus = ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("01|07", samptag.all$Sample_Tag), "LPS",
                              ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("02|08", samptag.all$Sample_Tag), "TNF",
                                     ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("03|09", samptag.all$Sample_Tag), "P3CSK",
                                            ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("04|10", samptag.all$Sample_Tag), "CpG",
                                                   ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("05|11", samptag.all$Sample_Tag), "PIC",
                                                          ifelse(grepl("^1_|^2_|^4_|5_|6_|7_|8_", samptag.all$Cell_Index) & grepl("06|12", samptag.all$Sample_Tag), "IFNb",
                                                                 ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("01|07", samptag.all$Sample_Tag), "LPS",
                                                                        ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("02|08", samptag.all$Sample_Tag), "TNF",
                                                                               ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("03|09", samptag.all$Sample_Tag), "P3CSK",
                                                                                      ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("04|10", samptag.all$Sample_Tag), "PIC",
                                                                                             ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("05|11", samptag.all$Sample_Tag), "IFNb",
                                                                                                    ifelse(grepl("9_|10_|14_", samptag.all$Cell_Index) & grepl("06|12", samptag.all$Sample_Tag), "Unstim",
                                                                                                           ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("01", samptag.all$Sample_Tag), "LPS",
                                                                                                                  ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("01|02", samptag.all$Sample_Tag), "TNF",
                                                                                                                         ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("03", samptag.all$Sample_Tag), "P3CSK",
                                                                                                                                ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("04", samptag.all$Sample_Tag), "CpG",
                                                                                                                                       ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("05", samptag.all$Sample_Tag), "PIC",
                                                                                                                                              ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("06", samptag.all$Sample_Tag), "IFNb",
                                                                                                                  ifelse(samptag.all$timept=="0hr","Unstim",NA)
                                                                                                                  ))))))))))))))))))

samptag.all$type = ifelse(grepl("^1_|^2_", samptag.all$Cell_Index) , "M2_IL4",
                          ifelse(grepl("^4_|5_|6_", samptag.all$Cell_Index) , "M1_IFNg",
                                 ifelse(grepl("7_|8_", samptag.all$Cell_Index) , "M0",
                                        ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06|07|08", samptag.all$Sample_Tag), "M2_IL4",
                                               ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("09|10", samptag.all$Sample_Tag), "M1_IFNg",
                                                      ifelse(grepl("^3_", samptag.all$Cell_Index) & grepl("11|12", samptag.all$Sample_Tag), "M0",
                                                             ifelse(grepl("9_", samptag.all$Cell_Index) & grepl("01|02|03|04|05|06", samptag.all$Sample_Tag), "BMDM1_Healthy",
                                                                    ifelse(grepl("9_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "BMDM1_Sick",
                                                                          ifelse(grepl("10_", samptag.all$Cell_Index)& grepl("01|02|03|04|05|06", samptag.all$Sample_Tag) , "BMDM2_WT",
                                                                                ifelse(grepl("10_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "BMDM2_MM",
                                                                                       ifelse(grepl("14_", samptag.all$Cell_Index)& grepl("01|02|03|04|05|06", samptag.all$Sample_Tag) , "PM_B6.WT",
                                                                                              ifelse(grepl("14_", samptag.all$Cell_Index) & grepl("07|08|09|10|11|12", samptag.all$Sample_Tag), "PM_B6.AireGW",
                                                                                       
                                                             NA
                                                      ))))))))))))
                              
samptag.all$batch = gsub("_..*", "", samptag.all$Cell_Index)
samptag.all$replicate = ifelse(grepl("0hr", samptag.all$timept)&grepl("08|10|12", samptag.all$Sample_Tag), "rep2", #0hrs done as replicates
                               ifelse(grepl("7_|8_", samptag.all$Cell_Index), "rep2", "rep1")) #replicates of M0


# annotate PMs Feb 2021--------------------------------------------------------------
samptag2021 = rbind(samptag15, samptag16)
samptag2021$timept = "3hr"
samptag2021$stimulus = ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("01|07", samptag2021$Sample_Tag), "LPS",
                                  ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("02|08", samptag2021$Sample_Tag), "TNF",
                                         ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("03|09", samptag2021$Sample_Tag), "PIC",
                                                ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("04|10", samptag2021$Sample_Tag), "P3CSK",
                                                       ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("05|11", samptag2021$Sample_Tag), "IFNb",
                                                              ifelse(grepl("15_|16_", samptag2021$Cell_Index) & grepl("06|12", samptag2021$Sample_Tag), "Unstim", NA
                                                                                                 ))))))
samptag2021$type = ifelse(grepl("^15_", samptag2021$Cell_Index) & grepl("01|02|03|04|05|06", samptag2021$Sample_Tag),"PM_B6.HFD",
                          ifelse(grepl("^15_", samptag2021$Cell_Index)& grepl("07|08|09|10|11|12", samptag2021$Sample_Tag) , "PM_B6.LFD",
                                 ifelse(grepl("^16_", samptag2021$Cell_Index) , "PM_B6.old",NA
                                 )))
                                        
samptag2021$batch = gsub("_..*", "", samptag2021$Cell_Index)
samptag2021$replicate = "rep1"   




samptag.all = rbind(samptag.all, samptag2019.all)
# write.table(samptag.all, "./output/samptag.all_cellannotations_metadata.txt",quote = F, sep = "\t",row.names = F)
samptag.all = read.delim("./output/samptag.all_cellannotations_metadata.txt")
samptag.all = rbind(samptag.all, samptag2021)
# write.table(samptag.all, "./output/samptag.all_cellannotations_metadata.txt",quote = F, sep = "\t",row.names = F)

#downsample
# set.seed(1)
# mrna1= mrna1[sample(nrow(mrna1),12000),]
# mrna2= mrna2[sample(nrow(mrna2),12000),]
# mrna3= mrna3[sample(nrow(mrna3),5000),]
# mrna4= mrna4[sample(nrow(mrna4),7000),]
# mrna5= mrna5[sample(nrow(mrna5),12000),]
# mrna.all = rbind(mrna1, mrna2, mrna3, mrna4, mrna5)

#remove Cell_indexes that are multiplet of undetermined-----------------------
na.omit(table(samptag.all$timept, samptag.all$stimulus, samptag.all$type))
na.omit(table(samptag.all$batch, samptag.all$stimulus))
na.omit(table(samptag.all$type, samptag.all$stimulus))
na.omit(table(samptag.all$batch, samptag.all$type))
remove= samptag.all$Cell_Index[is.na(samptag.all$timept)]

# clean data ------------------------------------------------------------------------------------------
head(mrna.all)[1:5]
rownames(mrna.all) = mrna.all$Cell_Index
data = mrna.all[,-1]

if(0){ #for 2019 samples
  na.omit(table(samptag2019.all$batch, samptag2019.all$stimulus))
  na.omit(table(samptag2019.all$timept, samptag2019.all$stimulus, samptag2019.all$batch))
  
  remove = samptag2019.all$Cell_Index[is.na(samptag2019.all$timept)]
  rownames(mrna2019.all) = mrna2019.all$Cell_Index
  data = mrna2019.all[,-1]
  data = data[!rownames(data) %in% remove,]
  data = data[ ,colSums((data == 0)) < nrow(data) ] #take out genes that are 0 across ~70K cells (40 of them)
  head(data)[1:5]
}


# gene_panel = gsub("_ENSMUST..*","", colnames(mrna.all)[-1])
# gene_panel = gsub("_NM_..*","", gene_panel)
# write.csv(gene_panel, "F://scRNAseq_macro/bulk_rnaseq/gene_panel_mouse.csv", quote = F, row.names = F)

# remove cells from data
data = data[!rownames(data) %in% remove,]
data = data[ ,colSums((data == 0)) < nrow(data) ] #take out genes that are 0 across ~70K cells (40 of them)
head(data)[1:5]

# subset the data by type----------------------------------------------------------------------
base_samps = samptag.all$Cell_Index[(samptag.all$timept=="0hr"|samptag.all$timept=="0.25hr")]
base_samps_0hr = samptag.all$Cell_Index[(samptag.all$timept=="0hr")]
M0samps = samptag.all$Cell_Index[(samptag.all$type=="M0")]
M1samps = samptag.all$Cell_Index[(samptag.all$type=="M1_IFNg")]
M2samps = samptag.all$Cell_Index[(samptag.all$type=="M2_IL4")]
BMDM1samps = samptag.all$Cell_Index[(samptag.all$type=="BMDM1_Healthy"|samptag.all$type=="BMDM1_Sick")]
BMDM2samps = samptag.all$Cell_Index[(samptag.all$type=="BMDM2_WT"|samptag.all$type=="BMDM2_MM")]

BMDMsamps = samptag.all$Cell_Index[(grepl("BMDM", samptag.all$type))]
M0samps_all = samptag.all$Cell_Index[(samptag.all$type=="M0")]

BMDM1.Hsamps = samptag.all$Cell_Index[(samptag.all$type=="BMDM1_Healthy")]
BMDM1.Ssamps = samptag.all$Cell_Index[(samptag.all$type=="BMDM1_Sick")]
BMDM2.WTsamps = samptag.all$Cell_Index[(samptag.all$type=="BMDM2_WT")]
BMDM2.MMsamps = samptag.all$Cell_Index[(samptag.all$type=="BMDM2_MM")]

data.base = data[ rownames(data) %in% base_samps,]
data.base0 = data[ rownames(data) %in% base_samps_0hr,]
data.M0 = data[ rownames(data) %in% M0samps,]
data.M1 = data[ rownames(data) %in% M1samps,]
data.M2 = data[ rownames(data) %in% M2samps,]
data.BMDM1 = data[ rownames(data) %in% BMDM1samps,]
data.BMDM2 = data[ rownames(data) %in% BMDM2samps,]

data.BMDM = data[ rownames(data) %in% BMDMsamps,]
data.M0all = data[ rownames(data) %in% M0samps_all,]

data.BMDM1.H = data[ rownames(data) %in% BMDM1.Hsamps,]
data.BMDM1.S = data[ rownames(data) %in% BMDM1.Ssamps,]
data.BMDM2.WT = data[ rownames(data) %in% BMDM2.WTsamps,]
data.BMDM2.MM = data[ rownames(data) %in% BMDM2.MMsamps,]


#summary of downsampled
# samptag.downsampled = samptag.all[samptag.all$Cell_Index %in% rownames(data),]
# (table(samptag.downsampled$timept, samptag.downsampled$stimulus))
# na.omit(table(samptag.downsampled$batch, samptag.downsampled$stimulus))
# na.omit(table(samptag.downsampled$batch, samptag.downsampled$timept))

#assess 0 distribution
res <- colSums(data ==0)/nrow(data)#get prop of 0 in any gene
hist(res, breaks = 30)
res <- rowSums(data ==0)/ncol(data )#get prop of 0 in any cell
hist(res, breaks = 30)
# # do PCA+=---------------------------------------------------------------
# 
# pca = prcomp(data, center = T, scale = T)
# pca_scores = pca$x
# pca_scores = cbind(Score = gsub("-", ".", rownames(pca_scores)),pca_scores)
# pca_loadings = pca$rotation
# pca_loadings = cbind(Loading = colnames(data),  pca_loadings)
# pca_evalues = pca$sdev
# 
# # name = "rhapsody_2019all"
# name = "rhapsodyy_L6"
# savename = paste0(out_dir,name, "_prcomp_scores.txt")
# write.table(pca_scores, savename, sep = "\t", row.names = FALSE, 
#             quote = FALSE)
# savename = paste0(out_dir,name, "_prcomp_loadings.txt")
# write.table(pca_loadings, savename, sep = "\t", row.names = FALSE, 
#             quote = FALSE)
# savename = paste0(out_dir,name, "_prcomp_sdev.txt")
# write.table(pca_evalues, savename, sep = "\t", row.names = FALSE, 
#             quote = FALSE)
# print(summary(pca))
# screeplot(pca)
# 
# # plot pca ---------------------------------------------------------
# 
# scores = read.delim("analysis_rhapsody/rhapsody_2019all_prcomp_scores.txt")
# scores = read.delim("analysis_rhapsody/rhapsody_L6_prcomp_scores.txt")
# plot_pca("analysis_rhapsody/rhapsody_2019all_prcomp_scores.txt", samptag.all$Cell_Index, samptag.all$Sample_Name)
# 
# file = "analysis_rhapsody/rhapsody_2019all_prcomp_scores.txt"
# PCx = "PC1"; PCy = "PC2"; labels = T; 
# info.name = samptag.all$Cell_Index; info.type = samptag.all$Sample_Name
# table = scores
# table$type = info.type[match(table$Score, info.name)]
# sdev = read.delim(paste0(gsub("scores.txt", "", file), "sdev.txt"))
# sdev$var = unlist(sdev^2)
# sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
# rownames(sdev) = paste0("PC", seq(1, nrow(sdev)))
# pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + geom_point(size = I(3), 
#                                                                   aes(color = factor(type))) + theme(legend.position = "right", 
#                                                                                                      plot.title = element_text(size = 30), legend.text = element_text(size = 22), 
#                                                                                                      legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
#                                                                                                      legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)), 
#                                                                                                      axis.text.y = element_text(margin = margin(l = -14))) + 
#   guides(color = guide_legend(title = "Type")) + labs(
#                                                       x = paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], 
#                                                                  "%)"), y = paste0(PCy, " (", sdev$pve[match(PCy, 
#                                                                                                              rownames(sdev))], "%)")) + theme_bw(base_size = 18) + 
#   if (labels == TRUE) {
#     geom_text(data = table, mapping = aes(label = Score), 
#               check_overlap = TRUE, size = 3)
#   }
# 
# pcx.y
# 
# # do tnse-----------------------------------------------------------------------
# rownames(scores) = scores$Score
# rownames(samptag.all) = samptag.all$Cell_Index
# tsne= tsne(t(scores[, c(2:11)]),labels=as.factor(samptag.all$Sample_Name))
# tsne
# # tsne$plot_env$labels = as.factor(samptag.all$Sample_Name)
# tsne$plot_env$dotsize = 1
# tsne
# 
# 
# 



# tempora package----
library(devtools)
install_github("BaderLab/Tempora")
library(Tempora)
