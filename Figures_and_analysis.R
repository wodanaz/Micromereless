# ------------------------------------------------------------------------------
# MS ID#: DEVELOP/2024/203152
# MS TITLE: Reprogramming of cells during embryonic transfating: overcoming a reprogramming block
# AUTHORS: Alejandro Berrio, Esther Miranda, Abdull J Massri, Anton Afanassiev,
#          Geoffrey Schiebinger, Gregory A Wray, and David R. McClay
# ------------------------------------------------------------------------------
# This R script contains the code and analysis used in the above manuscript.
# Created by: Alejandro Berrio
# Date: 11/06/2024
#
# Usage:
# - Ensure all required packages are installed.
# - Follow the instructions provided in the following file for running the scripts.
#
# License:
# This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License.
# To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to
# Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
#
# Contact:
# For any questions or issues, please contact alejo.berrio@duke.edu
# ------------------------------------------------------------------------------


library(Seurat) # Seurat Mapping was done in HARDAC (Duke's HPC), using Seurat version 4.1.2
library(gplots)
library(viridis)
library(xlsx)
library(Matrix)
library(tidyverse)
library(scCustomize)
library(cowplot)
library(CIDER)
library(pheatmap)

load(file="precomputed_control.Rda", verbose = T)
load(file="precomputed_micromereless.Rda", verbose = T)
load(file="precomputed_integrated.Rda", verbose = T)
pal <- viridis(n = 10, option = "D")

#DimPlot_scCustom(control_6to18,  label=F, group.by="Stage") +  theme_void() + NoLegend()
#DimPlot_scCustom(control_6to18,  label=T, group.by="seurat_clusters")  +  theme_void() + NoLegend()
#DimPlot_scCustom(control_6to18,  label=T, group.by="celltype") + theme_void() + NoLegend()

#DimPlot_scCustom(micromereless_6to18,  label=F, group.by="Stage") +  theme_void() + NoLegend()
#DimPlot_scCustom(micromereless_6to18,  label=T, group.by="seurat_clusters")  +  theme_void() + NoLegend()
#DimPlot_scCustom(micromereless_6to18,  label=T, group.by="celltype") + theme_void() + NoLegend()

#DimPlot_scCustom(mmless_integrated_sct,  label=F, group.by="Stage") +  theme_void() + NoLegend()
#DimPlot_scCustom(mmless_integrated_sct,  label=T, group.by="seurat_clusters")  +  theme_void() + NoLegend()
#DimPlot_scCustom(mmless_integrated_sct,  label=T, group.by="celltype") + theme_void() + NoLegend()


##########################################
### Figure 2
##########################################
# Figure 2.A Control 6-18 hpf by celltype

control  <- as.data.frame(Embeddings(control_6to18, "umap"))
head(control)
dim((control))

control$cells <- rownames(control)
control
colnames(control)

metadata <- control_6to18@meta.data 
metadata$cells <- rownames(metadata)

control2 <- merge(control, metadata, by = "cells", all.x = T)
head(control2)

control3 <- control2 %>% separate(Stage, c("time", "hpf"), sep = "_", remove = F)
control3$time <- as.numeric(control3$time)


# Compute the average value of time for each unique value in celltype
control3 %>% group_by(celltype) %>% summarize(avg_time = mean(time)) %>% 
  write_csv(file = sprintf("control/control_celltype_avg_times_%s.csv", Sys.Date()))

control3 %>% group_by(seurat_clusters) %>% summarize(avg_time = mean(time)) %>% 
  write_csv(file = sprintf("control/control_cluster_avg_times_%s.csv", Sys.Date()))



unique(control2$celltype)

ggplot() + 
  geom_point(data=subset(control2,  celltype == "Undefined"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="turquoise2",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "APD"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="midnightblue",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Neural"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="blueviolet",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Aboral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue4",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Oral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="dodgerblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Ciliary_Band"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deepskyblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Early_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightskyblue",  pch=20, alpha=1) + 
  
  geom_point(data=subset(control2,  celltype == "Gut"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#d45500",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff6600",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Early_Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff9955",  pch=20, alpha=1) + 
  
  
  geom_point(data=subset(control2,  celltype == "Blastocoelar"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#FFFF00",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Pigment"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff0000",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Mesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#877fb3",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Endomesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#c5795eff",  pch=20, alpha=1) + 
  
  geom_point(data=subset(control2,  celltype == "PGCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deeppink3",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "PMCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#00ff00",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  celltype == "Micromeres"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()

##############################################
# Figure 2.C Control 6-18 hpf by stage
ggplot() + 
  geom_point(data=subset(control2,  Stage == "6_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="#fde725",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "8_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#fdb42f",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "10_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ed7953",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "12_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#cc4778",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "14_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#9c179e",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "16_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#5c01a6",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  Stage == "18_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#0d0887",  pch=20, alpha=1) + 
  theme_bw() + 
  
  ggtitle("") +  theme_void()

################################################
# Figure 2.B Micromereless 6-18 hpf by celltype

mmless  <- as.data.frame(Embeddings(micromereless_6to18, "umap"))
head(mmless)
colnames(mmless)
mmless$cells <- rownames(mmless)
mmless
colnames(mmless)

metadata <- micromereless_6to18@meta.data 
metadata$cells <- rownames(metadata)

mmless2 <- merge(mmless, metadata, by = "cells", all.x = T)
head(mmless2)
dim(mmless2)


mmless3 <- mmless2 %>% separate(Stage, c("time", "mmless"), sep = "_", remove = F)
mmless3$time <- as.numeric(mmless3$time)


# Compute the average value of time for each unique value in celltype
mmless3 %>% group_by(celltype) %>% summarize(avg_time = mean(time)) %>% 
  write_csv(file = sprintf("mmless/mmless_celltype_avg_times_%s.csv", Sys.Date()))

mmless3 %>% group_by(seurat_clusters) %>% summarize(avg_time = mean(time)) %>% 
  write_csv(file = sprintf("mmless/mmless_cluster_avg_times_%s.csv", Sys.Date()))


ggplot() + 
  geom_point(data=subset(mmless2,  celltype == "Undefined"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="turquoise2",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "APD"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="midnightblue",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Neural"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="blueviolet",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Aboral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue4",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Oral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="dodgerblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Ciliary_Band"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deepskyblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Early_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightskyblue",  pch=20, alpha=1) + 
  
  geom_point(data=subset(mmless2,  celltype == "Gut"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#d45500",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff6600",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Early_Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff9955",  pch=20, alpha=1) + 
  
  
  geom_point(data=subset(mmless2,  celltype == "Blastocoelar"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#FFFF00",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Pigment"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff0000",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Mesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#877fb3",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Endomesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#c5795eff",  pch=20, alpha=1) + 
  
  geom_point(data=subset(mmless2,  celltype == "PGCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deeppink3",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "PMCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#00ff00",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  celltype == "Micromeres"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()


##############################################
# Figure 2.D Micromereless 6-18 hpf by stage
ggplot() + 
  geom_point(data=subset(mmless2,  Stage == "6_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="#fde725",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "8_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#fdb42f",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "10_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ed7953",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "12_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#cc4778",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "14_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#9c179e",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "16_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#5c01a6",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  Stage == "18_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#0d0887",  pch=20, alpha=1) + 
  theme_bw() + 
  
  ggtitle("") +  theme_void()


################################################
### Figure 3



integration  <- as.data.frame(Embeddings(mmless_integrated_sct, "umap"))
head(integration)
colnames(integration)
integration$cells <- rownames(integration)
integration
colnames(integration)

metadata <- mmless_integrated_sct@meta.data 
metadata$cells <- rownames(metadata)

integration2 <- merge(integration, metadata, by = "cells", all.x = T)
head(integration2)

unique(integration2$celltype)

################################################
# Figure 3.A Integration 6-18 hpf by experiment
ggplot() + 
  geom_point(data=subset(integration2, experiment == "CONTROL" ),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="black",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "MMLESS" ),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()

################################################
# Figure 3.B Integration 6-18 hpf Control Alone
ggplot() + 
  geom_point(data=subset(integration2, experiment == "MMLESS" ),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="white",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" ),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="black",  pch=20, alpha=1) + 
  theme_bw() + 
  ggtitle("") +  theme_void()

################################################
# Figure 3.C Integration 6-18 hpf Mmless Alone
ggplot() + 
  geom_point(data=subset(integration2, experiment == "CONTROL" ),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="white",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "MMLESS" ),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()

################################################
# Figure 3.D Integration 6-18 hpf by celltype

ggplot() + 
  #geom_point(data=integration2,  aes(x=(UMAP_1), y=(UMAP_2)), color="gray85", pch=20, alpha=0.1) +
  geom_point(data=subset(integration2,  celltype2 == "Undefined"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="turquoise2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "APD"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="midnightblue",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Neural"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="blueviolet",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Aboral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue4",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Oral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="dodgerblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Ciliary_Band"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deepskyblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Early_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightskyblue",  pch=20, alpha=1) + 
  
  geom_point(data=subset(integration2,  celltype2 == "Gut"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orange",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightgoldenrod",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Early_Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gold",  pch=20, alpha=1) + 
  
  
  geom_point(data=subset(integration2,  celltype2 == "Blastocoelar"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="springgreen4",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Pigment"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="plum",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Mesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="coral",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Endomesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="navajowhite4",  pch=20, alpha=1) + 
  
  geom_point(data=subset(integration2,  celltype2 == "PGCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deeppink3",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "PMCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="firebrick1",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Micromeres"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()



ggplot() + 
  geom_point(data=subset(integration2,  celltype2 == "Undefined"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="turquoise2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "APD"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="midnightblue",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Neural"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="blueviolet",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Aboral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue4",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Oral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype == "Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="dodgerblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Ciliary_Band"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deepskyblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Early_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightskyblue",  pch=20, alpha=1) + 
  
  geom_point(data=subset(integration2,  celltype2 == "Gut"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#d45500",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff6600",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Early_Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff9955",  pch=20, alpha=1) + 
  
  
  geom_point(data=subset(integration2,  celltype2 == "Blastocoelar"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#FFFF00",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Pigment"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ff0000",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Mesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#877fb3",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Endomesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#c5795eff",  pch=20, alpha=1) + 
  
  geom_point(data=subset(integration2,  celltype2 == "PGCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deeppink3",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "PMCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#00ff00",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  celltype2 == "Micromeres"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()



# Split control and mmless
head(integration2)

ggplot() + 
  #geom_point(data=integration2,  aes(x=(UMAP_1), y=(UMAP_2)), color="gray85", pch=20, alpha=0.1) +
  geom_point(data=subset(integration2, experiment == "CONTROL" &  celltype2 == "Undefined"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="turquoise2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" &  celltype2 == "APD"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="midnightblue",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Neural"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="blueviolet",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Aboral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue4",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Oral_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="royalblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="dodgerblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Ciliary_Band"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deepskyblue2",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Early_Ectoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightskyblue",  pch=20, alpha=1) + 
  
  #geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Gut"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orange",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" &  celltype2 == "Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="lightgoldenrod",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Early_Endoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gold",  pch=20, alpha=1) + 
  
  
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Blastocoelar"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="springgreen4",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Pigment"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="plum",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "Mesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="coral",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" &  celltype2 == "Endomesoderm"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="navajowhite4",  pch=20, alpha=1) + 
  
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "PGCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="deeppink3",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" & celltype2 == "PMCs"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="firebrick1",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2, experiment == "CONTROL" &  celltype2 == "Micromeres"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="orchid1",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()


ggplot() + 
  #geom_point(data=hpf,  aes(x=(UMAP_1), y=(UMAP_2)), color="gray85", pch=20, alpha=0.1) +
  geom_point(data=subset(integration2,  Stage == "6_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="#fde725",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "8_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#fdb42f",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "10_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#ed7953",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "12_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#cc4778",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "14_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#9c179e",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "16_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#5c01a6",  pch=20, alpha=1) + 
  geom_point(data=subset(integration2,  Stage == "18_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="#0d0887",  pch=20, alpha=1) + 
  theme_bw() + 
  
  ggtitle("") +  theme_void()



################################################
### Figure 4

DimPlot(micromereless_6to18,  label=T, group.by="seurat_clusters")  +  NoLegend()
DimPlot_scCustom(micromereless_6to18,  label=T, group.by="seurat_clusters")  +  theme_void() + NoLegend()


#######################################################
### Figure 5: The following scripts demonstrate how we 
###           made our box figures from our collected
###           data and statistic analyses

### FIGURE 5B
# Control vs Micromereless 
library(tidyverse)
library(hrbrthemes)
library(viridis)
# create a dataset
data <- data.frame(
  experiment=c( #rep("Control",3), 
    rep("Control",12), 
    rep("Mmless",30)),
  celltype=c( #rep("PMCs",1), rep("Blasto",1), rep("Pigment",1),
    rep("PMCs",4), rep("Blasto",4), rep("Pigment",4),
    rep("PMCs",10), rep("Blasto",10), rep("Pigment",10)),
  value=c( #c(70,  61,	71),           
    c(57,	36,	54,	57,	65,	52,	40,	50,	58,	43,	68,	51),
    c(71,	61,	48,	73,	57,	57,	52,	74,	60,	51,	4,	6,	18,	18,	1,	11,	0,	0,	0,	4,	0,	0,	0,	25,	1,	0,	0,	0,	0,	0))
)



my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")
experiment_order <- c("Control", "Mmless")
data <- transform(data, celltype = factor(celltype, levels = cell_order))
Control <- data[data$experiment %in% "Control", ]
Mmless  <-  data[data$experiment %in% "Mmless", ]


# Group by cell type and perform Kruskal-Wallis tests
results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment)$p.value
  )

as.data.frame(results)


ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Control vs Micromereless",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ factor(celltype, levels = cell_order), scales = "free_x", strip.position = "top") +  # Move facet titles to top
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) +   ylim(0, 100) +
  geom_text(data = as.data.frame(results), aes(x = 2.2, y = 100, label = paste("p =", round(p_value, 3)))) 


### FIGURE 6B
# Control vs Micromereless  with Brachyuri MO


data <- data.frame(
  experiment=c( rep("Control",9), 
                rep("Control Bra KD",15), 
                rep("Mmless Bra KD",12)#,
                #rep("Mmless Bra Lo KD",12)
  ),
  
  celltype=c( rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",5), rep("Blasto",5), rep("Pigment",5),
              rep("PMCs",4), rep("Blasto",4), rep("Pigment",4)#,
              #rep("PMCs",4), rep("Blasto",4), rep("Pigment",4)
  ),
  value=c( c(62,	52,	60, 60,	57,	56, 82,	52,	60 ),           
           c(53,	63,	52,	68,	81, 96,	38,	48,	35,	90, 57,	33,	45,	33,	60 ),
           c(1,	1,	0,	0, 1,	1,	2,	0, 0,	0,	0,	0)#,
           #c(41,	18,	44,	15, 6,	0,	3,	2, 0,	0,	0,	0)
  )
)




# Create a custom order for cell types
experiment_order <- c("Control", "Control Bra KD", "Mmless Bra KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))



results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )


subset_data1 <- subset(data, experiment %in% c('Control', 'Control Bra KD'))
subset_data2 <- subset(data, experiment %in% c('Control Bra KD', 'Mmless Bra KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Bra KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )



results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )



# Convert the results to a data frame
results_df <- as.data.frame(results)

# Print the results
print(results_df)



# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Bra KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 100)



results_sd1
results_sd2
results_sd3


### FIGURE 7C
# Control vs Micromereless  with Delta MO

data <- data.frame(
  experiment=c( rep("Control",9), 
                rep("Control Delta KD",9), 
                rep("Mmless Delta KD",15)),
  celltype=c( rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",5), rep("Blasto",5), rep("Pigment",5)),
  value=c( c(68,	64,	55 ,62 ,	59,	38, 63 ,	81,	67),           
           c(53,	65,	61, 43,	43,	52, 5,	6,	12),
           c(209,	165,	172,	114,	197,19,	46,	9	,2,	30,1,	12,	4,	0,	4)))


# Create a custom order for cell types
experiment_order <- c("Control", "Control Delta KD", "Mmless Delta KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Delta KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 210)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

subset_data1 <- subset(data, experiment %in% c('Control', 'Control Delta KD'))
subset_data2 <- subset(data, experiment %in% c('Control Delta KD', 'Mmless Delta KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Delta KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1
results_sd2
results_sd3

### FIGURE 7E ############################
# Control vs Micromereless  with Nodal MO

data <- data.frame(
  experiment=c( rep("Control",15), 
                rep("Control Nodal KD",9), 
                rep("Mmless Nodal KD",18)),
  celltype=c( rep("PMCs",5), rep("Pigment",5), rep("Blasto",5),
              rep("PMCs",3), rep("Pigment",3), rep("Blasto",3),
              rep("PMCs",6), rep("Pigment",6), rep("Blasto",6)),
  value=c( c(60,	57,	56,	53,	63,
             52,	40,	29,	113,	68,
             126,	41,	29,	63,	42),           
           c(67,	61,	84,
             84,	68,	81,
             18,	39,	43),
           c(123,	112,	118,	120,	77,	111,
             15,	19,	21,	34,	84,	38,
             23,	12,	8,	11,	23,	6))
)


experiment_order <- c("Control", "Control Nodal KD", "Mmless Nodal KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Nodal KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 150)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

subset_data1 <- subset(data, experiment %in% c('Control', 'Control Nodal KD'))
subset_data2 <- subset(data, experiment %in% c('Control Nodal KD', 'Mmless Nodal KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Nodal KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )


results_sd1
results_sd2
results_sd3

### FIGURE 7G
# Control vs Micromereless  with Lefty MO

data <- data.frame(
  experiment=c( rep("Control",9), 
                rep("Control Lefty KD",6), 
                rep("Mmless Lefty KD",21)),
  celltype=c( rep("PMCs",3), rep("Pigment",3), rep("Blasto",3),
              rep("PMCs",2), rep("Pigment",2), rep("Blasto",2),
              rep("PMCs",7), rep("Pigment",7), rep("Blasto",7)),
  value=c( c(50,	46,	41,
             43,	74,	59,
             43,	49,	44),           
           c(55,	65,	
             0,	4,
             65,	58 ),
           c(67,	41,	29,	36,	48,	34,	29,
             0,	0,	0,	0,	0,	0,	0,
             1,	0,	0,	0,	0,	3,	0))
)

# Create a custom order for cell types
experiment_order <- c("Control", "Control Lefty KD", "Mmless Lefty KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Lefty KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 100)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

subset_data1 <- subset(data, experiment %in% c('Control', 'Control Lefty KD'))
subset_data2 <- subset(data, experiment %in% c('Control Lefty KD', 'Mmless Lefty KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Lefty KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1
results_sd2
results_sd3




### FIGURE 7J
# Control vs Micromereless  with Delta Over-Expression

# create a dataset
data <- data.frame(
  experiment=c( 
    rep("Control",10), 
    rep("Mmless Delta+++",20),
    rep("Mmless Control",18),
    rep("Mmless Delta 6hr",16),
    rep("Mmless Delta 8hr",8),
    rep("Mmless Delta 10hr",12)),
  celltype=c( 
    rep("Blasto",5),  rep("Pigment",5),
    rep("Blasto",10),  rep("Pigment",10),
    rep("Blasto",9),  rep("Pigment",9),
    rep("Blasto",8),  rep("Pigment",8),
    rep("Blasto",4),  rep("Pigment",4),
    rep("Blasto",6),  rep("Pigment",6)),
  value=c(         
    c(45,	55,	57,	41,	43,	
      76,	41,	79,	54,	49),
    
    c(108,	124,	90,	143,	174,	104,	90,	82,	88,	60,
      116,	45,	133,	86,	97,	107	,180,	133,	79,	116),
    
    c(20,	25,	16,	13,	25,	20,	11,	37,	21,
      3,	0,	0,	0,	9,	0,	3,	0,	5),
    
    c(70,	46,	50,	30,	36,	35,	28,	31,
      6,	1,	8,37,	40,	33,	25,	13),
    
    c(25,	55,	8,	47, 
      8,	33,	4,	5),
    
    c(10,	17,	21,	12,	7,	27,
      1,	0	,4,	0	,2,	0)
  )
)



# Create a custom order for cell types
experiment_order <- c(    "Control",
                          "Mmless Delta+++",
                          "Mmless Control",
                          "Mmless Delta 6hr",
                          "Mmless Delta 8hr",
                          "Mmless Delta 10hr")



my_colors <- c("red",  "yellow")
cell_order <- c("Pigment", "Blasto")


data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Delta KDs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 200)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

kruskal.test( value ~ experiment , data = subset(data, celltype == "Blasto"))
kruskal.test( value ~ experiment , data = subset(data, celltype == "Pigment"))

subset_data1 <- subset(data, experiment %in% c('Control', "Mmless Delta+++"))
subset_data2 <- subset(data, experiment %in% c('Control', "Mmless Control"))
subset_data3 <- subset(data, experiment %in% c('Control', "Mmless Delta 6hr"))
subset_data4 <- subset(data, experiment %in% c('Control', "Mmless Delta 8hr"))
subset_data5 <- subset(data, experiment %in% c('Control', "Mmless Delta 10hr"))



results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd4 <- subset_data4 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd5 <- subset_data5 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1
results_sd2
results_sd3
results_sd4
results_sd5

### FIGURE 7L #######################################
# Control vs Micromereless  with Nodal MO + Delta MO

data <- data.frame(
  experiment=c( rep("Control",12), 
                rep("Control Nodal-Delta KD",9), 
                rep("Mmless Nodal-Delta KD",12)),
  celltype=c( rep("PMCs",4), rep("Blasto",4), rep("Pigment",4),
              rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",4), rep("Blasto",4), rep("Pigment",4)),
  value=c( c(56,	58,	47, 56,
             67,	47,	67, 38,
             32,	77,	48, 97),           
           c(93,	78,	94,
             37,	2, 38,
             0,	0,0 ),
           c(102,	104,	73,	106,	
             13,	2,	31,	18,	
             0,	0,	0,	0 ))
)


# Create a custom order for cell types
experiment_order <- c("Control", "Control Nodal-Delta KD", "Mmless Nodal-Delta KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Nodal-Delta KD",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 120)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results


subset_data1 <- subset(data, experiment %in% c('Control', 'Control Nodal-Delta KD'))
subset_data2 <- subset(data, experiment %in% c('Control Nodal-Delta KD', 'Mmless Nodal-Delta KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Nodal-Delta KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1
results_sd2
results_sd3



##############################################################
### SUPPLEMENTARY FIGURES 

# Figure S1: For triangle plots, see WOT pipelines
# Figure S1B  # Pigment cells vs PMCs' (Skeletogenic Cells [SC])

PMC_16hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_16_threshold_0.7_barcodes.txt", header = T)
PMC_14hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_14_threshold_0.7_barcodes.txt", header = T)
PMC_12hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_12_threshold_0.7_barcodes.txt", header = T)
PMC_10hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_10_threshold_0.7_barcodes.txt", header = T)
PMC_08hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_8_threshold_0.7_barcodes.txt", header = T)
PMC_06hpf_0.70 <- read.table("control/Waddington-OT/PMCs_hpf_6_threshold_0.7_barcodes.txt", header = T)

Pigment_16hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_16_threshold_0.7_barcodes.txt", header = T)
Pigment_14hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_14_threshold_0.7_barcodes.txt", header = T)
Pigment_12hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_12_threshold_0.7_barcodes.txt", header = T)
Pigment_10hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_10_threshold_0.7_barcodes.txt", header = T)
Pigment_08hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_8_threshold_0.7_barcodes.txt", header = T)
Pigment_06hpf_0.70 <- read.table("control/Waddington-OT/Pigment_hpf_6_threshold_0.7_barcodes.txt", header = T)



PMC_16hpf_0.70_seurat = subset(control_6to18,cells=PMC_16hpf_0.70$Cell_Barcode)
PMC_16hpf_0.70_umap <- as.data.frame(Embeddings(PMC_16hpf_0.70_seurat, "umap"))

PMC_14hpf_0.70_seurat = subset(control_6to18,cells=PMC_14hpf_0.70$Cell_Barcode)
PMC_14hpf_0.70_umap <- as.data.frame(Embeddings(PMC_14hpf_0.70_seurat, "umap"))

PMC_12hpf_0.70_seurat = subset(control_6to18,cells=PMC_12hpf_0.70$Cell_Barcode)
PMC_12hpf_0.70_umap <- as.data.frame(Embeddings(PMC_12hpf_0.70_seurat, "umap"))

PMC_10hpf_0.70_seurat = subset(control_6to18,cells=PMC_10hpf_0.70$Cell_Barcode)
PMC_10hpf_0.70_umap <- as.data.frame(Embeddings(PMC_10hpf_0.70_seurat, "umap"))

PMC_08hpf_0.70_seurat = subset(control_6to18,cells=PMC_08hpf_0.70$Cell_Barcode)
PMC_08hpf_0.70_umap <- as.data.frame(Embeddings(PMC_08hpf_0.70_seurat, "umap"))

PMC_06hpf_0.70_seurat = subset(control_6to18,cells=PMC_06hpf_0.70$Cell_Barcode)
PMC_06hpf_0.70_umap <- as.data.frame(Embeddings(PMC_06hpf_0.70_seurat, "umap"))




Pigment_16hpf_0.70_seurat = subset(control_6to18,cells=Pigment_16hpf_0.70$Cell_Barcode)
Pigment_16hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_16hpf_0.70_seurat, "umap"))

Pigment_14hpf_0.70_seurat = subset(control_6to18,cells=Pigment_14hpf_0.70$Cell_Barcode)
Pigment_14hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_14hpf_0.70_seurat, "umap"))

Pigment_12hpf_0.70_seurat = subset(control_6to18,cells=Pigment_12hpf_0.70$Cell_Barcode)
Pigment_12hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_12hpf_0.70_seurat, "umap"))

Pigment_10hpf_0.70_seurat = subset(control_6to18,cells=Pigment_10hpf_0.70$Cell_Barcode)
Pigment_10hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_10hpf_0.70_seurat, "umap"))

Pigment_08hpf_0.70_seurat = subset(control_6to18,cells=Pigment_08hpf_0.70$Cell_Barcode)
Pigment_08hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_08hpf_0.70_seurat, "umap"))

Pigment_06hpf_0.70_seurat = subset(control_6to18,cells=Pigment_06hpf_0.70$Cell_Barcode)
Pigment_06hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_06hpf_0.70_seurat, "umap"))


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "16_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_16hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_16hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "14_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_14hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_14hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "12_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_12hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_12hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "10_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_10hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_10hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "8_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_08hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_08hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=control2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(control2,  Stage == "6_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_06hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_06hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


#### Figure S1C 

PMC_16hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_16_threshold_0.7_barcodes.txt", header = T)
PMC_14hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_14_threshold_0.7_barcodes.txt", header = T)
PMC_12hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_12_threshold_0.7_barcodes.txt", header = T)
PMC_10hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_10_threshold_0.7_barcodes.txt", header = T)
PMC_08hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_8_threshold_0.7_barcodes.txt", header = T)
PMC_06hpf_0.70 <- read.table("mmless/Waddington-OT/PMCs_hpf_6_threshold_0.7_barcodes.txt", header = T)

Pigment_16hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_16_threshold_0.7_barcodes.txt", header = T)
Pigment_14hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_14_threshold_0.7_barcodes.txt", header = T)
Pigment_12hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_12_threshold_0.7_barcodes.txt", header = T)
Pigment_10hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_10_threshold_0.7_barcodes.txt", header = T)
Pigment_08hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_8_threshold_0.7_barcodes.txt", header = T)
Pigment_06hpf_0.70 <- read.table("mmless/Waddington-OT/Pigment_hpf_6_threshold_0.7_barcodes.txt", header = T)



PMC_16hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_16hpf_0.70$Cell_Barcode)
PMC_16hpf_0.70_umap <- as.data.frame(Embeddings(PMC_16hpf_0.70_seurat, "umap"))

PMC_14hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_14hpf_0.70$Cell_Barcode)
PMC_14hpf_0.70_umap <- as.data.frame(Embeddings(PMC_14hpf_0.70_seurat, "umap"))

PMC_12hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_12hpf_0.70$Cell_Barcode)
PMC_12hpf_0.70_umap <- as.data.frame(Embeddings(PMC_12hpf_0.70_seurat, "umap"))

PMC_10hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_10hpf_0.70$Cell_Barcode)
PMC_10hpf_0.70_umap <- as.data.frame(Embeddings(PMC_10hpf_0.70_seurat, "umap"))

#PMC_08hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_08hpf_0.70$Cell_Barcode)
#PMC_08hpf_0.70_umap <- as.data.frame(Embeddings(PMC_08hpf_0.70_seurat, "umap"))

#PMC_06hpf_0.70_seurat = subset(micromereless_6to18,cells=PMC_06hpf_0.70$Cell_Barcode)
#PMC_06hpf_0.70_umap <- as.data.frame(Embeddings(PMC_06hpf_0.70_seurat, "umap"))




Pigment_16hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_16hpf_0.70$Cell_Barcode)
Pigment_16hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_16hpf_0.70_seurat, "umap"))

#Pigment_14hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_14hpf_0.70$Cell_Barcode)
#Pigment_14hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_14hpf_0.70_seurat, "umap"))

#Pigment_12hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_12hpf_0.70$Cell_Barcode)
#Pigment_12hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_12hpf_0.70_seurat, "umap"))

#Pigment_10hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_10hpf_0.70$Cell_Barcode)
#Pigment_10hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_10hpf_0.70_seurat, "umap"))

#Pigment_08hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_08hpf_0.70$Cell_Barcode)
#Pigment_08hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_08hpf_0.70_seurat, "umap"))

#Pigment_06hpf_0.70_seurat = subset(micromereless_6to18,cells=Pigment_06hpf_0.70$Cell_Barcode)
#Pigment_06hpf_0.70_umap <- as.data.frame(Embeddings(Pigment_06hpf_0.70_seurat, "umap"))


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "16_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_16hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  geom_point(data=Pigment_16hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="plum",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "14_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_14hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  #geom_point(data=Pigment_14hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "12_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_12hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  #geom_point(data=Pigment_12hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "10_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  geom_point(data=PMC_10hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  #geom_point(data=Pigment_10hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "8_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="gray55",  pch=20, alpha=1) + 
  #geom_point(data=PMC_08hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  #geom_point(data=Pigment_08hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  
  theme_bw() + 
  ggtitle("") +  theme_void()


ggplot() + 
  geom_point(data=mmless2,  aes(x=(UMAP_1), y=(UMAP_2)), color="lightgrey", size=0.1, pch=20, alpha=0.5) +
  
  geom_point(data=subset(mmless2,  Stage == "6_HPF"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.1,  color="gray55",  pch=20, alpha=1) + 
  #geom_point(data=PMC_06hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="green",  pch=20, alpha=1) + 
  #geom_point(data=Pigment_06hpf_0.70_umap,  aes(x=(UMAP_1), y=(UMAP_2)), size=0.1,  color="red",  pch=20, alpha=1) + 
  
  theme_bw() + 
  ggtitle("") +  theme_void()


# Figure S2: Heatmap of similarity in integration

seu <- mmless_integrated_sct

p1 <- scatterPlot(seu, "pca",colour.by = "experiment", title = "PCA") 
p2 <- scatterPlot(seu, "umap",colour.by = "experiment", title = "UMAP") 
plot_grid(p1, p2)

seu$initial_cluster <- paste(seu$Stage, seu$experiment, sep = "_")

Idents(seu) <- "initial_cluster"



table(seu$initial_cluster)


ider <- getIDEr(seu, 
                group.by.var = "initial_cluster",
                batch.by.var = "experiment",
                downsampling.size = 35, 
                use.parallel = FALSE, verbose = FALSE)

celltimes <- c(
  "6_HPF"   ,
  "8_HPF"  ,
  "10_HPF"   ,
  "12_HPF"   ,
  "14_HPF"   ,
  "16_HPF"   ,
  "18_HPF"    )

idx1 <- paste0(celltimes, "_CONTROL")
idx2 <- paste0(celltimes, "_MMLESS")


pheatmap::pheatmap(
  ider[[1]][idx1, idx2],
  color = inferno(100),
  border_color = "gray",
  display_numbers = TRUE,
  cluster_rows = F,
  cluster_cols = F,
  width = 10,
  fontsize_number = 15,
  height = 10,
  cellwidth = 40,
  cellheight = 40
)

# Figure S3: Dotplot

##############################################################################################
## DOTPLOTS

Gene <- "SoxC"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
SoxC <- features[1]

Gene <- "Delta"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Delta <- features[1]

Gene <- "LVA-19158.t1:Th"

features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Th <- features[1]

Gene <- "Ngn"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Ngn <- features[1]


Gene <- "Trk"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Trk <- features[1]

Gene <- "Acsc"

features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Acsc <- features[1]

Gene <- "LVA-34741.t1:Brn1/2/4"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Brn <- features[1]

Gene <- "Insm"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Insm <- features[1]

## APD

Gene <- "FoxQ2"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
FoxQ2 <- features[1]

Gene <- "Hbn"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Hbn1 <- features[1]

Gene <- "Six3"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Six31 <- features[1]

Gene <- "Sfrp"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Sfrp <- features[1]

Gene <- "Zic"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Zic <- features[1]

Gene <- "Wnt5"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Wnt51 <- features[1]
Wnt52 <- features[2]

Gene <- "Hox7"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Hox7 <- features[1]

Gene <- "IrxA"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
IrxA <- features[1]

Gene <- "HesD"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
HesD <- features[1]

Gene <- "4258.t1:Nkx2-2"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Nkx2 <- features[1]

Gene <- "Msx"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Msx <- features[1]

Gene <- "Gsc"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Gsc <- features[1]

Gene <- "Chordin"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Chordin <- features[1]

Gene <- "Nodal"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Nodal <- features[1]

Gene <- "Lefty"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Lefty <- features[1]

Gene <- "LVA-22122.t1:Not"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Not <- features[1]


Gene <- "PaxB"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
PaxB <- features[1]


Gene <- "Univin"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Univin1 <- features[1]
Univin2 <- features[2]



Gene <- "Isl"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Isl <- features[1]

Gene <- "Btub1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Btub11 <- features[1]
Btub12 <- features[2]
Btub13 <- features[3]

Gene <- "Dach1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Dach1 <- features[1]


Gene <- "FoxG"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
FoxG <- features[1]

Gene <- "Dr-1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Dr <- features[1]

Gene <- "Hmx"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Hmx <- features[1]

Gene <- "Emx"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Emx <- features[1]

Gene <- "Onecut"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Onecut <- features[1]
Onecut2 <- features[2]


Gene <- "Z166"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Z166 <- features[1]

# Early Endoderm
Gene <- "Fzd1/2/7"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Fzd127 <- features[1]

Gene <- "Fzd5/8"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Fzd58 <- features[1]

Gene <- "SoxB"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
SoxB <- features[1]

Gene <- "He1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
He1 <- features[1]

Gene <- "Ost"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Ost <- features[1]

Gene <- "Smadip"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Smadip <- features[1]


# Gut
Gene <- "Endo16"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Endo16 <- features[1]


Gene <- "FoxA"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
FoxA <- features[2]

Gene <- "Bra"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Bra <- features[1]

Gene <- "Blimp"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Blimp <- features[1]
Blimp2 <- features[2]

Gene <- "Gatae"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Gatae <- features[1]
Gatae2 <- features[2]
Gatae3 <- features[3]

#Hh, , Npnt, Wnt8, Wnt1, Gat2,            # Endoderm
#LVA-1800.t1:Notch
#LVA-30544.t1:Hox11/13b
#LVA-26446.t1:Moap1L

#LVA-9056.t1:Nlr131

#LVA-32173.t1:Wnt5



Gene <- "LVA-30168.t1:Hh"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Hh <- features[1]

Gene <- "LVA-11480.t1:Npnt"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Npnt <- features[1]

Gene <- "Wnt8"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Wnt8 <- features[1]

Gene <- "Wnt1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Wnt16 <- features[1]
Wnt1 <- features[2]

Gene <- "LVA-3505.t1:FoxN2/3"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
FoxN <- features[1]


Gene <- "Eve"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Eve <- features[1]

Gene <- "Hox11"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Hox11 <- features[1]


# MESODERM
#, , Nlk, Erg, Sarmr, SCl, Elav 
Gene <- "GataC"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
GataC <- features[1]

Gene <- "Bpnt"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Bpnt <- features[1]

Gene <- "Nlk"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Nlk <- features[1]

Gene <- "Erg"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Erg <- features[1]

Gene <- "Scl"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Scl <- features[2]

Gene <- "Elav"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Elav <- features[1]


#, Tgif, FoxN,                             
# Endomesoderm
Gene <- "Ese"

features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Ese <- features[1]
Ese2 <- features[2]

Gene <- "Tgif"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Tgif <- features[1]


#           
# PGCs
Gene <- "LVA-26071.t1:Nanos2"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Nanos <- features[1]


# Micromeres
Gene <- "Jun"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Jun <- features[1]

Gene <- "LVA-5802.t1:Alx1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Alx1 <- features[1]

Gene <- "LVA-36687.t1:ElavL"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
ElavL <- features[1]

# PMCs
Gene <- "Tbr"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Tbr <- features[1]

Gene <- "LVA-10710.t1:Pks2"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Pks2 <- features[1]

Gene <- "Sm29"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Sm29 <- features[1]

Gene <- "LVA-28884.t1:Msp130-1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Msp130 <- features[1]

Gene <- "LVA-13092.t1:p58-a"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
P58 <- features[1]



# Pigment
Gene <- "Acaca"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Acaca <- features[1]

Gene <- "Pks1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Pks1 <- features[1]

Gene <- "Mif5"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Mif5 <- features[1]

Gene <- "LVA-4074.t1:Six1/2"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Six12 <- features[1]

Gene <- "LVA-26843.t1:Gcm"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Gcm <- features[1]



# Blastocoelar
Gene <- "Astacin4"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Astacin4 <- features[1]


Gene <- "LVA-13883.t1:Erg"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Erg <- features[1]

Gene <- "LVA-3189.t1:Prox1"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Prox <- features[1]

Gene <- "LVA-32108.t1:GataC"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
GataC <- features[1]

Gene <- "Ron"
features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
Ron <- features[1]


Celltype_features <- c(FoxQ2, Sfrp, Zic, Hbn1, Six31,              # APD
                       Delta, Acsc, Brn,SoxC,Th, Trk,   Insm,      # Neural
                       Emx, Hmx,FoxG,Dr,  Onecut,          # Ciliary Band
                       
                       IrxA,  Nkx2,HesD, Wnt51, Hox7,Msx, Eve,            # Aboral
                       Chordin,   Lefty,Nodal, Gsc,Not,            # Oral
                       Univin2, Isl, #Btub13,                       # Ectoderm
                       He1,Fzd127, Fzd58, SoxB,  Smadip,           # Early Ectoderm
                       
                       Endo16, FoxA, Bra, Blimp, Gatae3,            # Gut
                       Hh, Npnt, Wnt8, Wnt1, Wnt16, FoxN,            # Endoderm
                       # Early Endoderm
                       
                       Ese2, Tgif,                              # Endomesoderm
                       Bpnt, Nlk, Scl,            # Mesoderm
                       Nanos,                  # PGCs
                       Jun, Alx1, ElavL,                              # Micromeres
                       Tbr, Pks2, Sm29, Msp130, P58   ,       # PMCs
                       Acaca,Pks1,Mif5,Six12,Gcm,                     # Pigment
                       Astacin4, Erg, Prox, GataC,  Ron                     # Blastocoelar
)

unique(Celltype_features)

list1 <- unique(control_6to18@meta.data$celltype)

# Elements in list1 but not in list2
unique_to_list1 <- setdiff(list1, list2)

# Print the results
print(unique_to_list1)

head(control_6to18@meta.data)


Idents(control_6to18) <- "celltype"
# Set the order of the identities (clusters or cell types)
control_6to18@active.ident <- factor(control_6to18@active.ident, levels = c("APD" , "Neural"  ,"Ciliary_Band" ,"Aboral_Ectoderm" ,  "Oral_Ectoderm", "Ectoderm", "Early_Ectoderm" ,
                                                                            "Gut", "Endoderm"   ,"Early_Endoderm","Endomesoderm" , "Mesoderm"   ,"PGCs"   ,
                                                                            "Micromeres"     ,  "PMCs" , "Pigment"    , "Blastocoelar"  , "Undefined"  ))

DotPlot(control_6to18, features = Celltype_features, cluster.idents = F,   col.min = -3,
        col.max = 3) + RotatedAxis() + #coord_flip() +
  #scale_color_gradient2(high='red', mid='white', low='darkblue') + 
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))


# Figure S4: Delta expression in Controls vs Mmless UMAP
library(scCustomize)

n=10
pal <-inferno(n, alpha = 1, begin = 0, end = 1, direction = 1)

Gene <-"LVA-3923.t1:Delta"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

features <- rownames(micromereless_6to18[grep(Gene,rownames(micromereless_6to18))])
features
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()


# Figure S5A: msp130, pks1, astacin4 and irf4 expression in Controls vs Mmless UMAP

Gene <-"LVA-28884.t1:Msp130-1"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void() 
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void() 

Gene <-"LVA-2666.t1:Pks1"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void() 
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void() 

Gene <-"LVA-23276.t1:Astacin4"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

Gene <-"LVA-23179.t1:Irf4"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void() 


# Figure S5B: bra, blimp1, alx1 expression in Controls vs Mmless UMAP

Gene <-"LVA-34084.t1:Bra"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

Gene <-"LVA-33680.t1:Blimp1"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

Gene <-"LVA-5802.t1:Alx1"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

# Figure S5C: wnt1, wnt8 expression in Controls vs Mmless UMAP

Gene <-"LVA-12953.t1:Wnt1"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()

Gene <-"LVA-6974.t1:Wnt8"
features <- rownames(control_6to18[grep(Gene,rownames(control_6to18))])
features
FeaturePlot_scCustom(control_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()
FeaturePlot_scCustom(micromereless_6to18, features = features[1] ,colors_use = viridis_inferno_dark_high) +  
  theme_void()


# Figure S6: Blimp1 MO
data <- data.frame(
  experiment=c( rep("Control",12), 
                rep("Control Blimp KD",9), 
                rep("Mmless Blimp KD",15)),
  celltype=c( rep("PMCs",4), rep("Blasto",4), rep("Pigment",4),
              rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",5), rep("Blasto",5), rep("Pigment",5)),
  value=c( c(66,	60,	65,	67,	37,	34,	43,	34,	31,	36,	21,	45),           
           c(49,	40,	36,	45,	44,	32,	28,	35,	18),
           c(35,	73,	28,	55,	41,	16,	21,	5,	6,	33,	0,	0,	0,	0,	2))
)

experiment_order <- c("Control", "Control Blimp KD", "Mmless Blimp KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Blimp KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 100)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

subset_data1 <- subset(data, experiment %in% c('Control', 'Control Blimp KD'))
subset_data2 <- subset(data, experiment %in% c('Control Blimp KD', 'Mmless Blimp KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Blimp KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd3




results_sd1
results_sd2
results_sd3

###########################################33
# Figure S7: wnt1 MO

data <- data.frame(
  experiment=c( 
    rep("Control",6), 
    rep("Control Wnt1 KD",12),
    rep("Mmless Wnt1 KD",27)),
  celltype=c( 
    rep("PMCs",2), rep("Blasto",2), rep("Pigment",2),
    rep("PMCs",4), rep("Blasto",4), rep("Pigment",4),
    rep("PMCs",9), rep("Blasto",9), rep("Pigment",9)),
  value=c(          
    c(57,	64,	85,	82,	40,	50),
    c(51,	53,	44,	49,	43,	48,	47,	56,	55,	48,	60,	51),
    c(54,	109,	116,	80,	81,	124,	116,	62,	108,	24,	21,	30,	28,	20,	18,	25,	48,	35,	53,	26,	23,	23,	51,	17,	42,	21,	44))
)

# Create a custom order for cell types
experiment_order <- c("Control", "Control Wnt1 KD", "Mmless Wnt1 KD")

my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")
data <- transform(data, celltype = factor(celltype, levels = cell_order))

# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Wnt1 KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 125)


results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )


subset_data1 <- subset(data, experiment %in% c('Control', 'Control Wnt1 KD'))
subset_data2 <- subset(data, experiment %in% c('Control Wnt1 KD', 'Mmless Wnt1 KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Wnt1 KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1

results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd3
# Convert the results to a data frame
results_df <- as.data.frame(results)

# Print the results
print(results_df)

results_sd1
results_sd2
results_sd3


# Figure S8: wnt8 MO

data <- data.frame(
  experiment=c( rep("Control",12), 
                rep("Control Wnt8 KD",9), 
                rep("Mmless Wnt8 KD",15)),
  celltype=c( rep("PMCs",4), rep("Blasto",4), rep("Pigment",4),
              rep("PMCs",3), rep("Blasto",3), rep("Pigment",3),
              rep("PMCs",5), rep("Blasto",5), rep("Pigment",5)),
  value=c( c(66,	60,	65,	67,	37,	34,	43,	33,	31,	36,	21,	45),           
           c(62,	55,	60,	33,	38,	41,	31,	11,	33),
           c(15,	49,	36,	66,	33,	0,	40,	0,	5,	33,	0,	33,	0,	0,	2))
)


# Create a custom order for cell types
experiment_order <- c("Control", "Control Wnt8 KD", "Mmless Wnt8 KD")
my_colors <- c("green",  "red", "yellow")
cell_order <- c("PMCs", "Pigment", "Blasto")

data <- transform(data, celltype = factor(celltype, levels = cell_order))
data


# Create the boxplot with facets (without facet titles and legend title)
ggplot(data, aes(x = factor(experiment, levels = experiment_order), 
                 y = value, 
                 fill = factor(celltype, levels = cell_order))) +
  geom_boxplot() +
  labs(title = "Effect of Wnt8 KOs",
       x = "Cell Type",
       y = "Cell Count",
       fill = "") +
  theme_bw() +
  facet_wrap(~ celltype, scales = "free", strip.position = "top") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill = "white", colour = NA)) +
  scale_fill_manual(values = my_colors) + ylim(0, 100)




results <- data %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results

subset_data1 <- subset(data, experiment %in% c('Control', 'Control Wnt8 KD'))
subset_data2 <- subset(data, experiment %in% c('Control Wnt8 KD', 'Mmless Wnt8 KD'))
subset_data3 <- subset(data, experiment %in% c('Control', 'Mmless Wnt8 KD'))

results_sd1 <- subset_data1 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )
results_sd1


results_sd2 <- subset_data2 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd2

results_sd3 <- subset_data3 %>%
  group_by(celltype) %>%
  summarise(
    p_value = kruskal.test(value ~ experiment )$p.value
  )

results_sd1
results_sd2
results_sd3

# Figure S10: Umap Clusters

ggplot() + 
  geom_point(data=subset(control2,  seurat_clusters == "0"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#BEBADA",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "1"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#8DD3C7",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "2"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B3E2CD",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "3"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#66A61E",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "4"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E5D8BD",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "5"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#666666",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "6"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDCDAC",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "7"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FFFF99",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "8"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B15928",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "9"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F781BF",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "10"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB4443",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "11"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB8072",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "12"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#8DA0CB",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "13"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FFD92F",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "14"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#7FC97F",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "15"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A6CEE3",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "16"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A65628",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "17"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCEBC5",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "18"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F4CAE4",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "19"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCEBC5",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "20"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#332323",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "21"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CBD5E8",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "22"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A44F22",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "23"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F3422C",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "24"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E7298A",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "25"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB9A99",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "26"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#6A3D9A",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "27"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#DECBE4",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "28"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E41A1C",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "29"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCCCCC",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "30"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A6761D",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "31"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FC8D62",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "32"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E78AC3",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "33"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#099018",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "34"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDB462",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "35"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#4DAF4A",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "36"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#80B1D3",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "37"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F0027F",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "38"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDC086",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "39"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FCCDE5",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "40"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B3DE69",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "41"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#D9D9D9",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "42"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F2F2F2",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "43"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDC086",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "44"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FCCDE5",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "45"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B3DE69",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "46"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#D9D9D9",  pch=20, alpha=1) + 
  geom_point(data=subset(control2,  seurat_clusters == "47"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F2F2F2",  pch=20, alpha=1) + 
  ggtitle("") +  theme_void()



ggplot() + 
  geom_point(data=subset(mmless2,  seurat_clusters == "0"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#BEBADA",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "1"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#8DD3C7",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "2"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B3E2CD",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "3"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#66A61E",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "4"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E5D8BD",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "5"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#666666",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "6"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDCDAC",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "7"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FFFF99",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "8"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B15928",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "9"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F781BF",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "10"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB4443",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "11"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB8072",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "12"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#8DA0CB",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "13"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FFD92F",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "14"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#7FC97F",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "15"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A6CEE3",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "16"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A65628",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "17"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCEBC5",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "18"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F4CAE4",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "19"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCEBC5",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "20"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#332323",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "21"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CBD5E8",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "22"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A44F22",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "23"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F3422C",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "24"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E7298A",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "25"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FB9A99",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "26"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#6A3D9A",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "27"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#DECBE4",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "28"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E41A1C",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "29"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#CCCCCC",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "30"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#A6761D",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "31"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FC8D62",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "32"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#E78AC3",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "33"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#099018",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "34"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDB462",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "35"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#4DAF4A",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "36"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#80B1D3",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "37"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F0027F",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "38"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FDC086",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "39"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#FCCDE5",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "40"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#B3DE69",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "41"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#D9D9D9",  pch=20, alpha=1) + 
  geom_point(data=subset(mmless2,  seurat_clusters == "42"),  aes(x=(UMAP_1), y=(UMAP_2)),size=0.5,  color="#F2F2F2",  pch=20, alpha=1) + 
  
  theme_bw() + 
  
  ggtitle("") +  theme_void()















