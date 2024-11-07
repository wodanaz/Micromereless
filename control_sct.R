library(Seurat)
library(tidyverse)


#Lv6hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-6hpf_anno.csv", row.names=1)
#Lv8hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-8hpf_anno.csv", row.names=1)
#Lv10hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-10hpf_anno.csv", row.names=1)
#Lv12hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-12hpf_anno.csv", row.names=1)
#Lv14hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-14hpf_anno.csv", row.names=1)
#Lv16hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-16hpf_anno.csv", row.names=1)
#Lv18hpf <- read.csv("/data/wraycompute/alejo/singlecell/atlas_combo_6/count_tables/Lv-18hpf_anno.csv", row.names=1)

# COnvert to Seurat object

#Lv6hpf_seurat <- CreateSeuratObject(counts = Lv6hpf, project = "Lv6hpf", min.cells = 3, min.features = 200)
#Lv8hpf_seurat <- CreateSeuratObject(counts = Lv8hpf, project = "Lv8hpf", min.cells = 3, min.features = 200)
#Lv10hpf_seurat <- CreateSeuratObject(counts = Lv10hpf, project = "Lv10hpf", min.cells = 3, min.features = 200)
#Lv12hpf_seurat <- CreateSeuratObject(counts = Lv12hpf, project = "Lv12hpf", min.cells = 3, min.features = 200)
#Lv14hpf_seurat <- CreateSeuratObject(counts = Lv14hpf, project = "Lv14hpf", min.cells = 3, min.features = 200)
#Lv16hpf_seurat <- CreateSeuratObject(counts = Lv16hpf, project = "Lv16hpf", min.cells = 3, min.features = 200)
#Lv18hpf_seurat <- CreateSeuratObject(counts = Lv18hpf, project = "Lv18hpf", min.cells = 3, min.features = 200)

# ADDING STAGE ANNOTATIONS FOR EACH TIME SERIES LIBRARY

#Lv6hpf_seurat$Stage <- "6_HPF"
#Lv8hpf_seurat$Stage <- "8_HPF"
#Lv10hpf_seurat$Stage <- "10_HPF"
#Lv12hpf_seurat$Stage <- "12_HPF"
#Lv14hpf_seurat$Stage <- "14_HPF"
#Lv16hpf_seurat$Stage <- "16_HPF"
#Lv18hpf_seurat$Stage <- "18_HPF"


#control_6to18_raw <-merge(Lv6hpf_seurat, y= c(Lv8hpf_seurat,
#                                                    Lv10hpf_seurat,
#                                                    Lv12hpf_seurat,
#                                                    Lv14hpf_seurat,
#                                                    Lv16hpf_seurat,
#                                                    Lv18hpf_seurat),
#                                
#                                add.cell.ids = c("6_HPF", 
#                                                 "8_HPF",
#                                                 "10_HPF" , 
#                                                 "12_HPF", 
#                                                 "14_HPF",
#                                                 "16_HPF",
#                                                 "18_HPF"),
#                                project = "Micromereless")


######### Put custom order to the groups

#control_6to18_raw$Stage <- factor(control_6to18_raw$Stage,  levels = c('6_HPF','8_HPF', '10_HPF', '12_HPF', '14_HPF', 
#                                                                                   '16_HPF','18_HPF'))

#control_6to18_raw$experiment <- "CONTROL"

#save(control_6to18_raw ,file="raw_dataset_single_cell_times_atlas.Rda")

# SCTransform the data
load("raw_dataset_single_cell_times_atlas.Rda", verbose=T)

control_6to18_raw$orig.ident <- factor(control_6to18_raw$orig.ident,  levels = c('Lv6hpf','Lv8hpf','Lv10hpf', 'Lv12hpf','Lv14hpf','Lv16hpf','Lv18hpf'))

control_6to18_raw <- PercentageFeatureSet(control_6to18_raw, pattern = "\\b\\w*Rp[sl]\\w*\\b", col.name = "percent.Rb")

pdf(file = "violin_plot_pre.pdf")
VlnPlot(control_6to18_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3)
dev.off()

control_6to18_raw <- subset(control_6to18_raw, subset = nFeature_RNA > 200 & nCount_RNA < 10000 )

pdf(file = "violin_plot_post.pdf")
VlnPlot(control_6to18_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.Rb"), ncol = 3)
dev.off()


#####  variable.features.n = 2500 and 2750 for atlas


control_6to18 <- SCTransform(control_6to18_raw, verbose = T, variable.features.n = 2500,  vars.to.regress = "percent.Rb")
control_6to18 <- RunPCA(control_6to18,  npcs = 200, verbose = T)
control_6to18 <- RunUMAP(control_6to18,  dims = 1:100)
control_6to18 <- FindNeighbors(control_6to18, dims = 1:100)
control_6to18 <- FindClusters(control_6to18, resolution = 1)

control_6to18b <- SCTransform(control_6to18_raw, verbose = T, variable.features.n = 2750,  vars.to.regress = "percent.Rb")
control_6to18b <- RunPCA(control_6to18b,  npcs = 200, verbose = T)
control_6to18b <- RunUMAP(control_6to18b,  dims = 1:100)
control_6to18b <- FindNeighbors(control_6to18b, dims = 1:100)
control_6to18b <- FindClusters(control_6to18b, resolution = 1)


save(control_6to18, control_6to18b,  file="precomputed_control.Rda")
