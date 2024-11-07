library(Seurat)
library(tidyverse)

# SCTransform the data

load("precomputed_micromereless.Rda", verbose = T)
load("precomputed_control.Rda", verbose = T)


############################################################################
############################################################################
######   INTEGRATION 
# Running integration with dataset


object_list <- list(control_6to18, micromereless_6to18)
features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3500)
object_list_2 <- PrepSCTIntegration(object.list = object_list, anchor.features = features)
mmless.anchors <- FindIntegrationAnchors(object.list = object_list_2, normalization.method="SCT", anchor.features= features, reference = 1, dims = 1:50)
mmless_integrated_sct <- IntegrateData(anchorset = mmless.anchors, dims = 1:50, normalization.method="SCT" )
mmless_integrated_sct <- ScaleData(mmless_integrated_sct, verbose = FALSE)
mmless_integrated_sct <- RunPCA(mmless_integrated_sct, verbose = T)
mmless_integrated_sct <- RunUMAP(mmless_integrated_sct, reduction= "pca" , dims = 1:45)
mmless_integrated_sct <- FindNeighbors(mmless_integrated_sct,  dims = 1:45)


features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 4000)
object_list_2 <- PrepSCTIntegration(object.list = object_list, anchor.features = features)
mmless.anchors <- FindIntegrationAnchors(object.list = object_list_2, normalization.method="SCT", anchor.features= features, reference = 1, dims = 1:50)
mmless_integrated_sctb <- IntegrateData(anchorset = mmless.anchors, dims = 1:50, normalization.method="SCT" )
mmless_integrated_sctb <- ScaleData(mmless_integrated_sctb, verbose = FALSE)
mmless_integrated_sctb <- RunPCA(mmless_integrated_sctb, verbose = T)
mmless_integrated_sctb <- RunUMAP(mmless_integrated_sctb, reduction= "pca" , dims = 1:20)
mmless_integrated_sctb <- FindNeighbors(mmless_integrated_sctb,  dims = 1:20)

save(mmless_integrated_sct , mmless_integrated_sctb,  file="precomputed_integrated.Rda")

print('all done!')


