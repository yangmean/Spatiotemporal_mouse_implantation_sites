library("SCENIC")
library("Seurat")
load("20201203/UE-merge-20201221-UMAP-rename.RData")
exprMat <- as.matrix(GetAssayData(UE.seurat@assays$RNA, slot = "counts"))
cellInfo <- UE.seurat@meta.data
cellInfo <- UE.seurat@meta.data
dbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="mgi", dbDir="../../../00.learning/scenic/databases/", dbs=dbs, nCores = 40)
### 共表达网络
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
#run GENIE3
runGenie3(exprMat = exprMat, scenicOptions = scenicOptions)
