#monocle
library("monocle")
load("20201203/UE-decidual-20201202-UMAP.RData")
UE.merge.decidual <- subset(UE.merge.decidual, subset = seurat_clusters != 9 & seurat_clusters != 47 & seurat_clusters != 38)
UE.merge.decidual.monocle <- subset(UE.merge.decidual, subset = orig.ident %in% c("U55.E65", "U65", "U85.1", "U85.2", "U95", "U105", "U.E75"))
UE.data <- as(as.matrix(UE.merge.decidual.monocle@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.merge.decidual.monocle@meta.data)
fd <- data.frame(gene_short_name = row.names(UE.data), row.names = row.names(UE.data))
fd <- new("AnnotatedDataFrame", data = fd)
UE.cds <- newCellDataSet(UE.data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
#过滤
UE.cds <- estimateSizeFactors(UE.cds)
UE.cds <- estimateDispersions(UE.cds)
#统计在设定基因的最低表达量为0.1的情况下，每个细胞表达多少个基因，以及每个基因有多少细胞表达
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
#保留下来至少在100个细胞中有表达量的基因，然后进行后续分析
expressed_genes <- rownames(subset(fData(UE.cds), num_cells_expressed >= 100))
length(expressed_genes)#15920个基因
#不再过滤，直接进行拟时间分析
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
#根据样本时间挑选差异基因分析
clustering_DEG_genes <- differentialGeneTest(UE.cds[expressed_genes, ],
                                             fullModelFormulaStr = "~orig.ident",
                                             cores = 10)
clustering_DEG_genes <- write.table(clustering_DEG_genes, file = "UE-decidual-DEGs-sample.csv", sep = ",", quote = F, row.names = F)
UE.cds_ordering_genes <- clustering_DEG_genes$gene_short_name[order(clustering_DEG_genes$qval)][1:1000]
UE.cds <- setOrderingFilter(UE.cds, ordering_genes = UE.cds_ordering_genes)
UE.cds <- reduceDimension(UE.cds, method = "ICA")
UE.cds <- orderCells(UE.cds, num_paths = 1)
save(UE.cds, file = "20201203/UE-decidual-monocle-sample-1000-ICA-path1-20201224.RData")
load("20201203/UE-decidual-monocle-sample-1000-ICA-path1-20201224.RData")
plot_cell_trajectory(UE.cds, color_by = "orig.ident", cell_size = 0.5, show_branch_points = F)
plot_cell_trajectory(UE.cds, color_by = "cellname", cell_size = 0.5, show_branch_points = F)
plot_cell_trajectory(UE.cds, color_by = "Pseudotime", cell_size = 0.5, show_branch_points = F)
UE.merge.decidual$cellname <- UE.merge.decidual@active.ident
UE.merge.decidual <- subset(UE.merge.decidual, subset = seurat_clusters != 9 & seurat_clusters != 47 & seurat_clusters != 38)
UE.merge.decidual.monocle <- subset(UE.merge.decidual, subset = orig.ident %in% c("U55.E65", "U65", "U85.1", "U85.2", "U95", "U105", "U.E75"))
Pseudotime <- UE.cds$Pseudotime
names(Pseudotime) <- colnames(UE.merge.decidual.monocle)
UE.merge.decidual.monocle <- AddMetaData(UE.merge.decidual.monocle, metadata = Pseudotime, col.name = "pseudotime")
save(UE.merge.decidual.monocle, file = "20201203/UE-merge-20201224-UMAP-pseudotime.RData")
FeaturePlot(UE.merge.decidual.monocle, features = "pseudotime", cols = c("LightGoldenrod", "red"))



#variable gene 做monocle
UE.merge.decidual$cellname <- UE.merge.decidual@active.ident
UE.merge.decidual <- subset(UE.merge.decidual, subset = seurat_clusters != 9 & seurat_clusters != 47 & seurat_clusters != 38)
UE.merge.decidual.monocle <- subset(UE.merge.decidual, subset = orig.ident %in% c("U55.E65", "U65", "U85.1", "U85.2", "U95", "U105", "U.E75"))
UE.data <- as(as.matrix(UE.merge.decidual.monocle@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.merge.decidual.monocle@meta.data)
fd <- data.frame(gene_short_name = row.names(UE.data), row.names = row.names(UE.data))
fd <- new("AnnotatedDataFrame", data = fd)
UE.cds <- newCellDataSet(UE.data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
#过滤
UE.cds <- estimateSizeFactors(UE.cds)
UE.cds <- estimateDispersions(UE.cds)
#统计在设定基因的最低表达量为0.1的情况下，每个细胞表达多少个基因，以及每个基因有多少细胞表达
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
#保留下来至少在100个细胞中有表达量的基因，然后进行后续分析
expressed_genes <- rownames(subset(fData(UE.cds), num_cells_expressed >= 100))
length(expressed_genes)#15920个基因
#不再过滤，直接进行拟时间分析
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
UE.cds_ordering_genes <- VariableFeatures(UE.merge.decidual)
UE.cds <- setOrderingFilter(UE.cds, ordering_genes = UE.cds_ordering_genes)
UE.cds <- reduceDimension(UE.cds, method = "DDRTree")
UE.cds <- orderCells(UE.cds)
GM_state <- function(UE.cds){
  if(length(unique(pData(UE.cds)$State))>1){
    UE.start <- table(pData(UE.cds)$State, pData(UE.cds)$orig.ident)[, 3]
    return(as.numeric(names(UE.start)[which(UE.start==max(UE.start))]))
  }else{return(1)}
}
UE.cds <- orderCells(UE.cds, root_state = GM_state(UE.cds))
plot_cell_trajectory(UE.cds, color_by = "cellname", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = brewer.pal(12, "Paired"))
save(UE.cds, file = "20201203/UE-decidual-monocle-varaible-2000.RData")

plot_cell_trajectory(UE.cds, color_by = "orig.ident", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = brewer.pal(12, "Paired"))

#去掉U10.5时间点decidual
UE.merge.decidual$cellname <- UE.merge.decidual@active.ident
UE.merge.decidual <- subset(UE.merge.decidual, subset = seurat_clusters != 9 & seurat_clusters != 47 & seurat_clusters != 38)
UE.merge.decidual.monocle <- subset(UE.merge.decidual, subset = orig.ident %in% c("U55.E65", "U65", "U85.1", "U85.2", "U95", "U.E75"))
UE.data <- as(as.matrix(UE.merge.decidual.monocle@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.merge.decidual.monocle@meta.data)
fd <- data.frame(gene_short_name = row.names(UE.data), row.names = row.names(UE.data))
fd <- new("AnnotatedDataFrame", data = fd)
UE.cds <- newCellDataSet(UE.data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
#过滤
UE.cds <- estimateSizeFactors(UE.cds)
UE.cds <- estimateDispersions(UE.cds)
#统计在设定基因的最低表达量为0.1的情况下，每个细胞表达多少个基因，以及每个基因有多少细胞表达
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
#保留下来至少在100个细胞中有表达量的基因，然后进行后续分析
expressed_genes <- rownames(subset(fData(UE.cds), num_cells_expressed >= 100))
length(expressed_genes)#15920个基因
#不再过滤，直接进行拟时间分析
UE.cds <- detectGenes(UE.cds, min_expr = 0.1)
UE.cds_ordering_genes <- VariableFeatures(UE.merge.decidual)
UE.cds <- setOrderingFilter(UE.cds, ordering_genes = UE.cds_ordering_genes)
UE.cds <- reduceDimension(UE.cds, method = "DDRTree")
UE.cds <- orderCells(UE.cds)
GM_state <- function(UE.cds){
  if(length(unique(pData(UE.cds)$State))>1){
    UE.start <- table(pData(UE.cds)$State, pData(UE.cds)$orig.ident)[, 2]
    return(as.numeric(names(UE.start)[which(UE.start==max(UE.start))]))
  }else{return(1)}
}
UE.cds <- orderCells(UE.cds, root_state = GM_state(UE.cds))
save(UE.cds, file = "20201203/UE-decidual-monocle-varaible-2000-remove-U105.RData")
plot_cell_trajectory(UE.cds, color_by = "cellname", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = brewer.pal(12, "Paired"))

plot_cell_trajectory(UE.cds, color_by = "orig.ident", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = brewer.pal(12, "Paired"))

plot_cell_trajectory(UE.cds, color_by = "Pseudotime", show_state_number = F, cell_size = 1, show_branch_points = F)

apoptosis <- read.table("20201203/03.decidual/GO-term-score/Apoptosis-GO-0006915.csv", sep = ",", header = T)
pos.apoptosis <- read.table("20201203/03.decidual/GO-term-score/positive-regulation-apoptotic-pathway-GO-2001235.csv", sep = ",", header = T)
decidual <- read.table("20201203/03.decidual/GO-term-score/Decidualization-mouse-GO-0046697.csv", sep = ",", header = T)
mitochondion <- read.table("20201203/03.decidual/GO-term-score/mitochondrion-GO-0005739.csv", sep = ",", header = T)
UE.merge.decidual.monocle.data <- data.frame(as.matrix(UE.merge.decidual.monocle@assays$RNA@data), check.names = F)
UE.merge.decidual.monocle.data.apoptosis <- UE.merge.decidual.monocle.data[rownames(UE.merge.decidual.monocle.data) %in% as.character(apoptosis$Symbol), ]
UE.merge.decidual.monocle.data.apoptosis.agv <- apply(UE.merge.decidual.monocle.data.apoptosis, 2, mean)

UE.merge.decidual.monocle.data.pos.apoptosis <- UE.merge.decidual.monocle.data[rownames(UE.merge.decidual.monocle.data) %in% as.character(pos.apoptosis$Symbol), ]
UE.merge.decidual.monocle.data.pos.apoptosis.agv <- apply(UE.merge.decidual.monocle.data.pos.apoptosis, 2, mean)

UE.merge.decidual.monocle.data.decidual <- UE.merge.decidual.monocle.data[rownames(UE.merge.decidual.monocle.data) %in% as.character(decidual$Symbol), ]
UE.merge.decidual.monocle.data.decidual.agv <- apply(UE.merge.decidual.monocle.data.decidual, 2, mean)

UE.merge.decidual.monocle.data.mitochondion <- UE.merge.decidual.monocle.data[rownames(UE.merge.decidual.monocle.data) %in% as.character(mitochondion$Symbol), ]
UE.merge.decidual.monocle.data.mitochondion.agv <- apply(UE.merge.decidual.monocle.data.mitochondion, 2, mean)

UE.cds$apoptosis <- UE.merge.decidual.monocle.data.apoptosis.agv
UE.cds$pos.apoptosis <- UE.merge.decidual.monocle.data.pos.apoptosis.agv
UE.cds$decidual <- UE.merge.decidual.monocle.data.decidual.agv
UE.cds$mitochondion <- UE.merge.decidual.monocle.data.mitochondion.agv

plot_cell_trajectory(UE.cds, color_by = "apoptosis", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 0.27)
plot_cell_trajectory(UE.cds, color_by = "pos.apoptosis", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 0.3)
plot_cell_trajectory(UE.cds, color_by = "decidual", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 0.5)
plot_cell_trajectory(UE.cds, color_by = "mitochondion", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 0.3)



