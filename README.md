# Seurat4_Analysis


***#Combine data and Run SCT tranfrom in Seurat 3***

Current SCT v0.3.1 along with Seurat 4 is having issues in normal progression (<a href="https://github.com/satijalab/seurat/issues/3618">Seurat issue</a>). 
```
pag.combined<- merge(Pool1.1, y=c(Pool2.1, Pool3.1, Pool4.1, Pool5.1,Pool6.1,Pool7.1,Pool8.1,Pool9.1,Pool10.1,Pool11.1,Pool12.1,Pool13.1,Pool14.1,Pool15.1,Pool16.1, Pool1.2, Pool2.2, Pool3.2, Pool4.2, Pool5.2,Poo
l6.2,Pool7.2,Pool8.2,Pool9.2,Pool10.2,Pool11.2,Pool12.2,Pool13.2,Pool14.2,Pool15.2,Pool16.2),project = "AllSamplesCombined")
saveRDS (pag.combined, file="Batch1_2Samples.Rds")

pag.combined[["percent.mt"]] <- PercentageFeatureSet(pag.combined, pattern = "^MT-")
pag.combined <- subset(pag.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5200 & percent.mt < 10)
pag.combined <- SCTransform(pag.combined, vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS (pag.combined, file="Batch1_2Samples_SCTransform.Rds")
```
***Mapping to PBMC multimodal reference datasets from Seurat***

Reference Location: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat

```
reference <- LoadH5Seurat("Seurat4/PBMC_REF/pbmc_multimodal.h5seurat")

MyData <- readRDS("Batch1_2Samples_SCTransform.Rds")


anchors <- FindTransferAnchors(
  reference = reference,
  query = MyData,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)


MyData <- MapQuery(
  anchorset = anchors,
  query = MyData,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)


saveRDS (MyData, file="Combined_SCT_MappedPBMC.Rds")

p1 = DimPlot(MyData, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(MyData, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
pdf(file = "DIMPLOT.pdf")
p1 + p2
dev.off()
```
***Transfering/Labeling identies to the Query dataset***
```
Idents(MyData) <- 'predicted.celltype.l1'
```
***Creating Count table per sample based on Celltypes***

```
CountTable<-table(MyData@meta.data$predicted.celltype.l1, MyData@meta.data$SAMP)
write.table(CountTable,file="CountTable.txt",sep="\t")
```
