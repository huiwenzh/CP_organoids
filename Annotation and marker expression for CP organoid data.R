# using 40%, 12000 cells as input
# 12/02/2024
library(data.table)
library(Seurat)
sample1 <- read.table('sample1_20240212125241_TCM.tsv.gz',header = T,row.names = 1)
sample2 <- read.table('sample2_20240212142138_TCM.tsv.gz',header = T,row.names = 1)

all.genes <- intersect(rownames(sample1), rownames(sample2))
sample1.1 <- sample1[match(all.genes, rownames(sample1)),]
colnames(sample1.1) <- paste0(colnames(sample1.1),'_sample1')
sample2.1 <- sample2[match(all.genes, rownames(sample2)),]
colnames(sample2.1) <- paste0(colnames(sample2.1),'_sample2')
all <- cbind(sample1.1, sample2.1) 
meta <- data.frame(sample=c(rep('sample1',12000),rep('sample2',12000)))
rownames(meta) <- colnames(all)

# creating seurat object
all.obj <- CreateSeuratObject(all,min.cells = 3,min.features = 200,meta.data = meta)  
all.obj <- NormalizeData(all.obj)
all.obj <- FindVariableFeatures(all.obj, selection.method = "vst", nfeatures = 3000)
all.obj <- ScaleData(all.obj)
all.obj <- RunPCA(all.obj)
 
all.obj <- FindNeighbors(all.obj, dims = 1:20)
all.obj <- FindClusters(all.obj, resolution = 0.8)
 
# performing doublet removal analysis
library(DoubletFinder)
sweep.res.list1<- paramSweep(all.obj, PCs = 1:20, sct = FALSE)
sweep.stats1 <- summarizeSweep(sweep.res.list1 ,GT = FALSE)
bcmvn <- find.pK(sweep.stats1) # 0.01
homotypic.prop <- modelHomotypic(all.obj$RNA_snn_res.0.8)           
nExp_poi <- round(0.075*nrow(all.obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
all.obj <- doubletFinder(all.obj, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
all.obj <- doubletFinder(all.obj, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_1286", sct = FALSE)

#remove doublets 
all.sub <- subset(all.obj,DF.classifications_0.25_0.01_1061=='Singlet')
all.sub.dat <- all.sub[['RNA']]$counts
# remove overly MT genes
all.sub1 <- CreateSeuratObject(all.sub.dat,meta.data = all.sub@meta.data, min.cells = 10, min.features = 500)
all.sub1[["percent.mt"]] <- PercentageFeatureSet(all.sub1, pattern = "^MT-") 
VlnPlot(all.sub1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
all.sub1 <- subset(all.sub1, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20)
# Visualize QC metrics as a violin plot
all.sub1 <- NormalizeData(all.sub1)
all.sub1 <- FindVariableFeatures(all.sub1, selection.method = "vst", nfeatures = 2000)
all.sub1 <- ScaleData(all.sub1)
all.sub1 <- RunPCA(all.sub1)
all.sub1 <- FindNeighbors(all.sub1, dims = 1:30)
all.sub1 <- FindClusters(all.sub1, resolution = 0.6)
all.sub1 <- RunUMAP(all.sub1, dims = 1:30 )
 
DimPlot(all.sub1, reduction = "umap",label = T)
DimPlot(all.sub1, reduction = "umap",split.by = 'sample')
FeaturePlot(all.sub1,features = c('TTR','CP','KL','LMX1A','AQP1','MSX1','MSX2'))

#update marker 
Neurons <- c("RBFOX3" , "NEFM" ,  "SYP"   , "RELN","CUX1","BCL11B", "TBR1","FOXP2","NEFH"  ) #"MAP2" , "STAB2",
Progenitors <- c("SOX2" ,"NES","NOTCH1","HES1" ,"OCLN" ,"PAX6","CDH1" ,"FABP7","HOPX" ) # 
## sctype for first round annotation
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = all.sub1[['RNA']]$scale.data, scaled = TRUE,  gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(all.sub1@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(all.sub1@meta.data[all.sub1@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(all.sub1@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
sctype_scores$type[4] <- 'Choroid plexus'
sctype_scores$type[1] <- 'Choroid plexus'
sctype_scores$type[7] <- 'Cortical Hem'
sctype_scores$type[11]<- 'Neuroepithelial cells' # cancer cells
all.sub1@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  all.sub1@meta.data$customclassif[all.sub1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(all.sub1, reduction = "umap")
DimPlot(all.sub1, reduction = "umap",split.by = 'sample')
DimPlot(all.sub1, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        
DotPlot(all.sub1,features = c('TTR','KL','LMX1A','AQP1','MSX1','MSX2'),group.by = 'customclassif')

# find clusters markers for manual annotations
Cluster_marker = FindAllMarkers(all.sub1,only.pos=T,logfc.threshold=0.25, min.pct = 0.05)
Cluster_marker_main <- Cluster_marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)%>%
  slice_head(n = 50)

# visulisation for figures in main and supplementary
library(scCustomize)
library(ggplot2)
all.sub1$Predicted_label <- all.sub1$customclassif
FeaturePlot_scCustom(seurat_object = all.sub1, features = c('TTR','KL','LMX1A','AQP1','MSX1','MSX2'),keep.scale='all')
FeaturePlot_scCustom(seurat_object = all.sub1, features =Progenitors)#Neurons
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Predicted_label', figure_plot = TRUE,label = T )
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Predicted_label', figure_plot = TRUE,label = T,split.by = 'sample')
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='seurat_clusters', figure_plot = TRUE,label = T )
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='seurat_clusters', figure_plot = TRUE,label = T ,split.by = 'sample')

# Identify COVID genes
COVID <- c('ACE2','TMPRSS2','FURIN','NRP1','NRP2','CTSL','CTSB')
FeaturePlot_scCustom(seurat_object = all.sub1, features =COVID)
FeaturePlot_scCustom(seurat_object = all.sub1, features ='CD68')

#keep same scale for plotting marker expression
library(RColorBrewer)
FeaturePlot(all.sub1, features = c('TTR','KL','LMX1A','AQP1','MSX1','MSX2'),keep.scale='all',cols = c('#f0e690','#41049d'),order = T) #&  scale_colour_viridis_c(option = 'plasma') #CP

FeaturePlot(all.sub1, features = COVID,keep.scale='all',cols = c('#f0e690','#41049d'),order = T)

FeaturePlot(all.sub1, features = Progenitors,keep.scale='all',cols = c('#f0e690','#41049d'),order = T)

FeaturePlot(all.sub1, features = Neurons,keep.scale='all',cols = c('#f0e690','#41049d'),order = T)



Stacked_VlnPlot(all.sub1,features = c('TTR','AQP1','MSX1','MAP2','NEFM','SOX2','PAX6'),group.by = 'Predicted_label')

DotPlot(all.sub1,features = c(ct_marker,'MAP2'))
ct_marker <- c('S100B', 'TOP2A','AQP1','MSX1','UNC5B','RSPO3' ,'PAX6', 'GAD1', 'GAD2', 'SCN2A', 'SLC17A6', 'NEFM' ,'SOX2','PDGFRB','TTR')#'TTR',
# more specific neuron types
#sctype_scores$type[3] <- 'GABAergic neurons'


# final annotation 
all.sub1$Annotation <- all.sub1$seurat_clusters
levels(all.sub1$Annotation) <- c('Choroid plexus','Cortical hem','Choroid plexus','SOX2+ RG','Inhibitory neurons','Excitatory neurons','Excitatory neurons','Excitatory neurons','SOX2+ RG','Astro-progenitor','Unknown','HOPX+ RG','Ependymal','Proliferating NSC')
library(forcats)
fct_relevel(all.sub1$Annotation,'Cortical hem','Ependymal','Choroid plexus','Proliferating NSC','SOX2+ RG','HOPX+ RG','Astro-progenitor','Excitatory neurons','Inhibitory neurons','Unknown')
all.sub1$Annotation <- ordered(all.sub1$Annotation,levels  = c('Cortical hem','Ependymal','Choroid plexus','Proliferating NSC','SOX2+ RG','HOPX+ RG','Astro-progenitor','Excitatory neurons','Inhibitory neurons','Unknown'))

DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Annotation', figure_plot = TRUE,label = T )
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Annotation', figure_plot = TRUE,label = T ,split.by = 'sample')
Idents(all.sub1) <- 'Annotation'
DotPlot_scCustom(all.sub1,features = c('MSX1','AQP1','MCIDAS', 'TOP2A','SOX2','HOPX', 'SLC1A3','SLC17A6','GAD1',dot.min=20),group.by = 'Annotation') 


DotPlot_scCustom(all.sub1,features = c('MSX1','MCIDAS','AQP1','TOP2A','SOX2','PAX6','HOPX', 'SLC1A3','NEFL','SLC17A6','GAD1'),x_lab_rotate = T,group.by = 'Annotation',col.min=0.3,scale.min=0.2) #,'TTR'

Stacked_VlnPlot(all.sub1,features = c('MSX1','AQP1','MCIDAS','TOP2A','SOX2','PAX6','HOPX', 'SLC1A3','SLC17A6','GAD1')) #,'TTR'
ct_color <- c('#aae479','#64c987','#468c5e','#87b3da','#3c84c3','#485c96','#e9c848','#f95e00','#89290a','#545b59')
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Annotation', figure_plot = TRUE,colors_use = ct_color)
DimPlot_scCustom(seurat_object = all.sub1,  group.by ='Annotation', figure_plot = TRUE,colors_use = ct_color,split.by = 'sample')
FeaturePlot_scCustom(seurat_object = all.sub1, features =COVID)
VlnPlot_scCustom(all.sub1, features =COVID, colors_use =  ct_color, num_columns = 4)
D
VlnPlot_scCustom(all.sub1, features =c('NRP1','NRP2','CTSL','CTSB'), colors_use =  ct_color, pt.size = 0.01,num_columns = 2)
VlnPlot_scCustom(all.sub1, features =c('TTR','MAP2'), colors_use =  ct_color, pt.size = 0.01,num_columns = 2,plot_boxplot = T)

FeaturePlot(all.sub1,c('CLDN5','F11R',"TJP1","CLDN1",'CLDN2','CLDN3','OCLN'))
VlnPlot_scCustom(all.sub1, features =c('CLDN5','F11R',"TJP1","CLDN1",'CLDN2','CLDN3','OCLN'), colors_use =  ct_color, pt.size = 0.01,num_columns = 4)
VlnPlot(all.sub1, features =c('CLDN5','F11R',"TJP1","CLDN1",'CLDN2','CLDN3','OCLN'), stack = T)

DotPlot_scCustom(all.sub1, features =c('TTR','KL','LMX1A','AQP1','MSX1','MSX2'))
VlnPlot_scCustom(all.sub1, features =c('TTR','KL','LMX1A','AQP1','MSX1','MSX2'), colors_use =  ct_color, pt.size = 0.01,num_columns = 2)


DotPlot_scCustom(all.sub1, features =Neurons,x_lab_rotate = T)
DotPlot_scCustom(all.sub1, features =Progenitors,x_lab_rotate = T)
DotPlot_scCustom(all.sub1, features =c(Progenitors,Neurons),x_lab_rotate = T)
