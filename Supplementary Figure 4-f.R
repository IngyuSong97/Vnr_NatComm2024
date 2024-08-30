##### Supplementary Figure 4-f #####
# Levels of prd, GT28s, lipases, and DAG metabolism genes in male reproductive organ 
# 
# Single-nucleus RNA-seq data are retrieved from Fly Cell Atlas (https://flycellatlas.org).
# Four files (two loom files and two meta data csv file) are used in this script;
# 1. s_fca_biohub_testis_10x_ss2.loom
# 2. metadata.testis.csv (file name is manually defined)
# 3. s_fca_biohub_male_reproductive_glands_10x_ss2.loom
# 4. metadata.mrg.csv (file name is manually defined)
# 
# Figure generated in this script is further editted manually.


##### Required packages #####
library(Seurat)
library(dplyr)
library(ggplot2)
library(SCopeLoomR)
library(SeuratDisk)


##### Data loading and processing #####
# Testis data
# Loading Loom file
scd = open_loom("Directory_To_Data/s_fca_biohub_testis_10x_ss2.loom", mode="r+")
scd = get_dgem(scd)
scd = data.frame(gene = rownames(scd), scd, row.names = rownames(scd))
scd = scd[,-1]
build_loom(file.name = "testis.loom",
           dgem = scd,
           title = "testis",
           genome = "Dmel")
loom = Connect(filename = "testis.loom", mode = "r+")

# Loom to Seurat conversion
seurat = as.Seurat(loom)

# Metadata integration
celltypes = read.delim("Directory_To_Data/metadata.testis.csv", sep = ',', col.names = c('cellID', 'annotation'))
celltypes$cellID = gsub('-', '.', celltypes$cellID) # fitting cell IDs
meta = seurat@meta.data
meta$cellID = rownames(meta)
meta = left_join(meta, celltypes, by='cellID')
rownames(meta) = meta$cellID
meta = meta[,-6]
seurat@meta.data = meta
Idents(seurat) = "annotation"

# Data processing
seurat = NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat = ScaleData(seurat, features = rownames(seurat))
seurat.testis = seurat


# Male reproductive glands data
# Loading Loom file
scd = open_loom("Directory_To_Data/s_fca_biohub_male_reproductive_glands_10x_ss2.loom", mode="r+")
scd = get_dgem(scd)
scd = data.frame(gene = rownames(scd), scd, row.names = rownames(scd))
scd = scd[,-1]
build_loom(file.name = "mrg.loom",
           dgem = scd,
           title = "MRG",
           genome = "Dmel")
loom = Connect(filename = "mrg.loom", mode = "r+")

# Loom to Seurat conversion
seurat = as.Seurat(loom)

# Metadata integration
celltypes = read.delim("Directory_To_Data/metadata.mrg.csv", sep = ',', col.names = c('cellID', 'annotation'))
celltypes$cellID = gsub('-', '.', celltypes$cellID) # fitting cell IDs
meta = seurat@meta.data
meta$cellID = rownames(meta)
meta = left_join(meta, celltypes, by='cellID')
rownames(meta) = meta$cellID
meta = meta[,-6]
seurat@meta.data = meta
Idents(seurat) = "annotation"

# data processing
seurat = NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat = ScaleData(seurat, features = rownames(seurat))
seurat.mrg = seurat


# Merging two data
seurat = merge(seurat.testis, seurat.mrg)


##### Vizualization #####
# Melting cell types for simpler visualization
seurat@meta.data$annotation = gsub('early-mid elongation-stage ', '', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('early elongation stage ', '', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('mid-late elongation-stage ', '', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('cell 2', 'cell', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('cell 1', 'cell', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('branch a', 'branch', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub('branch b', 'branch', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub(' \\(2\\)', '', seurat@meta.data$annotation)
seurat@meta.data$annotation = gsub(' \\(1\\)', '', seurat@meta.data$annotation)
seurat = subset(seurat, idents = unique(Idents(seurat))[-c(3,9,29)])

# Order
order = c("late primary spermatocyte", "spermatocyte 7a", "spermatocyte", "spermatid",
          "adult tracheal cell", "hemocyte", "cyst stem cell", "early cyst cell",
          "head cyst cell", "cyst cell intermediate", "cyst cell branch",
          "spermatocyte cyst cell branch", "late cyst cell branch", "spermatocyte 2",
          "spermatocyte 1", "spermatocyte 0", "mid-late proliferating spermatogonia",
          "spermatogonium-spermatocyte transition", "spermatogonium",
          "secretory cell of the male reproductive tract", "seminal vesicle",
          "male gonad associated epithelium", "spermatocyte 3", "spermatocyte 4",
          "spermatocyte 5", "spermatocyte 6", "pigment cell", "epithelial cell",
          "male reproductive tract muscle", "anterior ejaculatory duct", "ejaculatory bulb",
          "ejaculatory bulb epithelium", "male accessory gland main cell",
          "male accessory gland secondary cell")
Idents(seurat) = factor(Idents(seurat), levels = unique(rev(order)))

# Target genes
MGDG = c('Ugt305A1','Ugt302C1','Ugt37D1','Ugt304A1','Ugt35E2','CG14512')
PLA = c('GXIVsPLA2','iPLA2-VIA','GIIIspla2','CG1309','nSMase', 'JMJD7',
        'PAPLA1','lama','CG14034','CG42237','CG14507','CG5966','sPLA2', 'CG3009')
Hit = c('Agpat1','Gk1','mino','inaE','Lsd-2')
target.genes = c('prd', MGDG, PLA, Hit)

# Plotting
p = DotPlot(seurat, features = target.genes) +
  labs(x='genes', y='cell types') +
  theme(axis.text.x=element_text(angle=90, hjust=1))
p