##### Supplementary Figure 2-a #####
# Transcriptional changes of six GT28 genes in male accessory glands (MAG)
# at 0, 1, and 2 days after eclosion
# 
# MAG bulk RNA sequencing data is deposited in NCBI SRA (PRJNA1140569).
# RNA-seq row data was process using kallisto to generate 'MAG_day0-3.txt' file.
# fbgn_fbtr_fbpp_expanded_fb_2022_03.tsv was download from FlyBase. (https://flybase.org/)
# 
# Figure generated in this script was further edited manually.


##### Required packages #####
require(dplyr)
require(gplots)

##### Data loading and processing #####
# Loading tpm values
dt = read.delim('Directory_To_Data/MAG_day0-3.txt', row.names = 1)
dt = dt[,grep('tpm', colnames(dt))]
colnames(dt) = gsub('_tpm', '', colnames(dt))

# Merging transcripts into genes
fbid = unique(read.delim('Directory_To_Data/fbgn_fbtr_fbpp_expanded_fb_2022_03.tsv', skip = 4)[,c(3,4,8)])
colnames(fbid) = c('FBgn','symbol','FBtr')
fbgnid = unique(fbid$FBgn)
gn2tr = lapply(fbgnid, function(x){
  fbid[grep(x,fbid$FBgn),3]
})
names(gn2tr) = fbgnid
trSum2gn = lapply(gn2tr, function(x){
  colSums(dt[which(rownames(dt) %in% x),])
})
dt = t(as.data.frame(trSum2gn))

# FBgnIDs to gene symbols
dt = data.frame(FBgn = rownames(dt), dt)
dt = left_join(unique(fbid[,1:2]), dt, by='FBgn')[,-1]
rownames(dt) = dt$symbol
dt = dt[,-1]

# GT28 genes
GT28 = c('Ugt305A1', 'Ugt302C1', 'Ugt35E2', 'Alg13', 'Ugt304A1', 'Ugt37D1')
GT28 = dt[GT28,]

# Mean value
GT28 = data.frame(row.names = rownames(GT28), day0 = rowMeans(GT28[,1:3]),
                  day1 = rowMeans(GT28[,4:6]), day2 = rowMeans(GT28[,7:9]))

# Scaling
GT28 = t(scale(t(GT28)))


##### Visualization #####
p = heatmap.2(GT28, trace="none",
              distfun = function(x){dist(x, method = 'euclidean')},
              hclustfun = function(x){hclust(x, method = 'average')},
              cexRow = 1.25, cexCol = 1.25, sepcolor = "black", key = T,
              key.title = "scaled expression",key.ylab = "",key.xlab = "", lhei = c(2,6),
              rowsep = c(0,dim(GT28)[1]), colsep = c(0,3,6,dim(GT28)[2]),
              margins = c(5, 10), srtCol = 90, sepwidth = c(0.005, 0.005),
              col=colorRampPalette(c("navy","white","red"))(256),
              scale = "none", density.info = "none",
              Rowv = T, Colv = F, dendrogram = 'row')