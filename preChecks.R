# Pre-analysis checks and supplemental figures for 
# Leveraging transcriptional signatures of diverse stressors for bumble bee conservation
## Authors: Gabriela M. Quinlan, Heather M. Hines, Christina M. Grozinger 

# Load the data and libraries 

rm(list=ls())

for (package in c("BiocManager", "DESeq2", "dplyr", "stringr", 
                  "ggplot2", "genefilter", "tidyr", "plotrix",
                  "patchwork", "ggVennDiagram", "randomForest",
                  "ggthemes", "ggh4x", "forcats", "UpSetR", "readxl",
                  "topGO", "GO.db", "AnnotationDbi", "pheatmap"
)) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


# Upload data ----
## lab data 
sampleLab <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS1", skip = 2))
countsLab <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS3", skip = 3))


# Format read-in data
rownames(countsLab) <- countsLab[,1]
rownames(sampleLab) <- sampleLab$sample_name

sampleLab <-  sampleLab %>%
  mutate(trt = factor(trt, levels= c("starve", "cold", "heat", "control", "sham", "ecoli"))) %>%
  mutate(trt = recode(trt, "control" = "Control","cold" = "Cold", "ecoli" = "Immune", 
                      "heat" = "Heat", "sham" = "Sham","starve" = "Starve")) %>%
  mutate(tissue = recode(tissue, "abdomen" = "Abdomen", "thorax" = "Thorax")) %>%
  mutate(colony = as.factor(colony), 
         trt = relevel(as.factor(trt), ref = "Control"),
         round = as.factor(round), 
         date = as.factor(date),
         tissue = as.factor(tissue)) 

countsLab <- countsLab[, -which(names(countsLab) %in% c("...1"))]

# separate tissues 
sampleLabTho <- sampleLab[which(sampleLab$tissue == "Thorax"),]
sampleLabAbd <- sampleLab[which(sampleLab$tissue == "Abdomen"),]

countTho <- countsLab[,grepl("THO", names(countsLab))]
countAbd <- countsLab[,grepl("ABD", names(countsLab))]

# Run Lab DESeq -----

dds <- DESeqDataSetFromMatrix(countData=countsLab, 
                              colData=sampleLab, 
                              design=~colony + tissue + date +  trt, tidy = F)

ddsTho <- DESeqDataSetFromMatrix(countData=countTho, 
                                 colData=sampleLabTho, 
                                 design=~colony + date +  trt, tidy = F)

ddsAbd <- DESeqDataSetFromMatrix(countData=countAbd, 
                                 colData=sampleLabAbd, 
                                 design=~colony + date +  trt, tidy = F)

# Pre-filter (remove rows with no or low counts)
keep_genes <- rowSums(counts(dds)) > 10
dds <- dds[ keep_genes, ]

keep_genesAbd <- rowSums(counts(ddsAbd)) > 10
ddsAbd <- ddsAbd[ keep_genesAbd, ]

keep_genesTho <- rowSums(counts(ddsTho)) > 10
ddsTho <- ddsTho[ keep_genesTho, ]

vsdata <- vst(dds, blind=FALSE)
vsdataAbd <- vst(ddsAbd, blind=FALSE)
vsdataTho <- vst(ddsTho, blind=FALSE)

# Sample distance -- as a check, visually check distances between samples  
## Heat map ----
sampleDists <- dist(t(assay(vsdata)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsdata$trt, vsdata$tissue, sep="-" )
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

sampleDistsTho <- dist(t(assay(vsdataTho)))
sampleDistMatrixTho <- as.matrix( sampleDistsTho )
rownames(sampleDistMatrixTho) <- paste( vsdataTho$trt, vsdataTho$colony, sep="-" )
colnames(sampleDistMatrixTho) <- NULL
pheatmap(sampleDistMatrixTho,
         clustering_distance_rows=sampleDistsTho,
         clustering_distance_cols=sampleDistsTho)

sampleDistsAbd <- dist(t(assay(vsdataAbd)))
sampleDistMatrixAbd <- as.matrix( sampleDistsAbd )
rownames(sampleDistMatrixAbd) <- paste( vsdataAbd$trt, vsdataAbd$colony,vsdataAbd$round, sep="-" )
colnames(sampleDistMatrixAbd) <- NULL
pheatmap(sampleDistMatrixAbd,
         clustering_distance_rows=sampleDistsAbd,
         clustering_distance_cols=sampleDistsAbd)


## PCA ----
pcaData <- plotPCA(vsdata, intgroup=c("trt", "tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
A <- ggplot(pcaData, aes(PC1, PC2, color=trt, shape=tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  guides(color=guide_legend(title="Trt"),
         shape=guide_legend(title="Tissue"), 
         linetype=guide_legend(title="Tissue")) +
  stat_ellipse(aes(linetype = tissue) ) + 
  theme_bw(base_size = 15) + 
  scale_colour_manual(breaks = c("Starve", "Cold", "Heat", "Control","Sham", "Immune"), values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7")) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

pcaDataAbd <- plotPCA(vsdataAbd, intgroup=c("trt", "colony"), returnData=TRUE)
percentVarAbd <- round(100 * attr(pcaDataAbd, "percentVar"))
B <- ggplot(pcaDataAbd, aes(PC1, PC2, color=trt, shape=colony, group=trt)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarAbd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarAbd[2],"% variance")) + 
  coord_fixed() + 
  guides(color=guide_legend(title="Trt"),
         shape=guide_legend(title="Colony")) +
  scale_colour_manual(breaks = c("Starve", "Cold", "Heat", "Control","Sham", "Immune"), values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7")) + 
  stat_ellipse() + 
  theme_bw(base_size = 15) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))


pcaDataTho <- plotPCA(vsdataTho, intgroup=c("trt", "colony"), returnData=TRUE)
percentVarTho <- round(100 * attr(pcaDataTho, "percentVar"))
C <- ggplot(pcaDataTho, aes(PC1, PC2, color=trt, shape=colony, group=trt)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarTho[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarTho[2],"% variance")) + 
  coord_fixed() + 
  guides(color=guide_legend(title="Trt"),
         shape=guide_legend(title="Colony")) +
  scale_colour_manual(breaks = c("Starve", "Cold", "Heat", "Control","Sham", "Immune"), values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7")) + 
  stat_ellipse() + 
  theme_bw(base_size = 15) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

(A / (C + B)) + plot_annotation(tag_levels = 'A') + 
  plot_layout(guides='collect') + 
  plot_layout(widths = c(2, 2, 2))
#ggsave("PCAall.tiff", device="tiff", width=10, height=8)

