# Leveraging transcriptional signatures of diverse stressors for bumble bee conservation
## Authors: Gabriela M. Quinlan, Heather M. Hines, Christina M. Grozinger 

# Load the data and libraries 

rm(list=ls())

for (package in c("BiocManager", "DESeq2", "dplyr", "stringr", 
                  "ggplot2", "genefilter", "tidyr", "plotrix",
                  "patchwork", "ggVennDiagram", "randomForest",
                  "ggthemes", "ggh4x", "forcats", "UpSetR", 
                  "topGO", "GO.db", "AnnotationDbi", "readxl"
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

## field data 
sampleField <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS2", skip = 2))
countsField <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS4", skip = 2))

# GO databases
bimp_to_dmel <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS11", skip = 3))

# read in table as short cut (code commented out below for how to produce from bimp_to_dmel)
Custom_GENE2GO <- ViSEAGO::Custom2GO("custom.txt") # make custom GO from bimp_to_dmel (code in myGO)


# Format read-in data
rownames(countsLab) <- countsLab[,1]
rownames(sampleLab) <- sampleLab$sample_name

sampleLab <-  sampleLab %>%
  mutate(colony = as.factor(colony), 
                      trt = relevel(as.factor(trt), ref = "control"),
                      round = as.factor(round), 
                      date = as.factor(date),
                      tissue = as.factor(tissue)) 

countsLab <- countsLab[, -which(names(countsLab) %in% c("...1"))]

# Remove very different sample (see preprocessing script)
countsLab <- countsLab[,-which(names(countsLab) %in% c("GQ93ABD"))] # remove outlier sample"GQ93ABD"
sampleLab <- sampleLab[which(sampleLab$sample_name != "GQ93ABD"),]

# separate tissues 
sampleLabTho <- sampleLab[which(sampleLab$tissue == "thorax"),]
sampleLabAbd <- sampleLab[which(sampleLab$tissue == "abdomen"),]

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

dds <- DESeq(dds)
ddsAbd <- DESeq(ddsAbd)
ddsTho <- DESeq(ddsTho)

# Summarise Lab Contrasts ----

deseqSum <- function (dat, control, stress, lfc, p) {
  df <- as.data.frame(results(dat, contrast = c("trt", stress, control),
                              lfcThreshold = lfc, alpha =p, altHypothesis="greaterAbs"))
  df$comp <- paste(control, stress, sep="-")
  df$tiss <- recode(deparse(substitute(dat)), "dds" = "all", "ddsAbd" = "abd", "ddsTho" = "tho")
  df$gene <- rownames(df)
  return(df)
}

# stringent
stringent <- rbind(deseqSum(dat = dds, control= "control", stress= "ecoli", lfc=0.5, p=0.05),
                   deseqSum(dat = dds, control= "sham", stress= "ecoli", lfc=0.5, p=0.05), 
                   deseqSum(dat = dds, control= "control", stress= "sham", lfc=0.5, p=0.05),
                   deseqSum(dat = dds, control= "control", stress= "heat", lfc=0.5, p=0.05), 
                   deseqSum(dat = dds, control= "control", stress= "cold", lfc=0.5, p=0.05), 
                   deseqSum(dat = dds, control= "control", stress= "starve", lfc=0.5, p=0.05), 
                   
                   deseqSum(dat = ddsTho, control= "control", stress= "ecoli", lfc=0.5, p=0.05),
                   deseqSum(dat = ddsTho, control= "sham", stress= "ecoli", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsTho, control= "control", stress= "sham", lfc=0.5, p=0.05),
                   deseqSum(dat = ddsTho, control= "control", stress= "heat", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsTho, control= "control", stress= "cold", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsTho, control= "control", stress= "starve", lfc=0.5, p=0.05), 
                   
                   deseqSum(dat = ddsAbd, control= "control", stress= "ecoli", lfc=0.5, p=0.05),
                   deseqSum(dat = ddsAbd, control= "sham", stress= "ecoli", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsAbd, control= "control", stress= "sham", lfc=0.5, p=0.05),
                   deseqSum(dat = ddsAbd, control= "control", stress= "heat", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsAbd, control= "control", stress= "cold", lfc=0.5, p=0.05), 
                   deseqSum(dat = ddsAbd, control= "control", stress= "starve", lfc=0.5, p=0.05))
rownames(stringent) <- 1:nrow(stringent)

# DEG trt compare ----

sigGeneList<- list(
  Sham = stringent[which(stringent$comp == "control-sham" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Immune = stringent[which(stringent$comp == "control-ecoli" & stringent$tiss == "tho" &  stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Heat = stringent[which(stringent$comp == "control-heat" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Cold = stringent[which(stringent$comp == "control-cold" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Starve = stringent[which(stringent$comp == "control-starve" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
)

common <- unique(c(intersect(sigGeneList$Sham, sigGeneList$Immune),
                   intersect(sigGeneList$Heat, sigGeneList$Cold), 
                   intersect(sigGeneList$Heat, sigGeneList$Immune),
                   intersect(sigGeneList$Immune, sigGeneList$Cold), 
                   intersect(sigGeneList$Heat, sigGeneList$Sham), 
                   intersect(sigGeneList$Sham, sigGeneList$Cold)))


shared <- stringent %>% 
  dplyr::filter(gene %in% common) %>%
  filter(padj < 0.05 &
           abs(log2FoldChange) > 0.5 & 
           tiss == "tho") %>%
  dplyr::select(comp, gene, log2FoldChange) %>%
  pivot_wider(names_from = comp, values_from = log2FoldChange) %>% 
  dplyr::select("gene", "control-heat", "control-cold", "control-ecoli", "control-sham") %>%
  rowwise() %>%
  mutate(mean = mean(c_across(`control-heat`:`control-sham` ), na.rm=T))

shared$nshare <- rowSums(!is.na(shared))

#write.csv(shared, "shared.csv")


topDEGTab <- stringent[which(stringent$tiss == "tho" &
                               stringent$padj < 0.05 &
                               abs(stringent$log2FoldChange) > 0.5), ] %>% 
  mutate(absLFC = abs(log2FoldChange)) %>%
  arrange(-absLFC) %>% 
  group_by(comp) %>%
  top_n(n = 15) %>% View()

#write.csv(topDEGTab, "topDEGtab.csv")


# DEG GO ----
# Commented out creation of custom DF from dmel_custom.txt (saved and uploaded) as Custom_GENE2GO
#dmel_go <- read.table("gene_association.fb", # read D. molanigaster fly base
#                     sep="\t",
#                    comment.char="!",
#                   na.strings=".",
#                  stringsAsFactors=FALSE,
#                 quote="", fill=FALSE)[,c(1,2,3,5,7)]

bimp_to_dmel <- bimp_to_dmel %>%
  rename(bimp.gene = biRef, 
         dmel.gene = dmFB)

# format dmel GO 
#names(dmel_go) <- c("taxid","gene_id","gene_symbol","GOID","evidence") # headers
#dmel_go$gene_symbol <- dmel_go$gene_id # idk why you'd do this. you already have this column 
#dmel_go <- dmel_go[dmel_go$GOID %in% keys(GO.db),] # Get only those w/ annotation; annotation maps (from Bioconductor)
#dmel_go <- dmel_go[!dmel_go$evidence=="ND",] # exclude when function unknown
#dmel_go <- dmel_go[complete.cases(dmel_go), ] # complete cases 
# dropped ~7.5k out of ~130k

#dmel_go <- dmel_go[dmel_go$gene_id %in% bimp_to_dmel$dmel.gene,] # include only those annontated D. molanigaster that also appear in A. mellifera 
# only about 90k
#write.table(dmel_go, "dmel_custom.txt",
#           sep="\t",quote=F,row.names=F) # save this (skip some steps) needs to be a table for next function. 

# with update -- some not in database
# newDat <- dmel_go[which(dmel_go$GOID %in% as.vector(keys(GO.db))),]
# write.table(newDat,"custom.txt",sep="\t",row.names=FALSE)

Custom_GENE2GO <- ViSEAGO::annotate("FB",Custom_GENE2GO)

# deg results for each tissue/ contrast 


# Run the GO enrichment analysis (function from Sean) 
GOEA <- function(DEG.object, contrast, tissue){
  geneList0 <- bimp_to_dmel$dmel.gene # list of D. molanigaster genes with B. impatiens analougs
  myInterestingGenes <- c(DEG.object$gene) # List of down and upregulated genes
  myInterestingGenes <- bimp_to_dmel[bimp_to_dmel$bimp.gene %in% myInterestingGenes, 2] # of the interesting genes, get FB equivilants
  geneList <- factor(as.integer(geneList0 %in% myInterestingGenes)) # which are interesting (call only those)
  names(geneList) <- geneList0
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, # specify biological 
                annot = annFUN.gene2GO, gene2GO = Custom_GENE2GO@BP)
  test <- runTest(GOdata, algorithm ="elim", statistic = "fisher")
  results <- GenTable(GOdata,elim=test,topNodes=length(test@score))
  results <- suppressWarnings(results[!is.na(as.numeric(as.character(results$elim))),]) # remove NA for elim 
  results <- suppressWarnings(results[!is.na(results$elim),]) # same but I guess just dif. format
  results$elim <- as.numeric(results$elim)
  results$percent.coverage <- results$Significant/results$Annotated
  results.sig <- results[results$elim < 0.05, ] # significance threshold 
  rank <- list()
  for(i in 1:length(results.sig$GO.ID)){
    rank[i] <- length(unlist(as.list(GOBPCHILDREN[results.sig$GO.ID[i]]))) # annotation to biological process
  }
  results.sig$rank <- unlist(rank)
  results.sig$contrast <- contrast # you supply these in the function
  results.sig$tissue <- tissue # you supply these in the function
  return(results.sig)
}

# comment out contrasts with no significant GO terms 

thoCvsCdGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "control-cold"),], "CvCd", "tho")
thoCvEGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "control-ecoli"),], "CvE", "tho")
thoCvHGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "control-heat"),], "CvH", "tho")
thoCvShGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "control-sham"),], "CvSh", "tho")
#thoCvStGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "control-starve"),], "CvSt", "tho")
thoSvEGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "tho" & stringent$comp == "sham-ecoli"),], "ShvE", "tho")

#abdCvsCdGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "control-cold"),], "CvCd", "abd")
abdCvEGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "control-ecoli"),], "CvE", "abd")
abdCvHGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "control-heat"),], "CvH", "abd")
abdCvShGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "control-sham"),], "CvSh", "abd")
#abdCvStGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "control-starve"),], "CvSt", "abd")
abdSvEGO <- GOEA(stringent[which(stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5 & stringent$tiss == "abd" & stringent$comp == "sham-ecoli"),], "ShvE", "abd")

## Format tables ----
percent <- function(x, digits = 2, format = "f", ...){
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

# all results (comment out if none)
goAll <- rbind(thoCvsCdGO, thoCvEGO, thoCvHGO, thoCvShGO, thoSvEGO, #thoCvStGO,
               abdCvEGO, abdCvHGO, abdCvShGO, abdSvEGO) #abdCvsCdGO, abdCvStGO

#separate table for each 
tab4 <- goAll[goAll$tissue=="abd",]
tab3 <- goAll[goAll$tissue=="tho",]

# how many total genes up + down regulated in each tissue
DEGs_A <- length(unique(c(#abdCvsCdGO[[1]],abdCvsCdGO[[2]],
  abdCvEGO[[1]],abdCvEGO[[2]],
  abdCvHGO[[1]],abdCvHGO[[2]], 
  abdCvShGO[[1]],abdCvShGO[[2]], 
abdSvEGO[[1]],abdSvEGO[[2]])))
#abdCvStGO[[1]],abdCvStGO[[2]]))) #
DEGs_T <- length(unique(c(thoCvsCdGO[[1]],thoCvsCdGO[[2]],
                          thoCvEGO[[1]],thoCvEGO[[2]],
                          thoCvHGO[[1]],thoCvHGO[[2]], 
                          thoCvShGO[[1]],thoCvShGO[[2]],
                          thoSvEGO[[1]],thoSvEGO[[2]])))#,
#thoCvStGO[[1]],thoCvStGO[[2]])))

# only those with were significant p < 0.05
tab3 <- tab3[tab3$elim<0.05,]
tab3 <- tab3[tab3$Annotated<(DEGs_T*0.5),] # choose small annotation values (less than 1/4 total DEG (up or down)) I changed this to 50% 
# small values specify the annotation is SPECIFIC (describes number of genes that are annotated for each GO)
tab3 <- tab3[order(-tab3$Annotated),] # decending order by annotated (lg at top)
row.names(tab3) <- NULL
tab3 <- tab3[1:20,]
tab3$percent.coverage <- percent(tab3$percent.coverage)
tab3$Coverage <- paste(paste(tab3$Significant,
                             tab3$Annotated,sep="/"),
                       paste("(",tab3$percent.coverage,")",sep=""),sep=" ")
#tab3 <- tab3[,c(9,1,2,11,6)]
#names(tab3) <- c("Contrast","GO.ID","Term","Coverage","p")
#tab3 <- tab3[order(-tab3$Contrast),]
tab3$p <- formatC(tab3$elim, format = "e", digits = 2) #tab3$p

tab4 <- tab4[tab4$elim<0.05,]
tab4 <- tab4[tab4$Annotated<(DEGs_A*0.5),] # I changed this to 50% 
tab4 <- tab4[order(-tab4$Annotated),]
row.names(tab4) <- NULL
tab4 <- tab4[1:20,]
tab4$percent.coverage <- percent(tab4$percent.coverage)
tab4$Coverage <- paste(paste(tab4$Significant,
                             tab4$Annotated,sep="/"),
                       paste("(",tab4$percent.coverage,")",sep=""),sep=" ")
#tab4 <- tab4[,c(9,1,2,11,6)]
#names(tab4) <- c("Contrast","GO.ID","Term","Coverage","p")
#tab4 <- tab4[order(-tab4$Contrast),] # is this correct? 
tab4$p <- formatC(tab4$elim, format = "e", digits = 2) # tab4$p

#write.csv(tab3, "thoGODif.csv", row.names = F)
#write.csv(tab4, "abdGODif.csv", row.names = F)

# Visualize Lab DEG  ----
## Venn Diagram ----
e <- ggVennDiagram(list(
  Sham = stringent[which(stringent$comp == "control-sham" & stringent$tiss == "abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Immune = stringent[which(stringent$comp == "control-ecoli" & stringent$tiss == "abd" &  stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Heat = stringent[which(stringent$comp == "control-heat" & stringent$tiss == "abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Cold = stringent[which(stringent$comp == "control-cold" & stringent$tiss == "abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Starve = stringent[which(stringent$comp == "control-starve" & stringent$tiss == "abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
), label=c("count")) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +  theme(legend.position = "none")


f <- ggVennDiagram(list(
  Sham = stringent[which(stringent$comp == "control-sham" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Immune = stringent[which(stringent$comp == "control-ecoli" & stringent$tiss == "tho" &  stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Heat = stringent[which(stringent$comp == "control-heat" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Cold = stringent[which(stringent$comp == "control-cold" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Starve = stringent[which(stringent$comp == "control-starve" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
), label=c("count")) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +  theme(legend.position = "none")


(f + e) + 
  plot_annotation(tag_levels = 'A') 

#ggsave("vennSupp.tif", device = "tiff", width=5, height=8)

## Upset plot ----

upset(fromList(list(
  Sham = stringent[which(stringent$comp == "control-sham" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Immune = stringent[which(stringent$comp == "control-ecoli" & stringent$tiss == "tho" &  stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Heat = stringent[which(stringent$comp == "control-heat" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Cold = stringent[which(stringent$comp == "control-cold" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"], 
  Starve = stringent[which(stringent$comp == "control-starve" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
)), order.by="freq", text.scale = c(2, 2, 2, 2, 2, 2))


## PCA for significant DEG --------
stringentSig <- stringent[which(stringent$tiss=="all" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
keep_genes <- rowSums(counts(dds)) > 10
keepStringent <- ifelse(names(keep_genes) %in% stringentSig & keep_genes==TRUE, TRUE, FALSE)
ddsStringent <- dds[ keepStringent, ] 

stringentSigAbd <- stringent[which(stringent$tiss=="abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
keep_genesAbd <- rowSums(counts(ddsAbd)) > 10
keepstringentAbd <- ifelse(names(keep_genesAbd) %in% stringentSigAbd & keep_genesAbd==TRUE, TRUE, FALSE)
ddsstringentAbd <- ddsAbd[ keepstringentAbd, ]

stringentSigTho <- stringent[which(stringent$tiss=="tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"]
keep_genesTho <- rowSums(counts(ddsTho)) > 10
keepstringentTho <- ifelse(names(keep_genesTho) %in% stringentSigTho & keep_genesTho==TRUE, TRUE, FALSE)
ddsstringentTho <- ddsTho[ keepstringentTho, ]

nsub1 <- sum( rowMeans( counts(ddsStringent, normalized=TRUE)) > 5 )
nsub2 <- sum( rowMeans( counts(ddsstringentAbd, normalized=TRUE)) > 5 )
nsub3 <- sum( rowMeans( counts(ddsstringentTho, normalized=TRUE)) > 5 )

vsdatastringent <- vst(ddsStringent, blind=FALSE, nsub=nsub1)
vsdataAbdstringent <- vst(ddsstringentAbd, blind=FALSE, nsub=nsub2)
vsdataThostringent <- vst(ddsstringentTho, blind=FALSE, nsub=nsub3)


pcaDatastringent <- plotPCA(vsdatastringent, intgroup=c("trt", "tissue"), returnData=TRUE)
percentVarstringent <- round(100 * attr(pcaDatastringent, "percentVar"))
pcaDatastringent$tissue <- factor(pcaDatastringent$tissue, levels = c("thorax", "abdomen"))
pcaDatastringent<- pcaDatastringent %>%
  mutate(stress = factor(trt, levels= c("starve", "cold", "heat", "control", "sham", "ecoli"))) %>%
  mutate(stress = recode(stress, "control" = "Control","cold" = "Cold", "ecoli" = "Immune", 
                         "heat" = "Heat", "sham" = "Sham","starve" = "Starve"))
allPCAstringent <- ggplot(pcaDatastringent, aes(PC1, PC2, shape=tissue)) + # color=stress, 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarstringent[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarstringent[2],"% variance")) + 
  #scale_colour_colorblind() + 
  scale_linetype_discrete(labels = c("Thorax", "Abdomen")) + 
  guides(#color=guide_legend(title="Trt"),
    shape=guide_legend(title="Tissue"), 
    linetype=guide_legend(title="Tissue")) +
  scale_shape_manual(values=c(16, 1), labels = c("Thorax", "Abdomen")) +
  stat_ellipse(aes(linetype = tissue) ) + 
  theme_bw(base_size = 15) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))  

pcaDataAbdstringent <- plotPCA(vsdataAbdstringent, intgroup=c("trt", "colony"), returnData=TRUE)
percentVarAbdstringent <- round(100 * attr(pcaDataAbdstringent, "percentVar"))
pcaDataAbdstringent<- pcaDataAbdstringent %>%
  mutate(stress = factor(trt, levels= c("starve", "cold", "heat", "control", "sham", "ecoli"))) %>%
  mutate(stress = recode(stress, "control" = "Control","cold" = "Cold", "ecoli" = "Immune", 
                         "heat" = "Heat", "sham" = "Sham","starve" = "Starve"))
abdPCAstringent <- ggplot(pcaDataAbdstringent, aes(PC1, PC2, color=stress, shape=colony, group=trt)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarAbdstringent[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarAbdstringent[2],"% variance")) + 
  scale_colour_manual(values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7"), guide="none") + 
  guides(color=guide_legend(title="Trt"),
         shape=guide_legend(title="Colony", order = 1)) +
  stat_ellipse(linetype="dashed") + 
  scale_shape_manual(values=c(0,5,2))+
  theme_bw(base_size = 15) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

pcaDataThostringent <- plotPCA(vsdataThostringent, intgroup=c("trt", "colony"), returnData=TRUE)
percentVarThostringent <- round(100 * attr(pcaDataThostringent, "percentVar"))
pcaDataThostringent<- pcaDataThostringent %>%
  mutate(stress = factor(trt, levels= c("starve", "cold", "heat", "control", "sham", "ecoli"))) %>%
  mutate(stress = recode(stress, "control" = "Control","cold" = "Cold", "ecoli" = "Immune", 
                         "heat" = "Heat", "sham" = "Sham","starve" = "Starve"))
thoPCAstringent <- ggplot(pcaDataThostringent, aes(PC1, PC2, fill = stress, color=stress, shape=colony, group=trt)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarThostringent[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarThostringent[2],"% variance")) + 
  #scale_colour_colorblind() + 
  scale_fill_manual(values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7")) + 
  scale_color_manual(values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7")) + 
  guides(color=guide_legend(title="Trt", order = 1),
         shape=guide_legend(title="Colony", order = 2)) +
  stat_ellipse() + 
  scale_shape_manual(values=c(22, 23, 24))+
  theme_bw(base_size = 15) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"), legend.position="none") 

(allPCAstringent / (thoPCAstringent + abdPCAstringent)) + plot_annotation(tag_levels = 'A') + 
  plot_layout(guides='collect') + 
  plot_layout(widths = c(2, 1, 2))
#ggsave("PCAstringent2.tiff", device="tiff", width=10, height=8)

## DEG heat map ----
degList <- unique(unlist(list(
  Sham = stringent[which(stringent$comp == "control-sham" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 3),"gene"], 
  Immune = stringent[which(stringent$comp == "control-ecoli" & stringent$tiss == "tho" &  stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 3),"gene"], 
  Heat = stringent[which(stringent$comp == "control-heat" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 3),"gene"], 
  Cold = stringent[which(stringent$comp == "control-cold" & stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 3),"gene"], 
  Starve = stringent[which(stringent$comp == "control-starve" & stringent$tiss == "tho" & stringent$padj < 0.05& abs(stringent$log2FoldChange) > 3),"gene"]
)))

geneOrder <- arrange(stringent[which(stringent$gene %in% degList & stringent$tiss == "tho" & stringent$comp == "control-ecoli"), ], log2FoldChange)$gene

#rfGenes <- c("LOC100746301", "LOC100744278", "LOC105680509", "LOC100744835", "LOC100747034", "LOC100748925", "LOC100746622", "LOC100747233", "LOC100741422", "LOC100743553")

stringent %>%
  dplyr::filter(stringent$gene %in% degList & stringent$tiss == "tho" & stringent$comp != "control-ecoli") %>% # & stringent$comp != "sham-ecoli"
  mutate(comp = recode(comp, "control-ecoli" = "Immune", 
                       "control-sham" = "Sham", 
                       "control-heat" = "Heat", 
                       "control-cold" = "Cold", 
                       "control-starve" = "Starve")) %>%
  mutate(comp = factor(comp, levels = c("Immune", "Sham", "Heat","Cold", "Starve")), 
         gene = factor(gene, levels = geneOrder)) %>%
  ggplot(aes(x= comp, y=gene, fill= log2FoldChange)) +
  geom_tile() + 
  scale_fill_gradient2(high="blue", low="red", mid= "white", name= "Log Fold Change") + 
  xlab("Stressor") + 
  ylab("Gene") + 
  theme(axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black")) + 
  theme_bw() +
  # if want to add the important genes that show up as highly DEG (>3) # unique(data[which(data$imp ==1),"gene"])
  annotate("rect", xmin=c(0.5, 0.5, 0.5), xmax=c(5.5, 5.5, 5.5), 
           ymin=match(c("LOC100747233", "LOC105680509", "LOC100746301"), geneOrder) + 0.5,
           ymax=match(c("LOC100747233", "LOC105680509", "LOC100746301"), geneOrder) - 0.5,
           colour="black", fill="transparent", size=0.5)

#ggsave("degInteract.tif", height = 15, width=5)


# Format Field Data ----

rownames(sampleField) <- sampleField[,"sample_name"]
rownames(countsField) <- countsField[,1]

countsField <- countsField[, -which(names(countsField) %in% c("...1"))]

ddsField <- DESeqDataSetFromMatrix(countData=countsField, 
                                   colData=sampleField, 
                                   design=~site + sampleDate, tidy = F) # colony + tissue + round +  trt

ddsSite<- DESeqDataSetFromMatrix(countData=countsField, 
                                 colData=sampleField, 
                                 design=~sampleDate + site, tidy = F) # colony + tissue + round +  trt

keep_genes <- rowSums(counts(ddsField)) > 10
ddsField <- ddsField[ keep_genes, ]

keep_genes <- rowSums(counts(ddsSite)) > 10
ddsSite <- ddsSite[ keep_genes, ]

# Run Field DESeq ----

ddsField <- DESeq(ddsField)
ddsSite <- DESeq(ddsSite)


# Summarise Field Contrasts ----
resSite <- as.data.frame(results(ddsSite, contrast = c("site", "Ridge", "Valley"), lfcThreshold = 0.5, alpha =0.05, altHypothesis="greaterAbs"))
dim(resSite[which(abs(resSite$log2FoldChange) > 0.5 & resSite$padj < 0.05),]) #0

resDate1 <- as.data.frame(results(ddsField, contrast = c("sampleDate", "7/28/23", "7/26/23"), lfcThreshold = 0.5, alpha =0.05, altHypothesis="greaterAbs"))
dim(resDate1[which(abs(resDate1$log2FoldChange) > 0.5 & resDate1$padj < 0.05),]) #43

resDate2 <- as.data.frame(results(ddsField, contrast = c("sampleDate", "7/28/23", "7/29/23"), lfcThreshold = 0.5, alpha =0.05, altHypothesis="greaterAbs"))
dim(resDate2[which(abs(resDate2$log2FoldChange) > 0.5 & resDate2$padj < 0.05),]) #14

resDate3 <- as.data.frame(results(ddsField, contrast = c("sampleDate", "7/28/23", "6/21/23"), lfcThreshold = 0.5, alpha =0.05, altHypothesis="greaterAbs"))
dim(resDate3[which(abs(resDate3$log2FoldChange) > 0.5 & resDate3$padj < 0.05),]) #2

## Supplemental check -- after heatwave vs. before 
resDate4 <- as.data.frame(results(ddsField, contrast = c("sampleDate", "7/29/23", "7/26/23"), lfcThreshold = 0.5, alpha =0.05, altHypothesis="greaterAbs"))
dim(resDate4[which(abs(resDate4$log2FoldChange) > 0.5 & resDate4$padj < 0.05),]) #0

resSite$gene <- rownames(resSite)
resDate1$gene <- rownames(resDate1)
resDate2$gene <- rownames(resDate2)
resDate3$gene <- rownames(resDate3)

deseq23 <- resSite[which(abs(resSite$log2FoldChange) > 0.5 & resSite$padj < 0.05),] %>%
  mutate(comp = "site") %>%
  full_join(resDate1[which(abs(resDate1$log2FoldChange) > 0.5 & resDate1$padj < 0.05),]) %>%
  mutate(comp = ifelse(is.na(comp), "pre_26", comp)) %>%
  full_join(resDate2[which(abs(resDate2$log2FoldChange) > 0.5 & resDate2$padj < 0.05),]) %>%
  mutate(comp = ifelse(is.na(comp), "post_29", comp)) %>%
  full_join(resDate3[which(abs(resDate3$log2FoldChange) > 0.5 & resDate3$padj < 0.05),]) %>%
  mutate(comp = ifelse(is.na(comp), "pre_6_21", comp)) 

#write.csv(deseq23, "DEseqResArb23.csv")

# Random Forest ----
## Identify important genes ----

impGenesFunc <- function(nreps = nreps, nTopGenes = nTopGenes, countDat = countDat, sampleDat = sampleDat) {
  
  countDatT <- dplyr::mutate_all(as.data.frame(t(countDat)), function(x) as.numeric(x))
  
  # number of times to run
  #nreps <- nreps
  #nTopGenes <- nTopGenes
  impGenes <- vector(mode='list', length=nreps)
  seeds <- vector(mode="numeric", length=nreps)
  
  ### start the random forest 
  
  for (thisrep in seq_len(nreps)){
    cat("rep", thisrep, "of", nreps)
    
    curT <- as.numeric(Sys.time())
    set.seed(curT)
    
    # Random forest model 
    thismodel <- randomForest(x=countDatT, y=as.factor(sampleDat$trt), 
                              num.trees=1000)
    
    # Save 20 most important genes 
    imp.tempAll <- abs(thismodel$importance[,])
    tAll <- order(imp.tempAll,decreasing=TRUE)
    gn.impAll <- names(imp.tempAll)[tAll]
    gn.10All <- gn.impAll[1:nTopGenes] 
    impGenes[[thisrep]] <- gn.10All
    
    seeds[[thisrep]] <- curT
    
  }
  
  
  impGeneSumAll <- dplyr::arrange(as.data.frame(table(unlist(impGenes))), -Freq) %>%
    mutate(prop = Freq/(nreps))
  seedsAll <- seeds
  
  return(impGeneSumAll)
  
}

impGenesCVFunc <- function(nreps = nreps, nfold = nfold, nTopGenes = nTopGenes, 
                           countDat = countDat, sampleDat = sampleDat) {
  
  countDatT <- dplyr::mutate_all(as.data.frame(t(countDat)), function(x) as.numeric(x))
  
  index.col1 <- which(sampleDat$colony=="1")
  index.col2 <- which(sampleDat$colony=="2")
  index.col3 <- which(sampleDat$colony=="3")
  
  impGenes <- vector(mode='list', length=nreps)
  
  impGenesSum <- vector(mode='list', length=nreps)
  
  
  ### start the random forest 
  
  for (thisrep in seq_len(nreps)){
    cat("REP", thisrep, "OF", nreps)
    
    curT <- as.numeric(Sys.time())
    set.seed(curT)
    
    # create a 10-way split:
    kvals.col1 <- cbind(index.col1, sample(seq_len(nfold), size=length(index.col1), replace=TRUE))
    kvals.col2 <- cbind(index.col2, sample(seq_len(nfold), size=length(index.col2), replace=TRUE))
    kvals.col3 <- cbind(index.col3, sample(seq_len(nfold), size=length(index.col3), replace=TRUE))
    
    kvals <- as.data.frame(rbind(kvals.col1,kvals.col2,kvals.col3))
    kvals <- kvals[order(kvals[,1]),][,2]
    
    impGenes <- vector(mode='list', length=nfold) 
    
    for(thisk in sort(unique(kvals))) {
      cat("fold", thisk, "of", nfold)
      trainsetSample <- sampleDat[kvals != thisk, ]
      testsetSample  <- sampleDat[kvals == thisk, ]
      
      trainsetCount <- countDatT[kvals != thisk, ]
      testsetCount  <- countDatT[kvals == thisk, ]
      
      # run training set on the full set of genes
      thismodel <- randomForest(x=trainsetCount, y=as.factor(trainsetSample$trt), 
                                num.trees=1000)
      
      
      imp.tempAll <- abs(thismodel$importance[,])
      tAll <- order(imp.tempAll,decreasing=TRUE)
      gn.impAll <- names(imp.tempAll)[tAll]
      gn.10All <- gn.impAll[1:nTopGenes] 
      impGenes[[thisk]] <- gn.10All
      
    }
    
    impGenesSum[[thisrep]] <- as.data.frame(sort(table(unlist(impGenes))))
    
  }
  
  impGeneSumAll <- data.frame(do.call(rbind, impGenesSum)) %>%
    group_by(Var1) %>%
    summarise(Freq = sum(Freq)) %>%
    arrange(-Freq) %>%
    ungroup() %>%
    mutate(prop = Freq / (nfold*nreps*nTopGenes))
  
  
  return(impGeneSumAll)
  
}

## Apply function ----

impAll <- impGenesCVFunc(nreps = 100, nTopGenes = 20, nfold=12,
                         countDat = countsLab, sampleDat = sampleLab)
impTho <- impGenesCVFunc(nreps = 100, nTopGenes = 10, nfold=12,
                         countDat = countTho, sampleDat = sampleLabTho)
impAbd <- impGenesCVFunc(nreps = 100, nTopGenes = 10, nfold=12,
                         countDat = countAbd, sampleDat = sampleLabAbd)

# Load important genes here directly from supplement so don't have to run entire thing above
# Also, because random this reproduces the results in the manuscript. 
important <- as.data.frame(read_excel("supplementbbTranscriptome23Oct.xlsx", sheet = "TableS8", skip = 5))

impAll <- important[which(important$tiss == "all"),]
impTho <- important[which(important$tiss == "tho"),]
impAbd <- important[which(important$tiss == "abd"),]


# Predict Field Stress & OOB Error ----
# This predicts the OOB error rates for each tissue and overall 
# for each of the 3 models (a full model using all the genes, 
  # a submodel using just the genes identified as important from the most recent RF run (model 1), 
  # and a top model based on the most important genes across 100 runs (identified in last section))


# We simultaneously predict on the field data, so the OOB predictions are accurate to the actual trees in each model

predictRF <- function(nreps = nreps, impGenes = impGenes, nTopGenes = nTopGenes, degs = degs,
                      countsDatLab = countsDatLab, countsDatField = countsDatField,
                      sampleDatLab = sampleDatLab) {
  
  impGenes <- as.data.frame(impGenes)
  topGenes <- impGenes[1:nTopGenes, 1] # will need to ensure the naming (not currently unique) is updated so you don't overwrite 
  
  countAll <- dplyr::mutate_all(as.data.frame(t(countsDatLab)), function(x) as.numeric(x))
  countsField <- dplyr::mutate_all(as.data.frame(t(countsDatField)), function(x) as.numeric(x))
  
  # number of times to run
  
  thisDEGpred <- vector(mode='list', length=nreps)
  thisDEGconf <- vector(mode='list', length=nreps)
  
  thisToppred <- vector(mode='list', length=nreps)
  thisTopconf <- vector(mode='list', length=nreps)
  
  thisSubpred <- vector(mode='list', length=nreps)
  thisSubconf <- vector(mode='list', length=nreps)
  
  thispred <- vector(mode='list', length=nreps)
  thisconf <- vector(mode='list', length=nreps)
  
  
  ### start the random forest 
  
  for (thisrep in seq_len(nreps)){
    cat("rep", thisrep, "of", nreps)
    
    curT <- as.numeric(Sys.time())
    set.seed(curT)
    
    # run training set on the full set of genes
    thismodel <- randomForest(x=countAll, y=as.factor(sampleDatLab$trt), 
                              num.trees=1000)
    thispred[[thisrep]]  <- predict(thismodel, countsField, type="prob")
    thisconf[[thisrep]]<- thismodel[[5]]
    
    # take 10 most important genes 
    imp.tempAll <- abs(thismodel$importance[,])
    tAll <- order(imp.tempAll,decreasing=TRUE)
    gn.impAll <- names(imp.tempAll)[tAll]
    gn.10All <- gn.impAll[1:nTopGenes] 
    tAll <- is.element(colnames(countAll),gn.10All)
    sig.esetAll <- countAll[,c(tAll)]
    
    # Subset of the full model 
    thisSubmodel <- randomForest(x=sig.esetAll, y=as.factor(sampleDatLab$trt) ,
                                 num.trees=1000)
    thisSubpred[[thisrep]]  <- predict(thisSubmodel, countsField, type="prob")
    thisSubconf[[thisrep]] <- thisSubmodel[[5]]
    
    # Top 10 most important genes from model assessment 
    
    thisTopmodel <- randomForest(x=countAll[,topGenes], y=as.factor(sampleDatLab$trt) ,
                                 num.trees=1000)
    thisToppred[[thisrep]]  <- predict(thisTopmodel, countsField, type="prob")
    thisTopconf[[thisrep]] <- thisTopmodel[[5]]
    
    # DEG's 
    
    thisDEGmodel <- randomForest(x=countAll[,unique(degs)], y=as.factor(sampleDatLab$trt) ,
                                 num.trees=1000)
    thisDEGpred[[thisrep]]  <- predict(thisDEGmodel, countsField, type="prob")
    thisDEGconf[[thisrep]] <- thisDEGmodel[[5]]

    
  }
  
  # make arrays 
  arrDEGpred <- array( unlist(thisDEGpred) , c(ncol(countsDatField),6,nreps) )
  arrToppred <- array( unlist(thisToppred) , c(ncol(countsDatField),6,nreps) )
  arrthisSubpred <- array( unlist(thisSubpred) , c(ncol(countsDatField),6,nreps) )
  arrthispred <- array( unlist(thispred) , c(ncol(countsDatField),6,nreps) )
  
  arrthisDEGconf <- array( unlist(thisDEGconf) , c(6,7,nreps) )
  arrthisTopconf <- array( unlist(thisTopconf) , c(6,7,nreps) )
  arrthisSubconf <- array( unlist(thisSubconf) , c(6,7,nreps) )
  arrthisconf <- array( unlist(thisconf) , c(6,7,nreps) )
  
  # clacluate mean and 95% CI 
  thisDEGpredMean <- apply(arrDEGpred, 1:2, mean)
  thisDEGpredLower <- apply(arrDEGpred, 1:2,  quantile, prob = 0.025)
  thisDEGpredUpper <- apply(arrDEGpred, 1:2,  quantile, prob = 0.975)
  
  thisDEGpredSum <- as.data.frame(cbind(thisDEGpredMean, thisDEGpredLower, thisDEGpredUpper))
  colnames(thisDEGpredSum) <- c(paste(colnames(thisDEGpred[[1]]), "mean", sep="_"), 
                                paste(colnames(thisDEGpred[[1]]), "lower", sep="_"), 
                                paste(colnames(thisDEGpred[[1]]), "upper", sep="_"))
  thisDEGpredSum$sample <- rownames(thisDEGpred[[1]])
  
  thisToppredMean <- apply(arrToppred, 1:2, mean)
  thisToppredLower <- apply(arrToppred, 1:2,  quantile, prob = 0.025)
  thisToppredUpper <- apply(arrToppred, 1:2,  quantile, prob = 0.975)
  
  thisToppredSum <- as.data.frame(cbind(thisToppredMean, thisToppredLower, thisToppredUpper))
  colnames(thisToppredSum) <- c(paste(colnames(thisToppred[[1]]), "mean", sep="_"), 
                                paste(colnames(thisToppred[[1]]), "lower", sep="_"), 
                                paste(colnames(thisToppred[[1]]), "upper", sep="_"))
  thisToppredSum$sample <- rownames(thisToppred[[1]])
  
  thisSubpredMean <- apply(arrthisSubpred, 1:2, mean)
  thisSubpredLower <- apply(arrthisSubpred, 1:2,  quantile, prob = 0.025)
  thisSubpredUpper <- apply(arrthisSubpred, 1:2,  quantile, prob = 0.975)
  
  thisSubpredSum <- as.data.frame(cbind(thisSubpredMean, thisSubpredLower, thisSubpredUpper))
  colnames(thisSubpredSum) <- c(paste(colnames(thisSubpred[[1]]), "mean", sep="_"), 
                                paste(colnames(thisSubpred[[1]]), "lower", sep="_"), 
                                paste(colnames(thisSubpred[[1]]), "upper", sep="_"))
  thisSubpredSum$sample <- rownames(thisSubpred[[1]])
  
  thispredMean <- apply(arrthispred, 1:2, mean)
  thispredLower <- apply(arrthispred, 1:2,  quantile, prob = 0.025)
  thispredUpper <- apply(arrthispred, 1:2,  quantile, prob = 0.975)
  
  thispredSum <- as.data.frame(cbind(thispredMean, thispredLower, thispredUpper))
  colnames(thispredSum) <- c(paste(colnames(thispred[[1]]), "mean", sep="_"), 
                             paste(colnames(thispred[[1]]), "lower", sep="_"), 
                             paste(colnames(thispred[[1]]), "upper", sep="_"))
  thispredSum$sample <- rownames(thispred[[1]])
  
  thisDEGconfMean <- apply(arrthisDEGconf, 1:2, mean)
  thisDEGconfLower <- apply(arrthisDEGconf, 1:2,  quantile, prob = 0.025)
  thisDEGconfUpper <- apply(arrthisDEGconf, 1:2,  quantile, prob = 0.975)
  
  thisDEGconfSum <- as.data.frame(cbind(thisDEGconfMean, thisDEGconfLower, thisDEGconfUpper))
  colnames(thisDEGconfSum) <- c(paste(colnames(thisDEGconf[[1]]), "mean", sep="_"), 
                                paste(colnames(thisDEGconf[[1]]), "lower", sep="_"), 
                                paste(colnames(thisDEGconf[[1]]), "upper", sep="_"))
  thisDEGconfSum$sample <- rownames(thisDEGconf[[1]])
  
  thisTopconfMean <- apply(arrthisTopconf, 1:2, mean)
  thisTopconfLower <- apply(arrthisTopconf, 1:2,  quantile, prob = 0.025)
  thisTopconfUpper <- apply(arrthisTopconf, 1:2,  quantile, prob = 0.975)
  
  thisTopconfSum <- as.data.frame(cbind(thisTopconfMean, thisTopconfLower, thisTopconfUpper))
  colnames(thisTopconfSum) <- c(paste(colnames(thisTopconf[[1]]), "mean", sep="_"), 
                                paste(colnames(thisTopconf[[1]]), "lower", sep="_"), 
                                paste(colnames(thisTopconf[[1]]), "upper", sep="_"))
  thisTopconfSum$sample <- rownames(thisTopconf[[1]])
  
  thisSubconfMean <- apply(arrthisSubconf, 1:2, mean)
  thisSubconfLower <- apply(arrthisSubconf, 1:2,  quantile, prob = 0.025)
  thisSubconfUpper <- apply(arrthisSubconf, 1:2,  quantile, prob = 0.975)
  
  thisSubconfSum <- as.data.frame(cbind(thisSubconfMean, thisSubconfLower, thisSubconfUpper))
  colnames(thisSubconfSum) <- c(paste(colnames(thisSubconf[[1]]), "mean", sep="_"), 
                                paste(colnames(thisSubconf[[1]]), "lower", sep="_"), 
                                paste(colnames(thisSubconf[[1]]), "upper", sep="_"))
  thisSubconfSum$sample <- rownames(thisSubconf[[1]])
  
  thisconfMean <- apply(arrthisconf, 1:2, mean)
  thisconfLower <- apply(arrthisconf, 1:2,  quantile, prob = 0.025)
  thisconfUpper <- apply(arrthisconf, 1:2,  quantile, prob = 0.975)
  
  thisconfSum <- as.data.frame(cbind(thisconfMean, thisconfLower, thisconfUpper))
  colnames(thisconfSum) <- c(paste(colnames(thisconf[[1]]), "mean", sep="_"), 
                             paste(colnames(thisconf[[1]]), "lower", sep="_"), 
                             paste(colnames(thisconf[[1]]), "upper", sep="_"))
  thisconfSum$sample <- rownames(thisconf[[1]])

  return(list(thisToppredSum, thisSubpredSum, thispredSum, thisDEGpredSum,
              thisTopconfSum, thisSubconfSum, thisconfSum, thisDEGconfSum))
  
}

## Apply function ----
# a little messy because added at revisions and didnt want to change the function and parameters too much
stringentAll <- unique(stringent[which(stringent$tiss == "all" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"])
stringentTho <- unique(stringent[which(stringent$tiss == "tho" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"])
stringentAbd <- unique(stringent[which(stringent$tiss == "abd" & stringent$padj < 0.05 & abs(stringent$log2FoldChange) > 0.5),"gene"])

# Actually apply function here 
predAll <- predictRF(nreps = 100, impGenes = impAll, nTopGenes = 20, degs = stringentAll,
          countsDatLab = countsLab, countsDatField = countsField,
          sampleDatLab = sampleLab)

predTho <- predictRF(nreps = 100, impGenes = impTho, nTopGenes = 10,
                     countsDatLab = countTho, countsDatField = countsField, degs = stringentTho,
                     sampleDatLab = sampleLabTho)

predAbd <- predictRF(nreps = 100, impGenes = impAbd, nTopGenes = 10,
                     countsDatLab = countAbd, countsDatField = countsField,degs = stringentAbd,
                     sampleDatLab = sampleLabAbd)

# Figures Random Forest ----
# read in data if not running full script at once 

## Confusion Matrix ----
errorFigure <- function (x) {
  title <- deparse(substitute(x))
  OOB <- round(mean(x$class.error_mean), 2)
  corCorrect <- max(rowSums(x[,1:6]))
  df <- x %>% 
  dplyr::select(control_mean:starve_mean, sample) %>%
  pivot_longer(cols = 1:6, names_to = "var", values_to = "corr") %>%
  mutate(corr = corr/ corCorrect) %>%
  separate(var, c("var", NA)) 
  
  ggplot(df, aes(x=sample, y=var, fill=corr)) + 
    geom_tile() + 
    xlab("Actual Stressor") + 
    ylab("Predicted Stressor") +  
    scale_x_discrete(labels=c("cold" = "Cold", "control" = "Control", "ecoli" = "Immune", 
                              "heat" = "Heat", "sham" = "Sham", "starve" = "Starve")) + 
    scale_y_discrete(labels=c("cold" = "Cold", "control" = "Control", "ecoli" = "Immune", 
                              "heat" = "Heat", "sham" = "Sham", "starve" = "Starve")) + 
    ggtitle(paste("OOB Error: ", OOB)) + 
    scale_fill_gradient(limits=c(0,1), low = "white", high = "red", name = "Mean Prop. \nClassified") +
    theme_bw(base_size = 15)+ 
    theme(axis.text=element_text(color="black")) 
  
}

errorFigure(predTho[[5]]) +  errorFigure(predAbd[[5]]) + #errorFigure(predAll[[5]]) + 
  plot_annotation(tag_levels = 'A') +   plot_layout(guides='collect') 

ggsave("confusedTop.tif", device = "tiff", width=17, height=5)

(errorFigure(predTho[[5]]) |  errorFigure(predTho[[6]]) | errorFigure(predTho[[8]]) | errorFigure(predTho[[7]]))  / 
  (errorFigure(predAbd[[5]]) | errorFigure(predAbd[[6]]) | errorFigure(predAbd[[8]]) | errorFigure(predAbd[[7]])) + # /
  #(errorFigure(predAll[[5]]) | errorFigure(predAll[[6]]) | errorFigure(predAll[[8]]) | errorFigure(predAll[[7]])) +
  plot_annotation(tag_levels = 'A') +   plot_layout(guides='collect') 

#ggsave("confusedSup.tif", device = "tiff", width=24, height=15)
#ggsave("confusedSup2.tif", device = "tiff", width=24, height=10)


## Field Predictions ----
predSummary <- function (x) {
  mod <- deparse(substitute(x))
  df <- x %>%
    pivot_longer(cols= 1:18, names_to="variable", values_to = "value") %>%
    separate(variable, sep="_", c("stress", "var")) %>% 
    pivot_wider(names_from = var, values_from=value) %>%
    mutate(tissueMod = mod)
  return(df)
}

# Don't actually need to run the whole thing beacuse we are only using pred from tho top
allPred <- predSummary(predTho[[4]]) %>% # DEG
  full_join(predSummary(predTho[[3]])) %>%  # all
  full_join(predSummary(predTho[[2]])) %>% # sub
  full_join(predSummary(predTho[[1]])) %>% # top
  full_join(predSummary(predAbd[[4]])) %>%
  full_join(predSummary(predAbd[[3]])) %>%
  full_join(predSummary(predAbd[[2]])) %>%
  full_join(predSummary(predAbd[[1]])) %>%
  full_join(predSummary(predAll[[4]])) %>%
  full_join(predSummary(predAll[[3]])) %>%
  full_join(predSummary(predAll[[2]])) %>%
  full_join(predSummary(predAll[[1]]))

### Prediction Figure ----
allPred2 <- allPred %>%
  filter(tissueMod == "predTho[[1]]") %>% # the model you want
  rename(extractID = sample) %>%
  full_join(samplesField) 

allPred3 <- allPred2 %>%
  group_by(extractID, sampleDate, site) %>%
  summarise(n=n()) %>%
  group_by(sampleDate, site) %>%
  mutate(beeID= 1:4) %>% 
  ungroup() %>%
  full_join(allPred2) %>%
  mutate(stress = recode(stress, "control" = "Control", "cold" = "Cold", "ecoli" = "Immune", 
                         "heat" = "Heat", "sham" = "Sham","starve" = "Starve"), 
         sampleDate = recode(sampleDate, "7/28/23" = "Heatwave\n (PM)", 
                             "6/21/23" = "Early Summer\n (AM)", 
                             "7/29/23" = "Post-Heatwave\n (AM)",
                             "7/26/23" = "Pre-Heatwave\n (AM)")) %>%
  mutate(sampleDate=fct_relevel(as.factor(sampleDate),c("Early Summer\n (AM)", "Pre-Heatwave\n (AM)",
                                                        "Heatwave\n (PM)", "Post-Heatwave\n (AM)"))) %>%
  mutate(stress = fct_relevel(as.factor(stress), c("Starve", "Cold", "Heat", "Control", "Sham", "Immune")))

f2c <- function(x) { (x-32) * (5/9) }

allPred4 <- allPred3 %>%
  full_join(samplesField[,c("extractID", "temperature_F")]) %>%
  mutate(tempC = f2c(temperature_F)) 

ggplot(allPred4[which(allPred4$mean > 0.15),], aes(x = factor(beeID), y= mean, color = factor(stress))) + 
  geom_point() +
  geom_point(aes(y = tempC/40), color = "#999999", pch= 8) +
  scale_y_continuous("Predicted Probability", sec.axis = sec_axis(~.*40, name = "Collection Temperature (°C)")) +
  xlab("Bee ID") + 
  #ylab("Predicted Probability") +  
  scale_color_manual(values= c("#000000" ,"#56B4E9","#E69F00" , "#0072B2", "#CC79A7"), name = "Stressor") + 
  geom_errorbar(aes(ymin=lower, ymax= upper)) + 
  facet_nested( ~ sampleDate + site , ) + 
  theme_bw() + 
  theme(axis.text=element_text(color="black")) 
ggsave("predictFieldSub.tif", device = "tiff", width=6, height=4) #height=4 for old figure

ggplot(allPred4, aes(x = factor(beeID), y= mean, color = factor(stress))) + 
  geom_point() +
  geom_point(aes(y = tempC/40), color = "#999999", pch= 8) +
  scale_y_continuous("Predicted Probability", sec.axis = sec_axis(~.*40, name = "Collection Temperature (°C)")) +
  xlab("Bee ID") + 
  #ylab("Predicted Probability") +  
  scale_color_manual(values= c("#000000" ,"#56B4E9","#E69F00" , "#009E73",  "#0072B2", "#CC79A7"), name = "Stressor") + 
  geom_errorbar(aes(ymin=lower, ymax= upper)) + 
  facet_nested( ~ sampleDate + site , ) + 
  theme_bw() + 
  theme(axis.text=element_text(color="black")) 
ggsave("predictFieldDone.tif", device = "tiff", width=8, height=6) #height=4 for old figure


# Compare Environment ----

#comp <- allPred2  %>% #allTop
  #filter(tissueMod == "thoTop") %>%
  #dplyr::rename(extractID = sample) %>%
  #dplyr::select(ends_with("_mean"), extractID) %>%
  #full_join(samplesField) %>% 
  #mutate(temp_c = f2c(temperature_F)) %>%
  #dplyr::select(extractID, sampleDate, time.of.day, site, temp_c, # temperature_F, 
   #             rh, heat_index, dew_point) # ,ends_with("_mean")


# site differences
# THIS IS THE UPDATED MODEL FOR REVISIONS
anova(lm(tempC ~ site * sampleDate, (allPred4[,c("tempC", "site", "sampleDate")])))

## site dif. figure 
ggplot(allPred4, aes(x= sampleDate, y=tempC)) + 
  geom_boxplot(aes(color= time.of.day, fill=site)) + 
  xlab("Sample Round") + 
  ylab("Temperature (°C)") + 
  scale_color_manual(labels = c('Morning','Evening'), name="Time of Day", values= c("grey1", "grey50")) + 
  scale_fill_manual(labels = c('Mountain (Forest)','Valley (Arboretum)'), name="Site", values= c("#B8E186", "#92C5DE")) + 
  scale_x_discrete(labels= c("Early Summer\n (AM)", "Pre-Heatwave\n (AM)",
                             "Heatwave\n (PM)", "Post-Heatwave\n (AM)")) +
  guides(fill = guide_legend(order = 1), 
         color = guide_legend(order = 2)) + 
  theme_bw(base_size = 15) + 
  theme(axis.text=element_text(color="black")) 
#ggsave("supSite.tif", device = "tiff", width=5, height=7)


# Heat
summary(lm(mean ~ site, allPred4[which(allPred4$stress == "Heat"),]))
summary(lm(mean ~ sampleDate, allPred4[which(allPred4$stress == "Heat"),]))
TukeyHSD(aov(mean ~ sampleDate, allPred4[which(allPred4$stress == "Heat"),]))
summary(lm(mean ~ time.of.day, allPred4[which(allPred4$stress == "Heat"),])) 
summary(lm(mean ~ tempC, allPred4[which(allPred4$stress == "Heat"),]))

# Cold
summary(lm(mean ~ site, allPred4[which(allPred4$stress == "Cold"),]))
summary(lm(mean ~ sampleDate, allPred4[which(allPred4$stress == "Cold"),]))
TukeyHSD(aov(mean ~ sampleDate, allPred4[which(allPred4$stress == "Cold"),]))
summary(lm(mean ~ time.of.day, allPred4[which(allPred4$stress == "Cold"),])) 
summary(lm(mean ~ tempC, allPred4[which(allPred4$stress == "Cold"),]))

# Starve
summary(lm(mean ~ site, allPred4[which(allPred4$stress == "Starve"),]))
summary(lm(mean ~ sampleDate, allPred4[which(allPred4$stress == "Starve"),]))
TukeyHSD(aov(mean ~ sampleDate, allPred4[which(allPred4$stress == "Starve"),]))
summary(lm(mean ~ time.of.day, allPred4[which(allPred4$stress == "Starve"),])) 


## Visualize ----

A <- ggplot(allPred2[which(allPred2$stress == "heat"),], aes(x= f2c(temperature_F), y=mean)) + 
  geom_point() + #aes(pch= site)
  annotate("text", x = 32, y = 0.5, label = "paste(R ^ 2, \" = 0.13, p = 0.02\")", parse = TRUE, size = 15/.pt) + 
  annotate("text", x = 32, y = 0.47, label = paste("F[\"1,30\"] ==", 5.78),  parse = TRUE, size = 15/.pt) + 
  geom_errorbar(aes(ymin=lower, ymax= upper)) + 
  xlab("Temperature (°C)") + 
  ylab("Probability Heat Stress") + 
  theme_bw(base_size = 15) + 
  geom_smooth(method="lm", formula=y ~ x,color="black") +
  theme(axis.text=element_text(color="black"))

B <- ggplot(allPred2[which(allPred2$stress == "heat"),], aes(x= sampleDate, y=mean)) + 
  annotate("text", x = 3.5, y = 0.5, label = "paste(R ^ 2, \" = 0.09, p = 0.13\")", parse = TRUE, size = 15/.pt) + 
  annotate("text", x = 3.5, y = 0.47, label = paste("F[\"3,28\"] ==", 2.04),  parse = TRUE, size = 15/.pt) + 
  geom_boxplot() + #aes(color= site)
  xlab("Sample Round") + 
  ylab("Mean Probability Heat Stress") + 
  theme_bw(base_size = 15) + 
  scale_x_discrete(labels= c("Early Summer\n (AM)", "Pre-Heatwave\n (AM)",
                             "Heatwave\n (PM)", "Post-Heatwave\n (AM)")) + #, guide = guide_axis(n.dodge = 2)
  geom_smooth(method="lm", formula=y ~ x,color="black") +
  theme(axis.text=element_text(color="black"))

C <- ggplot(allPred2[which(allPred2$stress == "cold"),], aes(x= f2c(temperature_F), y=mean)) + 
  geom_point() + 
  annotate("text", x = 32, y = 0.5, label = "paste(R ^ 2, \" = 0.29, p < 0.01\")", parse = TRUE, size = 15/.pt) + 
  annotate("text", x = 32, y = 0.47, label = paste("F[\"1,30\"] ==", 13.72),  parse = TRUE, size = 15/.pt) + 
  geom_errorbar(aes(ymin=lower, ymax= upper)) + 
  xlab("Temperature (°C)") + 
  ylab("Probability Cold Stress") + 
  theme_bw(base_size = 15) + 
  geom_smooth(method="lm", formula=y ~ x,color="black") +
  theme(axis.text=element_text(color="black"))


ggpubr::ggarrange(A, C, B, nrow = 1, ncol=3, labels = "AUTO")
ggsave("corField.tif", device = "tiff", width=16, height=7)

# Compare RF and DEG ----

impGenes <- function (x, nTopGenes) {
  tissue <- tolower(strsplit(deparse(substitute(x)), "imp")[[1]][2]) 
  df <- x[1:nTopGenes, ] %>%
    dplyr::rename(gene = Var1) %>% 
    mutate(tiss = tissue) %>%
    left_join(stringent) %>%
    filter(comp != "control-ecoli") %>%#, 
    # padj <0.01) %>%
    separate(comp, c("control", "comp"), "-") %>%
    mutate(log2FoldChange = round(log2FoldChange, 2)) %>%
    dplyr::select(tiss, gene, comp, Freq, log2FoldChange) %>%
    pivot_wider(names_from = comp, values_from = log2FoldChange)
  return(df)
}

# Table 4
tab4 <- impGenes(impAll, 20) %>%
  full_join(impGenes(impTho, 10)) %>%
  full_join(impGenes(impAbd, 10))

tab4 <- as.data.frame(tab4) #%>%
  #rename("Gene" = "gene") %>%
  #full_join(names)
  
#write.csv(tab4, "tab4.csv")

ggVennDiagram(list(
  #Pooled = tab4[which(tab4$tiss == "all"),"gene"], 
  Thorax = tab4[which(tab4$tiss == "tho"),"gene"], 
  Abdomen = tab4[which(tab4$tiss == "abd"),"gene"]),
  label=c("count")) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

