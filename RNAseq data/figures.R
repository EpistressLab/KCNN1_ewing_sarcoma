############ Packages ############

require(DESeq2)
library(biomaRt)
library(genefilter,quietly=TRUE)
library(RColorBrewer)
library(gplots)
library(fdrtool)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(pvclust)
library(circlize)
library(pathview)
library(DOSE)
library(fgsea)
library(IHW)
library(org.Hs.eg.db)
library(tximport)
library(tidyverse)

############ Arguments & data ############


INFILE <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/RNAseq_ewing/deseq2_kallisto/kallisto_conditions.tab"
OUTDIR <- "/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/projects/RNAseq_ewing/deseq2_kallisto"
ASSEMBLY <- "hg38"
LOGFCTHRESHOLD <- 1
FDRTHRESHOLD <- 0.05

##################################
# 			with kallisto 		#
##################################

sampleTable <- read.table(INFILE, header=TRUE,sep = '\t')

conditions<-factor( sampleTable[ , 3] )
uniq_conds <- unique(conditions)

files <- sampleTable[,2]
names(files) <- sampleTable[,1]
txi.kallisto <- tximport(files, type="kallisto", txOut=T)
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
dds <- DESeq(dds)

# Normalized counts matrix
matrix <- counts(dds,normalized=TRUE)
# Transforming raw counts
rld <- rlogTransformation(dds, blind=TRUE)
mrld <- assay(rld)

########## Boxplot genes ########

fpkm_exp <- fpkm(dds)
gene_selected <- "KCNN1"
count_norm <- fpkm_exp
transcripts_selected <- c("XM_011528004.2","NR_170374.1","NR_170373.1","NM_001386974.1","NM_001386975.1","NM_001386977.1","NM_001386976.1","NM_002248.5")
data_boxplot <- t(count_norm[which(rownames(count_norm) %in% transcripts_selected),])
data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

for (transcript in colnames(data_boxplot)){
  sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
  data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
}
data_boxplot <- data_boxplot_2
transcripts_nicknames <- data.frame('transcript'=c("XM_011528004.2","NR_170374.1","NR_170373.1","NM_001386974.1","NM_001386975.1","NM_001386977.1","NM_001386976.1","NM_002248.5"),'nickname'=c('H','G','F','E','D','C','B','A'))
data_boxplot <- merge(data_boxplot,transcripts_nicknames)

outPlot = c(paste(OUTDIR,"/KCNN1_expression_with_ASP14.pdf",sep=""))
pdf(outPlot, width=9.9, height=8)
ggplot(data_boxplot, aes(x = nickname, y = expression, fill=transcript)) +
  geom_boxplot(alpha=0.5) +
  theme_minimal() +
  guides(fill = "none", size="none") +
  xlab("KCNN1 transcripts") +
  ylab("Normalized expression of KCNN1 transcripts (FPKM)") +
  facet_wrap(~condition)
dev.off()

cells <- c('A673','ASP14_day0','ASP14_day7','EW24','EW3','MHHES1','RDES','SKES1','TC71')
for (cell in cells){
  outPlot = c(paste(c(OUTDIR,"/KCNN1_expression_",cell,".pdf"),collapse=""))
  pdf(outPlot, width=5.9, height=4.7)
  print(ggplot(subset(data_boxplot,condition==cell), aes(x = nickname, y = expression, fill=transcript)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal() +
    guides(fill = "none", size="none") +
    xlab("KCNN1 transcripts") +
    ylab("Normalized expression of KCNN1 transcripts (FPKM)") +
    ggtitle(cell) +
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
}

palette_ASP14 <- c("#ece0c2","#090093")
names(palette_ASP14) <- c('ASP14_day0','ASP14_day7')
outPlot = c(paste(c(OUTDIR,"/KCNN1_expression_ASP14.pdf"),collapse=""))
pdf(outPlot, width=5.9, height=4.7)
print(ggplot(subset(data_boxplot,condition %in% c('ASP14_day0','ASP14_day7')), aes(x = nickname, y = expression, fill=condition)) +
  geom_boxplot(alpha=0.65) +
  labs(fill="Condition") +
  theme_minimal() +
  guides(size="none") +
  scale_fill_manual(values=palette_ASP14, labels=c('ASP14 day 0','ASP14 day 7')) +
  xlab("KCNN1 transcripts") +
  ylab("Normalized expression of KCNN1 transcripts (FPKM)") +
  ggtitle("ASP14") +
  theme(plot.title = element_text(hjust = 0.5)))
dev.off()

######## ASP14 comparison ########
condCol="condition"
logFCthreshold=1
AdjPValthreshold=0.05
GenesInFig=50
bootstrap=FALSE
nboot=30

cond1 <- "ASP14_day0"
cond2 <- "ASP14_day7"

comp <- paste(cond1,cond2,sep="__vs__")

sampleAnnot <- as.data.frame(sampleTable[,3])
rownames(sampleAnnot) <- sampleTable[,2]
colnames(sampleAnnot) <- "condition"
samples <- rownames(sampleAnnot)

res <- results(dds,contrast=c(condCol,cond1,cond2),independentFiltering=T)
res$meanInComp <- rowMeans(mrld[,sampleAnnot[,condCol]%in% c(cond1,cond2)])
res$padj <- adj_pvalues(ihw(pvalues = res$pvalue,covariates = res$meanInComp,alpha=AdjPValthreshold))

DE <- data.frame(res)

DE.sel <- list()
DE.sel$up <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange > logFCthreshold),]
DE.sel$down <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange < -logFCthreshold),]

DE.sel$isDE=rbind(DE.sel$up,DE.sel$down)
DE.sel$notDE=DE[setdiff(rownames(DE),rownames(DE.sel$isDE)),]

DE$DE="NONE"
DE[rownames(DE.sel$up),"DE"]="UP"
DE[rownames(DE.sel$down),"DE"]="DOWN"
DE$DE=factor(DE$DE,levels=c("DOWN","NONE","UP"))
DE$transcript <- rownames(DE)
correspondance_symbol <- read.table("/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/genomes/hg38/table_symbol_refseq.tsv",sep="\t",h=T)
DE <- merge(DE,correspondance_symbol)
write.table(DE, paste(paste(OUTDIR,comp,sep="/"),"DE.tsv",sep="_"),sep="\t",row.names=F,quote=F)


DE$DE_label <- NA
list_genes <- c("KCNN1")
DE$DE_label[which(DE$gene %in% list_genes)] <- DE$transcript[which(DE$gene %in% list_genes)]
DE_palette = c("firebrick3", "gray60", "palegreen3")
names(DE_palette) = c("DOWN","NONE","UP")

pdf(paste(c(OUTDIR,"/",comp,"_Volcano_with_names.pdf"),collapse=""),width=7.8,height=7.8)
  ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), color=DE, label=DE_label)) + 
  geom_point() +
  theme_minimal() +
  geom_text_repel(color="black", min.segment.length=0, nudge_x=2) +
  scale_color_manual(values=DE_palette)
dev.off()


