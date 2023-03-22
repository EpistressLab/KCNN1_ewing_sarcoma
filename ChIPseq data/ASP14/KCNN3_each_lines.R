
library(ggplot2)
library(gggenes)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

##################### Functions ###############

make_cov_plot <- function(data_cov,start,end){
  plot <- ggplot(data_cov, aes(x=position)) +
  theme_minimal() +
  theme(legend.position="none") %+replace% theme(
    axis.ticks.y = element_blank(),
    axis.title.y=element_text(size=7),
    axis.title.x=element_text(size=7)) +
  xlim(start,end) +
  xlab("Position on chr1") +
  geom_line(aes(y = reads_all), linewidth=0.5) +
  geom_area(aes(y = reads_all), linewidth=0.05,linetype = 0, fill="#bac5ff", alpha=0.5, outline.type="full")
  return(plot)
}

make_macs2_plot <- function(macs2_data,start,end,min_gradient,max_gradient){
    #macs2_data <- read.table(macs2_file, h=F, sep="\t")
    plot_macs2 <- ggplot(macs2_data) +
        geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
        xlim(start,end) +
        theme_minimal() +
        theme(legend.position="none") %+replace% theme(
            axis.ticks.y = element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_text(size=7),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        scale_fill_gradient(
            low = "#bac5ff",
            high = "#3b0e62",
            limits=c(min_gradient,max_gradient))
    return(plot_macs2)
}

make_macs2_legend <- function(macs2_data,start,end,min_gradient,max_gradient){
    #macs2_data <- read.table(macs2_file, h=F, sep="\t")
    plot_macs2 <- as_ggplot(get_legend(ggplot(macs2_data) +
        geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
        xlim(start,end) +
        theme_minimal() + 
        theme(legend.title = element_text(hjust = 0.5, size=8)) +
        scale_fill_gradient(
            low = "#bac5ff",
            high = "#3b0e62",
            limits=c(min_gradient,max_gradient)) +
        labs(fill=bquote(atop(bold("Macs2"),"\n-log10(qvalue)")))))
    return(plot_macs2)
}

###################### H3K4me3 KCNN3 plots ############

####### GTF file
#gtf <- read.table("/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/genomes/hg38/hg38_refseq_with_symbols_and_transcripts.gtf", sep="\t", h=F)
gtf <- read.table("~/bird/bone_epigenetics/anais/genomes/hg38/GCF_000001405.40_GRCh38.p14_genomic.gtf", sep="\t", h=F)
gtf <- gtf %>% separate('V9',into=c('a','b','c','d','e','f','g','h','i'),sep="; ",remove=T) %>% separate('a',into=c("gene_id","gene"),sep=" ") %>% separate('b',into=c("transcript_id","transcript"),sep=" ")
gtf <- gtf[,-c(2,9,11,13:19)]

gtf_KCNN3 <- subset(gtf,gene=="KCNN3")
gtf_KCNN3_transcripts <- subset(gtf_KCNN3, V3=='transcript')
gtf_KCNN3_transcripts <- gtf_KCNN3_transcripts[,-c(2,5,7)]
colnames(gtf_KCNN3_transcripts) <- c("chr","start_transcript","end_transcript","strand","gene","transcript")


gtf_KCNN3 <- subset(gtf_KCNN3, V3 %in% c("exon","CDS"))
gtf_KCNN3 <- merge(gtf_KCNN3,gtf_KCNN3_transcripts)
gtf_KCNN3$orientation <- 1
gtf_KCNN3$orientation[which(gtf_KCNN3$V7 == '-')] <- 0
gtf_KCNN3$transcript <- factor(gtf_KCNN3$transcript, levels = c("NM_170782.3","NM_001365838.1","NM_001365837.1","NM_001204087.2","NM_002249.6"))


plot_annot <- ggplot(gtf_KCNN3,aes(xmin=start_transcript, xmax=end_transcript, y=transcript, forward = orientation)) +
    geom_gene_arrow(fill="white") +
    geom_subgene_arrow(data=subset(gtf_KCNN3, V3=="exon"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=0.6, fill='coral1') +
    geom_subgene_arrow(data=subset(gtf_KCNN3, V3=="CDS"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=1, fill='coral4') +
    theme_genes() +
    geom_gene_label(aes(label = transcript)) +
    theme(legend.position="none") %+replace% theme(
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
    xlim(154685000,154880000) +
    ylab("KCNN3 transcripts")

######## Coverage plots

# H3K27ac
ASP14_d0_H3K27ac_cov <- read.table("./Samples/ASP14_d0_H3K27ac/ASP14_d0_H3K27ac_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d0_H3K27ac <- make_cov_plot(ASP14_d0_H3K27ac_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 0"),"Depth coverage")))
ASP14_d0_H3K27ac_macs2 <- subset(read.table("./Samples/ASP14_d0_H3K27ac/macs2/ASP14_d0_H3K27ac_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d7_H3K27ac_cov <- read.table("./Samples/ASP14_d7_H3K27ac/ASP14_d7_H3K27ac_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d7_H3K27ac <- make_cov_plot(ASP14_d7_H3K27ac_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 7"),"Depth coverage")))
ASP14_d7_H3K27ac_macs2 <- subset(read.table("./Samples/ASP14_d7_H3K27ac/macs2/ASP14_d7_H3K27ac_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")

all_peaks_H3K27ac <- rbind(ASP14_d0_H3K27ac_macs2,ASP14_d7_H3K27ac_macs2)
all_peaks_H3K27ac <- subset(all_peaks_H3K27ac, V2>=154685000 & V3 <= 154880000)
min_gradient <- floor(min(all_peaks_H3K27ac$V9))
max_gradient <- ceiling(max(all_peaks_H3K27ac$V9))

legend_macs2_H3K27ac <- make_macs2_legend(ASP14_d0_H3K27ac_macs2,154685000,154880000,min_gradient,max_gradient)
plot_macs2_ASP14_d0_H3K27ac <- make_macs2_plot(ASP14_d0_H3K27ac_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 0"),"Macs2")))
plot_macs2_ASP14_d7_H3K27ac <- make_macs2_plot(ASP14_d7_H3K27ac_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 7"),"Macs2")))

png("plots/KCNN3_ASP14_H3K27ac.png", width = 20, height = 25, units = "cm", res = 300)
plot_grid(plot_grid(plot_ASP14_d0_H3K27ac,plot_macs2_ASP14_d0_H3K27ac,
            plot_ASP14_d7_H3K27ac,plot_macs2_ASP14_d7_H3K27ac,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=5, rel_heights = c(rep(c(0.25, 0.10),2),0.3)),
        ggdraw(),legend_macs2_H3K27ac, ncol=3, rel_widths=c(0.85,0.05,0.15))
dev.off()

# FLI1
ASP14_d7_FLI1_cov <- read.table("./Samples/ASP14_d7_FLI1/ASP14_d7_FLI1_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d7_FLI1 <- make_cov_plot(ASP14_d7_FLI1_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 7"),"Depth coverage")))
ASP14_d7_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d7_FLI1/macs2/ASP14_d7_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d10_FLI1_cov <- read.table("./Samples/ASP14_d10_FLI1/ASP14_d10_FLI1_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d10_FLI1 <- make_cov_plot(ASP14_d10_FLI1_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 10"),"Depth coverage")))
ASP14_d10_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d10_FLI1/macs2/ASP14_d10_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d11_FLI1_cov <- read.table("./Samples/ASP14_d11_FLI1/ASP14_d11_FLI1_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d11_FLI1 <- make_cov_plot(ASP14_d11_FLI1_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 11"),"Depth coverage")))
ASP14_d11_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d11_FLI1/macs2/ASP14_d11_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d14_FLI1_cov <- read.table("./Samples/ASP14_d14_FLI1/ASP14_d14_FLI1_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d14_FLI1 <- make_cov_plot(ASP14_d14_FLI1_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 14"),"Depth coverage")))
ASP14_d14_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d14_FLI1/macs2/ASP14_d14_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")
ASP14_d17_FLI1_cov <- read.table("./Samples/ASP14_d17_FLI1/ASP14_d17_FLI1_KCNN3_coverage.tsv", h=T, sep="\t")
plot_ASP14_d17_FLI1 <- make_cov_plot(ASP14_d17_FLI1_cov,154685000,154880000) + ylab(bquote(atop(bold("Day 17"),"Depth coverage")))
ASP14_d17_FLI1_macs2 <- subset(read.table("./Samples/ASP14_d17_FLI1/macs2/ASP14_d17_FLI1_peaks.narrowPeak", h=F, sep="\t"),V1=="chr1")

all_peaks_FLI1 <- rbind(rbind(rbind(rbind(ASP14_d7_FLI1_macs2,ASP14_d10_FLI1_macs2),ASP14_d11_FLI1_macs2),ASP14_d14_FLI1_macs2),ASP14_d17_FLI1_macs2)
all_peaks_FLI1 <- subset(all_peaks_FLI1, V2>=154685000 & V3 <= 154880000)
min_gradient <- floor(min(all_peaks_FLI1$V9))
max_gradient <- ceiling(max(all_peaks_FLI1$V9))

legend_macs2_FLI1 <- make_macs2_legend(ASP14_d7_FLI1_macs2,154685000,154880000,min_gradient,max_gradient)
plot_macs2_ASP14_d7_FLI1 <- make_macs2_plot(ASP14_d7_FLI1_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 7"),"Macs2")))
plot_macs2_ASP14_d10_FLI1 <- make_macs2_plot(ASP14_d10_FLI1_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 10"),"Macs2")))
plot_macs2_ASP14_d11_FLI1 <- make_macs2_plot(ASP14_d11_FLI1_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 11"),"Macs2")))
plot_macs2_ASP14_d14_FLI1 <- make_macs2_plot(ASP14_d14_FLI1_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 14"),"Macs2")))
plot_macs2_ASP14_d17_FLI1 <- make_macs2_plot(ASP14_d17_FLI1_macs2,154685000,154880000,min_gradient,max_gradient) + ylab(bquote(atop(bold("Day 17"),"Macs2")))

png("plots/KCNN3_ASP14_FLI1.png", width = 20, height = 35, units = "cm", res = 300)
plot_grid(plot_grid(plot_ASP14_d7_FLI1,plot_macs2_ASP14_d7_FLI1,
            plot_ASP14_d10_FLI1,plot_macs2_ASP14_d10_FLI1,
            plot_ASP14_d11_FLI1,plot_macs2_ASP14_d11_FLI1,
            plot_ASP14_d14_FLI1,plot_macs2_ASP14_d14_FLI1,
            plot_ASP14_d17_FLI1,plot_macs2_ASP14_d17_FLI1,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=11, rel_heights = c(rep(c(0.18, 0.14),5),0.2)),
        ggdraw(),legend_macs2_H3K27ac, ncol=3, rel_widths=c(0.85,0.05,0.15))
dev.off()