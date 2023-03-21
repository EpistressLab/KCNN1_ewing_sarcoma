
library(ggplot2)
library(gggenes)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

##################### Functions ###############

get_density_pos <- function(pos,data) {
  density_entry <- data[data$V2 <= pos & data$V3 > pos, "V4"]
  if (length(density_entry) == 1) return(density_entry)
  else return(NA)
}

make_df_cov <- function(file_path,chrom,start,end){
    data_density <- subset(read.table(file_path, sep="\t", h=F),V1 == chrom)
    all_density <- sapply(start:end, get_density_pos, data=data_density)
    return(data.frame("chromosome"=rep(chrom,length(start:end)), position=start:end, reads_all=all_density))
}

make_cov_plot <- function(data_cov,start,end,color_fill){
  plot <- ggplot(data_cov, aes(x=position)) +
  theme_minimal() +
  theme(legend.position="none") %+replace% theme(
    axis.ticks.y = element_blank(),
    axis.title.y=element_text(size=7),
    axis.title.x=element_text(size=7)) +
  xlim(start,end) +
  xlab("Position on chr1") +
  geom_line(aes(y = reads_all), linewidth=0.5) +
  geom_area(aes(y = reads_all), linewidth=0.05,linetype = 0, fill=color_fill, alpha=0.5, outline.type="full")
  return(plot)
}

# make_macs2_plot <- function(macs2_data,start,end,min_gradient,max_gradient){
#     #macs2_data <- read.table(macs2_file, h=F, sep="\t")
#     plot_macs2 <- ggplot(macs2_data) +
#         geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
#         xlim(start,end) +
#         theme_minimal() +
#         theme(legend.position="none") %+replace% theme(
#             axis.ticks.y = element_blank(),
#             axis.text.y=element_blank(),
#             axis.title.y=element_text(size=7),
#             axis.line.x = element_blank(),
#             axis.ticks.x = element_blank(),
#             axis.text.x = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank()) +
#         scale_fill_gradient(
#             low = "#bac5ff",
#             high = "#3b0e62",
#             limits=c(min_gradient,max_gradient))
#     return(plot_macs2)
# }

# make_macs2_legend <- function(macs2_data,start,end,min_gradient,max_gradient){
#     #macs2_data <- read.table(macs2_file, h=F, sep="\t")
#     plot_macs2 <- as_ggplot(get_legend(ggplot(macs2_data) +
#         geom_rect(aes(xmin=V2,xmax=V3,ymin=0,ymax=1,fill=V9),alpha=0.7) +
#         xlim(start,end) +
#         theme_minimal() + 
#         theme(legend.title = element_text(hjust = 0.5, size=8)) +
#         scale_fill_gradient(
#             low = "#bac5ff",
#             high = "#3b0e62",
#             limits=c(min_gradient,max_gradient)) +
#         labs(fill=bquote(atop(bold("Macs2"),"\n-log10(qvalue)")))))
#     return(plot_macs2)
# }

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

# plasmid control
MSC_control_FLI1_data <- make_df_cov("./Samples/MSC1_control_FLI1/bedgraphs/MSC1_control_FLI1_hg38.bedGraph","chr1",154685000,154880000)
MSC_control_FLI1 <- make_cov_plot(MSC_control_FLI1_data,154685000,154880000,"#85d1e4") + ylab(bquote(atop(bold("Plasmid control"),"FLI1")))
MSC_control_H3K27ac_data <- make_df_cov("./Samples/MSC1_control_H3K27ac/bedgraphs/MSC1_control_H3K27ac_hg38.bedGraph","chr1",154685000,154880000)
MSC_control_H3K27ac <- make_cov_plot(MSC_control_H3K27ac_data,154685000,154880000,"#85d1e4") + ylab(bquote(atop(bold("Plasmid control"),"H3K27ac")))
MSC_control_input_data <- make_df_cov("./Samples/MSC1_control_input/bedgraphs/MSC1_control_input_hg38.bedGraph","chr1",154685000,154880000)
MSC_control_input <- make_cov_plot(MSC_control_input_data,154685000,154880000,"#85d1e4") + ylab(bquote(atop(bold("Plasmid control"),"Input")))

plot_control <- plot_grid(MSC_control_FLI1 ,MSC_control_H3K27ac,MSC_control_input, align = "v", axis="tb", ncol=1, nrow=3, rel_heights = c(0.33,0.33,0.33))
png("./plots/KCNN3_MSC1_plasmControl.png", width = 20, height = 35, units = "cm", res = 300)
plot_control
dev.off()

# plasmid EWSFLI1
MSC_EWSFLI1_FLI1_data <- make_df_cov("./Samples/MSC1_EWSFLI1_FLI1/bedgraphs/MSC1_EWSFLI1_FLI1_hg38.bedGraph","chr1",154685000,154880000)
MSC_EWSFLI1_FLI1 <- make_cov_plot(MSC_EWSFLI1_FLI1_data,154685000,154880000,"#ec93a6") + ylab(bquote(atop(bold("Plasmid EWSFLI1"),"FLI1")))
MSC_EWSFLI1_H3K27ac_data <- make_df_cov("./Samples/MSC1_EWSFLI1_H3K27ac/bedgraphs/MSC1_EWSFLI1_H3K27ac_hg38.bedGraph","chr1",154685000,154880000)
MSC_EWSFLI1_H3K27ac <- make_cov_plot(MSC_EWSFLI1_H3K27ac_data,154685000,154880000,"#ec93a6") + ylab(bquote(atop(bold("Plasmid EWSFLI1"),"H3K27ac")))
MSC_EWSFLI1_input_data <- make_df_cov("./Samples/MSC1_EWSFLI1_input/bedgraphs/MSC1_EWSFLI1_input_hg38.bedGraph","chr1",154685000,154880000)
MSC_EWSFLI1_input <- make_cov_plot(MSC_EWSFLI1_input_data,154685000,154880000,"#ec93a6") + ylab(bquote(atop(bold("Plasmid EWSFLI1"),"Input")))

plot_EWSFLI1 <- plot_grid(MSC_EWSFLI1_FLI1 ,MSC_EWSFLI1_H3K27ac,MSC_EWSFLI1_input, align = "v", axis="tb", ncol=1, nrow=3, rel_heights = c(0.33,0.33,0.33))
png("./plots/KCNN3_MSC1_plasmEWSFLI1.png", width = 20, height = 35, units = "cm", res = 300)
plot_EWSFLI1
dev.off()


# plasmid FLI1
MSC_FLI1_FLI1_data <- make_df_cov("./Samples/MSC1_FLI1_FLI1/bedgraphs/MSC1_FLI1_FLI1_hg38.bedGraph","chr1",154685000,154880000)
MSC_FLI1_FLI1 <- make_cov_plot(MSC_FLI1_FLI1_data,154685000,154880000,"#f49b70") + ylab(bquote(atop(bold("Plasmid FLI1"),"FLI1")))
MSC_FLI1_H3K27ac_data <- make_df_cov("./Samples/MSC1_FLI1_H3K27ac/bedgraphs/MSC1_FLI1_H3K27ac_hg38.bedGraph","chr1",154685000,154880000)
MSC_FLI1_H3K27ac <- make_cov_plot(MSC_FLI1_H3K27ac_data,154685000,154880000,"#f49b70") + ylab(bquote(atop(bold("Plasmid FLI1"),"H3K27ac")))
MSC_FLI1_input_data <- make_df_cov("./Samples/MSC1_FLI1_input/bedgraphs/MSC1_FLI1_input_hg38.bedGraph","chr1",154685000,154880000)
MSC_FLI1_input <- make_cov_plot(MSC_FLI1_input_data,154685000,154880000,"#f49b70") + ylab(bquote(atop(bold("Plasmid FLI1"),"Input")))

plot_FLI1 <- plot_grid(MSC_FLI1_FLI1 ,MSC_FLI1_H3K27ac,MSC_FLI1_input, align = "v", axis="tb", ncol=1, nrow=3, rel_heights = c(0.33,0.33,0.33))
png("./plots/KCNN3_MSC1_plasmFLI1.png", width = 20, height = 35, units = "cm", res = 300)
plot_FLI1
dev.off()

png("./plots/KCNN3_MSC1.png", width = 20, height = 35, units = "cm", res = 300)
plot_grid(MSC_control_FLI1 ,MSC_control_H3K27ac,
            MSC_FLI1_FLI1 ,MSC_FLI1_H3K27ac,
            MSC_EWSFLI1_FLI1 ,MSC_EWSFLI1_H3K27ac,
            plot_annot, align = "v", axis="tb", ncol=1, nrow=7, rel_heights = c(rep(c(0.125, 0.125),3),0.25))
dev.off()