
library(ggplot2)
library(gggenes)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

##################### Functions ###############

make_cov_plot <- function(data_cov,start,end,max_y){
  plot <- ggplot(data_cov, aes(x=position)) +
  theme_minimal() +
  theme(legend.position="none") %+replace% theme(
    axis.ticks.y = element_blank(),
    axis.title.y=element_text(size=7),
    axis.title.x=element_text(size=7)) +
  xlim(start,end) + 
  ylim(0,max_y) +
  xlab("Position on chr19") +
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

###################### H3K4me3 KCNN1 plots ############

####### GTF file
#gtf <- read.table("/LAB-DATA/BiRD/shares/Phy-OS/bone_epigenetics/anais/genomes/hg38/hg38_refseq_with_symbols_and_transcripts.gtf", sep="\t", h=F)
gtf <- read.table("~/bird/bone_epigenetics/anais/genomes/hg38/GCF_000001405.40_GRCh38.p14_genomic.gtf", sep="\t", h=F)
gtf <- gtf %>% separate('V9',into=c('a','b','c','d','e','f','g','h','i'),sep="; ",remove=T) %>% separate('a',into=c("gene_id","gene"),sep=" ") %>% separate('b',into=c("transcript_id","transcript"),sep=" ")
gtf <- gtf[,-c(2,9,11,13:19)]

gtf_KCNN1 <- subset(gtf,gene=="KCNN1")
gtf_KCNN1_transcripts <- subset(gtf_KCNN1, V3=='transcript')
gtf_KCNN1_transcripts <- gtf_KCNN1_transcripts[,-c(2,5,7)]
colnames(gtf_KCNN1_transcripts) <- c("chr","start_transcript","end_transcript","strand","gene","transcript")


gtf_KCNN1 <- subset(gtf_KCNN1, V3 %in% c("exon","CDS"))
gtf_KCNN1 <- merge(gtf_KCNN1,gtf_KCNN1_transcripts)
gtf_KCNN1$orientation <- 1
gtf_KCNN1$orientation[which(gtf_KCNN1$V7 == '-')] <- 0
gtf_KCNN1$transcript <- factor(gtf_KCNN1$transcript, levels = c("XM_011528004.2","NR_170374.1","NR_170373.1","NM_001386974.1","NM_001386975.1","NM_001386977.1","NM_001386976.1","NM_002248.5"))
transcripts_nicknames <- data.frame('transcript'=c("XM_011528004.2","NR_170374.1","NR_170373.1","NM_001386974.1","NM_001386975.1","NM_001386977.1","NM_001386976.1","NM_002248.5"),'nickname'=c('H','G','F','E','D','C','B','A'))
gtf_KCNN1 <- merge(gtf_KCNN1,transcripts_nicknames)
gtf_KCNN1$nickname <- factor(gtf_KCNN1$nickname, levels=c('H','G','F','E','D','C','B','A'))

png("KCNN1_transcripts.png", width = 20, height = 12, units = "cm", res = 300)
ggplot(gtf_KCNN1,aes(xmin=start_transcript, xmax=end_transcript, y=nickname, forward = orientation)) +
    geom_gene_arrow(fill="white") +
    geom_subgene_arrow(data=subset(gtf_KCNN1, V3=="exon"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=0.6, fill='coral1') +
    geom_subgene_arrow(data=subset(gtf_KCNN1, V3=="CDS"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=1, fill='coral4') +
    theme_genes() +
    geom_gene_label(aes(label = transcript)) +
    theme(legend.position="none") %+replace% theme(
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
    xlim(17940000,18010200) +
    xlab("KCNN1 transcripts")
dev.off()

plot_annot <- ggplot(gtf_KCNN1,aes(xmin=start_transcript, xmax=end_transcript, y=transcript, forward = orientation)) +
    geom_gene_arrow(fill="white") +
    geom_subgene_arrow(data=subset(gtf_KCNN1, V3=="exon"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=0.6, fill='coral1') +
    geom_subgene_arrow(data=subset(gtf_KCNN1, V3=="CDS"),aes(xsubmin=V4,xsubmax=V5,fill=gene),color="black", alpha=1, fill='coral4') +
    theme_genes() +
    geom_gene_label(aes(label = nickname)) +
    theme(legend.position="none") %+replace% theme(
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
    xlim(17940000,18010200) +
    ylab("KCNN1 transcripts")

######## Coverage plots

cells <- c('A673','CHLA10','CHLA25','EW1','EW3','EW7','EW22','EW24','MHHES1','MIC','POE','RDES','RH1','SKES1','SKNMC','TC32','TC71','TC106')
GGAA_coords <- data.frame('start'=c(17966024),'end'=c(17966083))

for (cell in cells){
    ## Data reading
    H3K4me3_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K4me3_KCNN1_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K4me3_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K4me3_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr19")
    H3K27ac_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K27ac_KCNN1_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K27ac_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K27ac_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr19")
    H3K27me3_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_H3K27me3_KCNN1_coverage.tsv"),collapse=""), h=T, sep="\t")
    H3K27me3_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_H3K27me3_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr19")
    transcription_factor_cov <- read.table(paste(c("./Samples/",cell,"/",cell,"_transcription_factor_KCNN1_coverage.tsv"),collapse=""), h=T, sep="\t")
    transcription_factor_macs2 <- subset(read.table(paste(c("./Samples/",cell,"/macs2/",cell,"_transcription_factor_peaks.narrowPeak"),collapse=""), h=F, sep="\t"),V1=="chr19")

    ## For a plot with all sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,H3K27me3_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_H3K27ac_cov <- make_cov_plot(H3K27ac_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27ac"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_H3K27me3_cov <- make_cov_plot(H3K27me3_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27me3"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,17940000,18010200,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    else {plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(rbind(rbind(H3K4me3_macs2,transcription_factor_macs2),H3K27ac_macs2),H3K27me3_macs2)
    all_peaks <- subset(all_peaks, V2>=17940000 & V3 <= 18010200)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))
    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_H3K27ac <- make_macs2_plot(H3K27ac_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27ac"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_H3K27me3 <- make_macs2_plot(H3K27me3_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,17940000,18010200,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    else {plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    png(paste(c("plots/KCNN1_",cell,"_all_seq.png"),collapse=""), width = 20, height = 30, units = "cm", res = 300)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_H3K27ac_cov, plot_macs2_H3K27ac,
                plot_H3K27me3_cov, plot_macs2_H3K27me3,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=10, rel_heights = c(0.04,rep(c(0.12, 0.07),4),0.2)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()

    ## For a plot with only transcription factor, H3K4me3 and H3K27ac sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_H3K27ac_cov <- make_cov_plot(H3K27ac_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K27ac"),"Depth coverage")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,17940000,18010200,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    else {plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(rbind(H3K4me3_macs2,transcription_factor_macs2),H3K27ac_macs2)
    all_peaks <- subset(all_peaks, V2>=17940000 & V3 <= 18010200)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))

    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_H3K27ac <- make_macs2_plot(H3K27ac_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K27ac"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,17940000,18010200,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    else {plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    png(paste(c("plots/KCNN1_",cell,"_three_seq.png"),collapse=""), width = 20, height = 28, units = "cm", res = 300)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_H3K27ac_cov, plot_macs2_H3K27ac,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=8, rel_heights = c(0.06,rep(c(0.15, 0.08),3),0.25)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()

    ## For a plot with only transcription factor and H3K4me3 sequencings
    # Coverage plots
    max_cov_plot <- max(c(H3K4me3_cov$reads_all,H3K27ac_cov$reads_all,transcription_factor_cov$reads_all))
    plot_H3K4me3_cov <- make_cov_plot(H3K4me3_cov,17940000,18010200,max_cov_plot) + 
        ylab(bquote(atop(bold("H3K4me3"),"Depth coverage"))) +
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_transcription_factor_cov <- make_cov_plot(transcription_factor_cov,17940000,18010200,max_cov_plot)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=max_cov_plot),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("ERG"),"Depth coverage")))}
    else {plot_transcription_factor_cov <- plot_transcription_factor_cov + ylab(bquote(atop(bold("FLI1"),"Depth coverage")))}

    # Macs2 plots
    all_peaks <- rbind(H3K4me3_macs2,transcription_factor_macs2)
    all_peaks <- subset(all_peaks, V2>=17940000 & V3 <= 18010200)
    min_gradient <- floor(min(all_peaks$V9))
    max_gradient <- ceiling(max(all_peaks$V9))

    legend_macs2 <- make_macs2_legend(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient)
    plot_macs2_H3K4me3 <- make_macs2_plot(H3K4me3_macs2,17940000,18010200,min_gradient,max_gradient) + 
        ylab(bquote(atop(bold("H3K4me3"),"Macs2")))+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    plot_macs2_transcription_factor <- make_macs2_plot(transcription_factor_macs2,17940000,18010200,min_gradient,max_gradient)+
        geom_rect(data=GGAA_coords, aes(xmin=start,xmax=end,ymin=0,ymax=1),inherit.aes=F, fill="#f2c025", alpha=0.5)
    if (cell %in% c('CHLA25','EW3','TC106')){plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("ERG"),"Macs2")))}
    else {plot_macs2_transcription_factor <- plot_macs2_transcription_factor + ylab(bquote(atop(bold("FLI1"),"Macs2")))}

    # Merging plots
    png(paste(c("plots/KCNN1_",cell,"_two_seq.png"),collapse=""), width = 20, height = 25, units = "cm", res = 300)
    print(plot_grid(plot_grid(ggplot() + annotate("text", x = 1, y = 1, size=5, label=cell) + theme_void(),
                plot_H3K4me3_cov,plot_macs2_H3K4me3,
                plot_transcription_factor_cov,plot_macs2_transcription_factor,
                plot_annot, align = "v", axis="tb", ncol=1, nrow=5, rel_heights = c(0.05,rep(c(0.25, 0.075),2),0.3)),
            ggdraw(),legend_macs2, ncol=3, rel_widths=c(0.85,0.05,0.15)))
    dev.off()
}

