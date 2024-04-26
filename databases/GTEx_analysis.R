######### Packages ########

library(ggplot2)
library(tidyverse)

######### Get data #########

palette <- c("#93b8fd","red")
names(palette) <- c("FALSE","TRUE")

exp <- read.table("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep="\t", h=T)
samples <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t", h=T, quote="")[,c(1,7)]

exp <- as.data.frame(t(exp[which(exp$Description == "KCNN1"),]))
exp$SAMPID <- gsub("\\.","-",rownames(exp))
exp <- exp[-c(1,2),]
exp$log2TPM_KCNN1 <- log(as.numeric(exp[,1])+1, 2)
exp <- merge(exp,samples)

###### Adding ewing sarcoma data #####

exp_CCLE <- read.table("../CCLE_Depmap_22Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", h=T)
samples_CCLE <- read.table("../CCLE_Depmap_22Q4/Model.csv", sep=",", h=T)[,c(1,23)]

colnames(exp_CCLE) <- sapply(strsplit(colnames(exp_CCLE), "..", fixed=T),'[[',1)
colnames(exp_CCLE)[1] <- "ModelID"

exp_CCLE <- exp_CCLE[,colnames(exp_CCLE) %in% c("ModelID","KCNN1")]
data_exp_CCLE <- merge(exp_CCLE,samples_CCLE)
colnames(data_exp_CCLE) <- c("SAMPID","log2TPM_KCNN1","SMTSD")
data_exp <- rbind(exp[,c(1,3,4)],data_exp_CCLE[which(data_exp_CCLE$SMTSD == "Ewing Sarcoma"),])
data_exp$is_ewing <- (data_exp$SMTSD == "Ewing Sarcoma")

pdf("GTEx_boxplot_KCNN1_with_ewing.pdf",height=10,width=5.7)
ggplot(data_exp, aes(x=log2TPM_KCNN1, y=reorder(SMTSD,log2TPM_KCNN1,decreasing=F), fill=is_ewing)) +
    geom_boxplot() +
    xlab(expression(atop(italic(KCNN1)~" expression","log2(TPM+1)"))) +
    ylab("Tissue Site Detail field") +
     scale_fill_manual(values=palette) +
    theme_minimal() +
    ggtitle("GTEx version 8 and\nEwing Sarcoma DepMap 22Q4") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(fill="none")
dev.off()