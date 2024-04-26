######### Packages ########

library(ggplot2)
library(tidyverse)

######### Get data #########

palette <- c("#93b8fd","red")
names(palette) <- c("FALSE","TRUE")

exp <- read.table("OmicsExpressionProteinCodingGenesTPMLogp1.csv", sep=",", h=T)
samples <- read.table("Model.csv", sep=",", h=T)[,c(1,23)]

colnames(exp) <- sapply(strsplit(colnames(exp), "..", fixed=T),'[[',1)
colnames(exp)[1] <- "ModelID"

exp <- exp[,colnames(exp) %in% c("ModelID","KCNN1")]
data_exp <- merge(exp,samples)

######## Plot #########

data_exp$is_ewing <- (data_exp$OncotreePrimaryDisease == "Ewing Sarcoma")

pdf("CCLE_boxplot_KCNN1.pdf",height=10,width=5.9)
ggplot(data_exp, aes(x=KCNN1, y=reorder(OncotreePrimaryDisease,KCNN1,decreasing=F), fill=is_ewing)) +
    geom_boxplot() +
    xlab(expression(atop(italic(KCNN1)~" expression","log2(TPM+1)"))) +
    ylab("Primary disease") + 
    scale_fill_manual(values=palette) +
    theme_minimal() +
    guides(fill="none") +
    ggtitle("DepMap 22Q4")
dev.off()