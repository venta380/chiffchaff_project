setwd("/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files")
library("SNPRelate")
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
vcf.fn<-"/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files/common.vcf_used_for_structure_PCA.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table("/proj/b2013146/nobackup/analysis/Faheema_analyis/vcftools/VCF_files/pop.group", header = F)


#ccm_pca<-snpgdsPCA(genofile)
#names(ccm_pca)
head(cbind(sample.id, pop_code))

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
names(snpset)
head(snpset$chr1)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=1)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sample.id = pca$sample.id
pop = factor(pop_code$V1)

tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code$V1)[match(pca$sample.id, sample.id)],colour=factor(pop_code$V2)[match(pca$sample.id, sample.id)] , EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab)

pdf(file= './chiff_PCA_leg' ,onefile=T,paper='A4', )
plot(tab$EV2~tab$EV1, col=alpha(as.character(tab$colour), 0.8), xlab="Eigenvector 2", ylab="Eigenvector 1",pch=19,cex = 2, cex.lab=2.5, par(mar =c(5,6,5.5,5.5)))
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
legend("topright",xpd = TRUE,legend=c("P. c. abietinus (allopatry)","P. c. abietinus (sympatric)","P. tristis (sympatry)","P. tristis (allopatry)"),text.font=3, pch=19, col=alpha(c("#FFE200","#285fba", "#28BA57", "#bf5509"), 0.8), cex=2)
dev.off()

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.character(tab$colour), labels=lbls, pch=19)

