#install.packages(c("fields","RColorBrewer","mapplots"))
#paste -d "" 1 2 3 4 5 6 7 8 9 10 22 23 24 25 26 27 28 29 30 31 18 19 20 21 32 33 34 35 36 37 39 38 11 12 13 14 15 16 17 > ../plink.geno

source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
install.packages("RColorBrewer")

library(RColorBrewer)

library(LEA)

jBrewColors <- brewer.pal(n = 10, name = "Paired")
V3<-c("Pcol01","Pcol02","Pcol03","Pcol04","Pcol06","Pcol07","Pcol09","Pcol10","Pcol11","Pcol12","Pcol13","Pcol14","Pcol15","Pcol16","Pcol17","Pcol18","Pcol19","Pcol20","Pcol21","Pcol22","Pcol23","Pcol24","Pcol25","Pcol26","Pcol27","Pcol29","Pcol31","Pcol32","Pcol33","Pcol34","Pcol35","Pcol37","Pcol38","Pcol39","Pcol41","Pcol42","Pcol47","Pcol510","Pcol704")


pdf(file= './chiff_struct' ,onefile=T,paper='A4', )

jBrewColors <- c('#FFE200','#bf5509')
obj.snmf = snmf("plink.geno", K = 1:8, ploidy = 2, entropy = T,alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)


obj.snmf = snmf("plink.geno", K = 2,  project = "new")
  qmatrix1 = Q(obj.snmf, K = 2)
plot1=barplot(t(qmatrix1), col = jBrewColors, border = NA, ylab = "K=2")

dev.off()

######################################end of ploting###############

obj.snmf1 = snmf("plink.geno", K = 3,  project = "new")
  qmatrix1 = Q(obj.snmf1, K = 3)
   #qmat2 <- qmatrix1[order(qmatrix1[,1]),]
plot1=barplot(t(qmatrix1), col = jBrewColors, border = NA, ylab = "K=3")


obj.snmf2 = snmf("plink.geno", K = 4,  project = "new")
  qmatrix2 = Q(obj.snmf2, K = 4)
  #qmat2 <- qmatrix1[order(qmatrix1[,1]),]
plot1=barplot(t(qmatrix2), col = jBrewColors, border = NA, ylab = "K=4")


obj.snmf3 = snmf("plink.geno", K = 5,  project = "new")
  qmatrix3 = Q(obj.snmf3, K = 5)
  #qmat2 <- qmatrix1[order(qmatrix1[,1]),]
plot1=barplot(t(qmat2), col = jBrewColors, border = NA, ylab = "K=5")

obj.snmf4 = snmf("plink.geno", K = 6,  project = "new")
  qmatrix4 = Q(obj.snmf4, K = 6)
  qmat2 <- qmatrix1[order(qmatrix1[,1]),]
plot1=barplot(t(qmat2), col = jBrewColors, border = NA, ylab = "K=6")

obj.snmf5 = snmf("plink.geno", K = 7,  project = "new")
  qmatrix5 = Q(obj.snmf5, K = 7)
plot5=barplot(t(qmatrix5), col = jBrewColors, border = NA, ylab = "K=7")

obj.snmf6 = snmf("plink.geno", K = 6,  project = "new")
  qmatrix6 = Q(obj.snmf6, K = 6)
plot6=barplot(t(qmatrix6), col = jBrewColors, border = NA, ylab = "K=8")


obj.snmf7 = snmf("plink.geno", K = 9,  project = "new")
  qmatrix7 = Q(obj.snmf7, K = 9)
plot7=barplot(t(qmatrix7), col = jBrewColors, border = NA, ylab = "K=9")


