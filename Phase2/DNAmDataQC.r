## this script loads DNAm data from idat files
## requires a sample sheet to identify which sample is which 
## most of the output is saved to the folder QC
## NB fully methylated control is labelled as SampleID = "Blank"

library(methylumi)
library(wateRmelon)

setwd("")
idatPath<-"450K/iDATs"

## load sample sheet and create column with basename to match to idats
sampleSheet<-read.csv("", stringsAsFactors = FALSE)

sampleSheet<-sampleSheet[which(is.na(sampleSheet$MethArray_Chip) == FALSE),]
sampleSheet<-cbind(as.character(paste(sampleSheet$MethArray_Chip, sampleSheet$MethArray_ChipPosition, sep = "_")), sampleSheet)
colnames(sampleSheet)[1]<-"Basename"
sampleSheet$Basename<-as.character(sampleSheet$Basename)

## load data from idats
mset450k <- methylumIDAT(sampleSheet$Basename, idatPath=idatPath)
save(mset450k, file = ".Rdata") ## not essential to save but helpful if you need to revisit the QC at a later date.


### check intensities
m_intensities<-methylated(mset450k)
u_intensities<-unmethylated(mset450k)
## summarise to single value for each sample
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

write.csv(cbind(sampleSheet, M.median, U.median), "QC/MedianInsensitiesPerSample.csv")

chips<-unique(unlist(strsplit(sampleSheet$Basename, "_"))[seq(from = 1, to = 2*length(sampleSheet$Basename), by = 2)])

# to identify outlier batchs, summarise intensities by chip
chip.M.median<-aggregate(M.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(m_intensities), by = 2)]), FUN = median)
chip.U.median<-aggregate(U.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(u_intensities), by = 2)]), FUN = median)
chip.M.quantile<-aggregate(M.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(m_intensities), by = 2)]), FUN = quantile, probs = 0.66)
chip.U.quantile<-aggregate(U.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(u_intensities), by = 2)]), FUN = quantile, probs = 0.66)

write.csv(cbind(chip.M.median, chip.U.median, chip.M.quantile, chip.U.quantile), "QC/SummaryIntensitiesPerChip.csv")

pdf("QC/Plots/M_U_Intensities_RankedByIntensity.pdf", paper = "a4r", width = 0, height = 0)
par(mfcol = c(2, 4))
for(each in chips[order(chip.U.quantile[,2])]){
	index<-grep(each, colnames(m_intensities))
	boxplot(m_intensities[,index],  names = gsub(paste(each, "_", sep = ""), "", colnames(m_intensities)[index]), las = 2, main = each, ylab = "M intensity", xlab = "")
	abline(h = 3000, col = "red")
	boxplot(u_intensities[,index],  names = gsub(paste(each, "_", sep = ""), "", colnames(u_intensities)[index]), las = 2, main = each, ylab = "U intensity", xlab = "")
	abline(h = 3000, col = "red")
}
dev.off()


## in plot below samples are coloured by plate and the fully methylated control is highlighted a different plotting character
pdf("QC/Plots/Scatterplot_SampleIntensity.pdf")
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity")
abline(v = 2500, col = "red") ## suggested threshold for 450K arrays
hist(U.median, xlab = "Median U intensity")
abline(v = 2500, col = "red")
par(mfrow = c(1,1))
plot(M.median, U.median, xlab = "Median M intensity", ylab = "Median U intensity", col = factor(sampleSheet$MethArray_Plate), pch = c(16,17)[factor(sampleSheet$SampleID == "Blank")])
abline(v = 2500, col = "red")
abline(h = 2500, col = "red") 
legend("topleft", levels(factor(sampleSheet$MethArray_Plate)), col = palette(), pch = 16)
dev.off()

#Check fully methylated controls are where they should be
pdf("QC/Plots/FullyMethylatedControls.pdf")
par(mfrow = c(1,2))
boxplot(m_intensities[,which(sampleSheet$SampleID == "Blank")], names = sampleSheet$Plate[which(sampleSheet$SampleID == "Blank")], las = 2,  ylab = "M intensity", xlab = "", cex = 0.7)
abline(h = 3000, col = "red")
boxplot(u_intensities[,which(sampleSheet$SampleID == "Blank")],  names = sampleSheet$Plate[which(sampleSheet$SampleID == "Blank")], las = 2, ylab = "U intensity", xlab = "", cex = 0.7)
abline(h = 3000, col = "red")
dev.off()

# pull out fully methylated control based on ratio of intensities
ratio.intens<-M.median/U.median
write.csv(sampleSheet[which(ratio.intens > 2),], "QC/SamplesThatLookLikeFullyMethylatedControls.csv")

# exclude fully methylated controls from the rest of the QC
sampleSheet<-sampleSheet[-which(ratio.intens > 2),]
mset450k<-mset450k[,-which(ratio.intens > 2)]

# exclude samples with substandard intensity values
U.median<-U.median[match(sampleSheet$Basename, names(U.median))]
M.median<-M.median[match(sampleSheet$Basename, names(M.median))]

sampleSheet<-sampleSheet[-which(U.median < 2500 | M.median < 2500),]
mset450k<-mset450k[,-which(U.median < 2500 | M.median < 2500)]

### check bisulfite conversion
greenChannel<-intensitiesByChannel(QCdata(mset450k))$Cy3
redChannel<-intensitiesByChannel(QCdata(mset450k))$Cy5

### these are fully methylated controls so closer to 1/100% the better the BS conversion.
### NOTE due to differnet chemistrys type I conversion controls and type II converstion controls need to be handled differently.
greenChannel<-greenChannel[grep("BS.Conversion", rownames(greenChannel)),]
redChannel<-redChannel[grep("BS.Conversion", rownames(redChannel)),]
BScon1<-rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),])
BScon2<-redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),]

### for type II probes red channel is M and green channel is U so easy conversion into beta values
BScon2.betas<-redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),]/(redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),] + greenChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),])

### for type I both methylated and unmethylated in the same channel
BScon1.betas<-rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),])/(rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),]) + rbind(greenChannel[c("BS.Conversion.I.U1", "BS.Conversion.I.U2", "BS.Conversion.I.U3"),], redChannel[c("BS.Conversion.I.U4", "BS.Conversion.I.U5", "BS.Conversion.I.U6"),]))

BScon.betas<-rbind(BScon1.betas, BScon2.betas)
BScon.beta.median<-apply(BScon.betas, 2, median)*100

pdf("QC/Plots/Histogram_MedianBSConversion.pdf", width = 18, height = 8)
hist(BScon.beta.median, xlab = "Median % BS conversion")
dev.off()

sampleSheet<-cbind(sampleSheet, BScon.beta.median)

## exclude samples with incomplete bisulfite conversion
sampleSheet<-sampleSheet[-which(BScon.beta.median < 80),]
mset450k<-mset450k[,-which(BScon.beta.median < 80)]


### check predicted sex from data with that provided in the sample sheet 

## extract betas matrix
betas<-betas(mset450k)
sampleSheet<-sampleSheet[match(colnames(betas), sampleSheet$Basename),]
load("AllProbeIlluminaAnno.Rdata")
probeAnnot<-probeAnnot[rownames(betas),]

### X chr
probeAnnot.x<-probeAnnot[which(probeAnnot$CHR == "X"),]
betas.x<-betas[rownames(probeAnnot.x),]

d<-dist(t(betas.x))
fit<-cmdscale(d, k = 2)

pdf("QC/Plots/MDS_SexChromosomes_GenderChecks.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(1,2))
sex_palette<-cbind(c("M", "F"), c("blue", "magenta"))
sex_col<-sex_palette[match(sampleSheet$Sex, sex_palette[,1]),2]
sex_col[which(sampleSheet$Pool_ID == "neg")]<-"black"
plot(fit[,1], fit[,2], xlab = "MDS Co-ordinate 1", ylab = "MDS Co-ordinate 2", pch = 18, col = sex_col, main = "X Chromosome")
legend("topleft", horiz = TRUE, c("M", "F"), col = c("blue", "magenta"), pch = 18)

gender.x.chr<-vector(length = nrow(sampleSheet))
gender.x.chr[which(fit[,1] < 0)]<-"M"
gender.x.chr[which(fit[,1] > 0)]<-"F"


### Y chr
probeAnnot.y<-probeAnnot[which(probeAnnot$CHR == "Y"),]
betas.y<-betas[rownames(probeAnnot.y),]

d<-dist(t(betas.y))
fit<-cmdscale(d, k = 2)

plot(fit[,1], fit[,2], xlab = "MDS Co-ordinate 1", ylab = "MDS Co-ordinate 2", pch = 18, col = sex_col, main = "Y Chromosome")
legend("topleft", horiz = TRUE, c("M", "F"), col = c("blue", "magenta"), pch = 18)

dev.off()

gender.y.chr<-vector(length = nrow(sampleSheet))
gender.y.chr[which(fit[,1] > 0)]<-"F"
gender.y.chr[which(fit[,1] < 0)]<-"M"

sampleSheet<-cbind(sampleSheet, gender.x.chr, gender.y.chr)


## createfiles for upload to online age calculator. Crashes if you upload too many samples atonce so separate into separate files 
write.csv(betas[,1:138], "MethDataForAgeCalculator_1.csv")
write.csv(betas[,139:275], "MethDataForAgeCalculator_2.csv")
write.csv(betas[,276:413], "MethDataForAgeCalculator_3.csv")
write.csv(betas[,414:552], "MethDataForAgeCalculator_4.csv")
write.csv(betas[,552:nrow(betas)], "MethDataForAgeCalculator_5.csv")


pheno<-sampleSheet
female.indictor<-rep(0, nrow(pheno))
female.indictor[which(pheno$Sex == "F")]<-1
tissue<-rep("Blood", nrow(pheno))
 
pheno<-cbind(pheno[, c("Basename", "Age_at_Collection")], female.indictor, tissue)
colnames(pheno)<-c("SampleID", "Age", "Female", "Tissue")

write.csv(pheno[1:138,], "SampleDataForAgeCalculator_1.csv")
write.csv(pheno[139:275,], "SampleDataForAgeCalculator_2.csv")
write.csv(pheno[276:413,], "SampleDataForAgeCalculator_3.csv")
write.csv(pheno[414:552,], "SampleDataForAgeCalculator_4.csv")
write.csv(pheno[552:nrow(betas),], "SampleDataForAgeCalculator_5.csv")

## load output from age calculator
horvath<-read.csv("MethDataForAgeCalculator_1.output.csv", stringsAsFactors = FALSE, row.names = 1)
horvath<-rbind(horvath, read.csv("MethDataForAgeCalculator_2.output.csv", stringsAsFactors = FALSE, row.names = 1))
horvath<-rbind(horvath, read.csv("MethDataForAgeCalculator_3.output.csv", stringsAsFactors = FALSE, row.names = 1))
horvath<-rbind(horvath, read.csv("MethDataForAgeCalculator_4.output.csv", stringsAsFactors = FALSE, row.names = 1))
horvath<-rbind(horvath, read.csv("MethDataForAgeCalculator_5.output.csv", stringsAsFactors = FALSE, row.names = 1))

horvath<-horvath[match(sampleSheet$Basename, gsub("X", "", rownames(horvath))),]

## check predicted tissue

tmp<-table(horvath$predictedTissue)
pdf("QC/Plots/Boxplot_ProbabilityBloodTissues_QCedSample.pdf", width = 18, height = 8)
par(mfrow = c(2,4))
for(i in c(13:19, 49)){
	boxplot(horvath[,i] ~ horvath$predictedTissue, ylab = colnames(horvath)[i], xlab = "predictedTissue", names = paste(names(tmp), "\n n=", tmp, sep = ""), ylim = c(0,0.7))
}
dev.off()

### remove sample predicted as something other than blood

sampleSheet<-sampleSheet[-which(horvath$predictedTissue == "Saliva"),] ## extend list here as required
mset450k<-mset450k[,-which(horvath$predictedTissue == "Saliva")]


### compare to SNP chip data
geno<-read.table("Genotypes/450K_rsProbes.raw", header = TRUE, stringsAsFactors = FALSE)
geno<-geno[match(sampleSheet$GWAS_ID.1, geno$IID),]

betas<-betas(mset450k)

a<-NULL
for(i in 7:ncol(geno)){
	snp<-unlist(strsplit(colnames(geno)[i], "_"))[1]
	a<-append(a, grep(snp, rownames(betas)))
	}
	
meth.sub<-betas[a,]

### first check direction of minor alleles
cors<-vector(length = length(a))
pdf("QC/Plots/Raw_GenoAgainstMeth_PerSNP.pdf")
for(each in 1:nrow(meth.sub)){
	cors[each]<-cor(geno[,each+6], meth.sub[each,], use = "pairwise.complete.obs")
	plot(geno[,each+6], meth.sub[each,], xlab = "Genotype", ylab = "Methylation", main = rownames(meth.sub)[each], pch = 16,xlim = c(0,2), ylim = c(0,1), cex = 0.7)
	}
dev.off()

### change minor allele in genotype data if negative correlation
for(each in which(cors < 0)){
	geno[,each+6]<-(2-geno[,each+6])
}


pdf("QC/Plots/Raw_CompareMethylationwithGenotypeData_PerSample.pdf")
for(i in 1:ncol(meth.sub)){
	plot(as.numeric(geno[i,-c(1:6)]), meth.sub[,i], main = sampleSheet$GWAS_ID.1[match(colnames(meth.sub)[i], sampleSheet$Basename)], xlab = "Genotype", ylab = "Methylation", xlim = c(0,2), ylim = c(0,1), pch = 16, cex = 0.8)
	}
dev.off()

## after visualising the plots exclude samples that do not match
mset450k<-mset450k[,-which(sampleSheet$GWAS_ID.1 == c())]
sampleSheet<-sampleSheet[-which(sampleSheet$GWAS_ID.1 == c()),]


### run pfilter
mset450k.pf<-pfilter(mset450k)


### prenormalisation PCA on autosomal sites
meth<-betas(mset450k.pf)

probeAnnot<-probeAnnot[rownames(meth),]
probeAnnot<-probeAnnot[which(probeAnnot$CHR != "X"),]
probeAnnot<-probeAnnot[which(probeAnnot$CHR != "Y"),]
meth<-meth[rownames(probeAnnot),]

fit.pre <- princomp(meth[which(rowSums(is.finite(meth)) == ncol(meth)),], cor=TRUE)

pdf("QC/Plots/PCA_Prenormalisation_ColouredByPhenotype.pdf", height = 9, width = 9)
par(mfrow = c(2,2))
for(i in 2:9){
	plot(fit.pre$x[,1], fit.pre$x[,i], col = factor(sampleSheet$Sex), xlab = paste("PC1 : ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[1], 2), "%", sep = ""), ylab = paste("PC", i,": ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[i], 2), "%", sep = ""), pch = 16, main = "Coloured by Sex")

}
legend("bottomright", levels(factor(sampleSheet$Sex)), col = palette()[1:length(levels(factor(sampleSheet$Sex)))], pch = 16)
for(i in 2:9){
	plot(fit.pre$x[,1], fit.pre$x[,i], col = factor(sampleSheet$Status), xlab = paste("PC1 : ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[1], 2), "%", sep = ""), ylab = paste("PC", i,": ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[i], 2), "%", sep = ""), pch = 16, main = "Coloured by Disease Status")

}
legend("bottomright", levels(factor(sampleSheet$Status)), col = palette()[1:length(levels(factor(sampleSheet$Status)))], pch = 16)
for(i in 2:9){
	plot(fit.pre$x[,1], fit.pre$x[,i], col = factor(sampleSheet$Sentrix_ID), xlab = paste("PC1 : ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[1], 2), "%", sep = ""), ylab = paste("PC", i,": ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[i], 2), "%", sep = ""), pch = 16, main = "Coloured by Chip")

}
legend("bottomright", levels(factor(sampleSheet$Sentrix_ID)), col = palette()[1:length(levels(factor(sampleSheet$Sentrix_ID)))], pch = 16)
for(i in 2:9){
	plot(fit.pre$x[,1], fit.pre$x[,i], col = factor(sampleSheet$chip.location), xlab = paste("PC1 : ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[1], 2), "%", sep = ""), ylab = paste("PC", i,": ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[i], 2), "%", sep = ""), pch = 16, main = "Coloured by Chip Position")

}
legend("bottomright", levels(factor(sampleSheet$chip.location)), col = palette()[1:length(levels(factor(sampleSheet$chip.location)))], pch = 16)

age<-c(seq(21, 66, by = 1))
colfunc <- colorRampPalette(c("green", "red"))
colfunc(length(age))->gradient_colours
cbind(gradient_colours, age)->age_colours
match(round(as.numeric(as.character(sampleSheet$Age_at_Collection))), age_colours[,"age"])->test
age_colours[test,1]->cols_bar

for(i in 2:9){
	plot(fit.pre$x[,1], fit.pre$x[,i], col = cols_bar, xlab = paste("PC1 : ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[1], 2), "%", sep = ""), ylab = paste("PC", i,": ", signif((fit.pre$sdev^2 / sum(fit.pre$sdev^2)*100)[i], 2), "%", sep = ""), pch = 16, main = "Coloured by Reported Age")

}
dev.off()


## run dasen normalisation
mset450k.pf.dasen<-dasen(mset450k.pf)

## add cell composition variables and DNAmAge to sampleSheet
cellTypes<-c("Mono","Gran","NK","CD4T","CD8T","Bcell","PlasmaBlast", "CD8pCD28nCD45RAn","CD8.naive","CD4.naive")
sampleSheet<-cbind(sampleSheet, horvath[,cellTypes])

pdf("QC/Plots/AgeagainstDNAmAge_QCedSample.pdf")
plot(sampleSheet$Age_at_Collection, horvath$DNAmAge, xlab = "Reported Age at Collection", ylab = "DNAmAge", pch = 16, type = "n")
polygon(c(0,0,120,120), c(6,-6,114,127), col = "skyblue")
polygon(c(0,0,120,120), c(3,-3,117,123), col = "gray")
points(sampleSheet$Age_at_Collection, horvath$DNAmAge, pch = 16)
dev.off()

save(mset450k.pf.dasen, sampleSheet, file = "Normalised.rdat")