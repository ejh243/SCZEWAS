### Script to perform EWAS of schizophrenia case control status

library(wateRmelon)

setwd()
load() ### load normalised QCed R object

betas<-betas(mset450k.pf.dasen) ### extract betas matrix

pheno<-pheno[match(colnames(betas), pheno$Basename),]
crosshyb<-read.csv("/mnt/data1/450K_reference/CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("/mnt/data1/450K_reference/SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_EuAf_within10bpSBE == ""),]

betas<-betas[-grep("rs", rownames(betas)),]


### load DNA methylation smoking scores
smokingScores<-read.csv("", row.names = 1)
smokingScores<-smokingScores[match(pheno$Basename, smokingScores$Basename),]


res<-matrix(data = NA, nrow = nrow(betas), ncol = 25)
colnames(res)<-c("Beta", "SE", "P-value", "Anova:P", "Age:P", "Age:beta", "Sex:P", "Sex:beta:M",  "Smoking:P", "Smoking:beta", "CD8.naive:P", "CD8.naive:betas","CD8pCD28nCD45RAn:P", "CD8pCD28nCD45RAn:Beta","PlasmaBlast:P", "PlasmaBlast:Beta","CD4T:P","CD4T:Beta", "NK:P", "NK:Beta", "Mono:P","Mono:Beta", "Gran:P", "Gran:Beta", "Chip:P")
rownames(res)<-rownames(betas)
for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ pheno$STATUS + pheno$DNAmAge + factor(pheno$Sex) + smokingScores$scores_combined + pheno$CD8.naive + pheno$CD8pCD28nCD45RAn + pheno$PlasmaBlast + pheno$CD4T + pheno$NK + pheno$Mono + pheno$Gran + factor(pheno$Sentrix_ID))
	res[i,c(1,6,8,10,12,14,16,18,20,22,24)]<-coefficients(model)[2:12]
	res[i,c(4,25)]<-anova(model)[c("pheno$STATUS","factor(pheno$Sentrix_ID)"), "Pr(>F)"]
	res[i,2]<-summary(model)$coefficients[2,2]
	res[i,c(3,5,7,9,11,13,15,17,19,21,23)]<-summary(model)$coefficients[2:12,4]

}

### annotate
load("") ### load illumina probe annotation
probeAnnot<-probeAnnot[rownames(res),]
res<-cbind(res, probeAnnot[,c("CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP")])

smokingProbes<-read.csv("", stringsAsFactors = FALSE) ### load smoking EWAS results
smokingProbes.adj<-read.csv("", stringsAsFactors = FALSE) ### load smoking EWAS results (adjusted for cellular composition

res<-cbind(res, smokingProbes$EuropeanPvalue[match(rownames(res), smokingProbes$ProbeID)],smokingProbes$EuropeanEffectSize[match(rownames(res), smokingProbes$ProbeID)],smokingProbes.adj$EuropeanPvalue[match(rownames(res), smokingProbes.adj$ProbeID)],smokingProbes.adj$EuropeanEffectSize[match(rownames(res), smokingProbes.adj$ProbeID)])

colnames(res)[-c(1:29)]<-c("SmokingP", "SmokingEffect", "SmokingPCellCompAdj", "SmokingEffectCellCompAdj")
res<-res[order(res[,3]),]

### save output
write.csv(res, "")

