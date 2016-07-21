### Script to perform EWAS of schizophrenia polygenic risk score

### Copyright (C) 2016  Eilis Hannon

##    This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License along
##    with this program; if not, write to the Free Software Foundation, Inc.,
##    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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

### load polygenic risk scores
scores<-read.table("", header = TRUE, row.names = 1)
scores<-scores[match(pheno$V2, scores$IID),]

res<-matrix(data = NA, nrow = nrow(betas), ncol = 3)
colnames(res)<-c("Beta", "SE", "P-value")
rownames(res)<-rownames(betas)
for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ scores$SCORE + pheno$DNAmAge + factor(pheno$Sex) + smokingScores$scores_combined + pheno$CD8.naive + pheno$CD8pCD28nCD45RAn + pheno$PlasmaBlast + pheno$CD4T + pheno$NK + pheno$Mono + pheno$Gran + factor(pheno$Sentrix_ID))
	res[i,1]<-coefficients(model)[2]
	res[i,2]<-summary(model)$coefficients[2,2]
	res[i,3]<-summary(model)$coefficients[2,4]

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

