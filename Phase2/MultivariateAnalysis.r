## Script to take CpGs associated with schizophrenia and repeat anlaysis with covariate for TRS.
## meta-analysis of both SCZ vs CON and TRS vs non-TRS effects

library(meta)
library(grid)
library(wateRmelon)

## load data for each cohort and resave pheno dataframe as separate R object called pheno.<cohortName>
setwd("UCL")
load("Normalised.rdat")
pheno.ucl<-pheno
betas.ucl<-betas

setwd("Aberdeen")
load("Normalised.rdat")
pheno.aber<-pheno
betas.aber<-betas

setwd("IOP")
load("Normalised.rdat")
pheno.iop<-pheno
betas.iop<-betas

setwd("Sweden")
load("Normalised.rdat")
pheno.swedish<-pheno
betas.swedish<-betas


## load schizophrenia case control EWAS results
res.meta<-read.csv("", stringsAsFactors = FALSE, row.names = 1)
## identify significant DMPs
index<-which(res.meta$All_P_Random < 9e-8) 

tabmeta<-matrix(data = NA, nrow = length(index), ncol = 21)
rownames(tabmeta)<-rownames(res.meta)[index]
colnames(tabmeta)<-c("N_cohorts", paste("SCZ", c("Effect_Fixed", "Effect_SE_Fixed", "P_Fixed", "Effect_Random", "Effect_SE_Random","P_Random", "tau", "I2", "Q", "Het P")), paste("TRS", c("Effect_Fixed", "Effect_SE_Fixed", "P_Fixed", "Effect_Random", "Effect_SE_Random","P_Random", "tau", "I2", "Q", "Het P")))

## perform EWAS in each cohort testing SCZ vs CON and TRS vd non-TRS controls
for(i in index){
	cohortCoef<-matrix(data = NA, ncol = 4, nrow = 4)
	probeID<-rownames(res.meta)[i]
	if(probeID %in% rownames(betas.ucl)){	
		model.ucl<-lm(betas.ucl[probeID,] ~ pheno.ucl$cloz + pheno.ucl$STATUS + pheno.ucl$DNAmAge + factor(pheno.ucl$Sex) + pheno.ucl$smokingScore + pheno.ucl$CD8.naive + pheno.ucl$CD8pCD28nCD45RAn + pheno.ucl$PlasmaBlast + pheno.ucl$CD4T + pheno.ucl$NK + pheno.ucl$Mono + pheno.ucl$Gran + factor(pheno.ucl$Sentrix_ID))
		cohortCoef[1,]<-c(summary(model.ucl)$coefficients[2,c(1,2)], summary(model.ucl)$coefficients[3,c(1,2)])
	} 
	if(probeID %in% rownames(betas.aber)){	
		model.aber<-lm(betas.aber[probeID,] ~ pheno.aber$cloz + pheno.aber$Status + pheno.aber$DNAmAge + factor(pheno.aber$gender.x.chr) + pheno.aber$DNAmSmokingScore + pheno.aber$CD8.naive + pheno.aber$CD8pCD28nCD45RAn + pheno.aber$PlasmaBlast + pheno.aber$CD4T + pheno.aber$NK + pheno.aber$Mono + pheno.aber$Gran + factor(pheno.aber$MethArray_Chip))
		cohortCoef[2,]<-c(summary(model.aber)$coefficients[2,c(1,2)], summary(model.aber)$coefficients[3,c(1,2)])
	}
	if(probeID %in% rownames(betas.iop)){	
		model.iop<-lm(betas.iop[probeID,] ~ pheno.iop$ClozapineStatus + pheno.iop$DiseaseStatus + pheno.iop$DNAmAge + factor(pheno.iop$gender.x.chr) + pheno.iop$smokingScore + pheno.iop$CD8.naive + pheno.iop$CD8pCD28nCD45RAn + pheno.iop$PlasmaBlast + pheno.iop$CD4T + pheno.iop$NK + pheno.iop$Mono + pheno.iop$Gran + factor(pheno.iop$MethArray_Chip))
		cohortCoef[3,]<-c(summary(model.iop)$coefficients[2,c(1,2)], summary(model.iop)$coefficients[3,c(1,2)])
	} 
	if(probeID %in% rownames(betas.swedish)){
		model.swedish<-lm(betas.swedish[probeID,] ~ pheno.swedish$cloz + factor(pheno.swedish$Status, levels = c("control", "case")) + factor(pheno.swedish$sex) + pheno.swedish$ageatsampling + factor(pheno.swedish$CHIP.ID) + pheno.swedish$Gran + pheno.swedish$Mono + pheno.swedish$Bcell + pheno.swedish$NK + pheno.swedish$CD4T + pheno.swedish$CD8T + pheno.swedish$SmokingScore)
		cohortCoef[4,]<-c(summary(model.swedish)$coefficients[2,c(1,2)], summary(model.swedish)$coefficients[3,c(1,2)])
		
	}
	out.cloz<-metagen(cohortCoef[,1], cohortCoef[,2])
	out.scz<-metagen(cohortCoef[,3], cohortCoef[,4])

	tabmeta[probeID,1:21]<-c(sum(!is.na(cohortCoef[,1])), out.scz$TE.fixed,out.scz$seTE.fixed,out.scz$pval.fixed, out.scz$TE.random, out.scz$seTE.random,out.scz$pval.random, out.scz$tau, out.scz$I2, out.scz$Q,1-pchisq(out.scz$Q, out.scz$df.Q), out.cloz$TE.fixed,out.cloz$seTE.fixed,out.cloz$pval.fixed, out.cloz$TE.random, out.cloz$seTE.random,out.cloz$pval.random, out.cloz$tau, out.cloz$I2, out.cloz$Q,1-pchisq(out.cloz$Q, out.cloz$df.Q))
}

## classify results based on SCZ and TRS covariate significance
type<-rep("n.s.", nrow(tabmeta))
## sites significant both SCZ analysis and TRS analysis
type[which(tabmeta$SCZ.P_Random < 0.05 & tabmeta$Cloz.P_Random < 0.05)]<-"SCZ & TRS"
## sites only significant in TRS analysis and non-significant SCZ effect
type[which(tabmeta$SCZ.P_Random > 0.05 & tabmeta$Cloz.P_Random < 0.05)]<-"TRS only"
## sites only significant in SCZ analysis and non-significant TRS effect
type[which(tabmeta$SCZ.P_Random < 0.05 & tabmeta$Cloz.P_Random > 0.05)]<-"SCZ only"
consistentDir<-sign(tabmeta$SCZ.Effect_Random) == sign(tabmeta$Cloz.Effect_Random)

tabmeta<-cbind(tabmeta, type, consistentDir)

write.csv(tabmeta, "")
