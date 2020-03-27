## script to test for differences in cellular composition between cases and controls controlling for smoking, age and sex

library(meta)
library(grid)
library(wateRmelon)


setwd("")

## load data for each cohort and resave pheno dataframe as separate R object called pheno.<cohortName>
setwd("UCL")
load("Normalised.rdat")
pheno.ucl<-pheno

setwd("Aberdeen")
load("Normalised.rdat")
pheno.aber<-pheno

setwd("Twins")
load("Normalised.rdat")
pheno.twins<-pheno

setwd("IOP")
load("Normalised.rdat")
pheno.iop<-pheno

setwd("Dublin")
load("Normalised.rdat")
pheno.dublin<-pheno

setwd("EUGEI")
load("Normalised.rdat")
pheno.eugei<-pheno

setwd("Sweden")
load("Normalised.rdat")
pheno.swedish<-pheno

## for all summary tables order as follows: UCL, Aberdeen, Twins, IOP, Dublin, EUGEI, Swedish
cohortNames<-c("UCL", "Aberdeen", "Twins", "IoPPN", "Dublin", "EUGEI", "Sweden")
cellTypes<-c("Mono","Gran","NK","CD4T","CD8T","Bcell","PlasmaBlast", "CD8pCD28nCD45RAn","CD8.naive","CD4.naive")

## meta analysis across cohorts
tab.cellcomp.meta<-matrix(data = NA, nrow = 10, ncol = 11)
rownames(tab.cellcomp.meta)<-cellTypes
colnames(tab.cellcomp.meta)<-c("N_cohorts", "Effect_Fixed", "Effect_SE_Fixed", "P_Fixed", "Effect_Random", "Effect_SE_Random","P_Random", "tau", "I2", "Q", "Het P")

pdf("", height = 8, width = 12)
for(each in cellTypes[c(1:6)]){
	t.ucl<-lm(pheno.ucl[,each] ~ pheno.ucl$STATUS + pheno.ucl$smokingScore + pheno.ucl$Age + pheno.ucl$Sex)
	t.aber<-lm(pheno.aber[,each] ~ pheno.aber$Status + pheno.aber$DNAmSmokingScore + pheno.aber$Age + pheno.aber$gender.x.chr)
	t.twins<-lm(pheno.twins[,each] ~ pheno.twins$Affection + pheno.twins$scores_combined_A + pheno.twins$Age + pheno.twins$Female)
	t.iop<-lm(pheno.iop[,each] ~ pheno.iop$DiseaseStatus + pheno.iop$smokingScore + pheno.iop$AgeAtBloodCollection + pheno.iop$gender.x.chr)
	t.dublin<-lm(pheno.dublin[,each] ~ factor(pheno.dublin$Status, levels = c("control", "case")) + pheno.dublin$smokingScore + pheno.dublin$Age + pheno.dublin$gender.x.chr)
	t.eugei<-lm(pheno.eugei[,each] ~ factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case")) + pheno.eugei$SmokingScore + pheno.eugei$Age + pheno.eugei$Sex)
	t.swedish<-lm(pheno.swedish[,each] ~ factor(pheno.swedish$Status, levels = c("control", "case")) + pheno.swedish$SmokingScore + pheno.swedish$ageatsampling + pheno.swedish$sex)
	
	meanDiff<-c(summary(t.ucl)$coefficients[2,"Estimate"],
	summary(t.aber)$coefficients[2,"Estimate"],
	summary(t.twins)$coefficients[2,"Estimate"],
	summary(t.iop)$coefficients[2,"Estimate"],
	summary(t.dublin)$coefficients[2,"Estimate"],
	summary(t.eugei)$coefficients[2,"Estimate"],
	summary(t.swedish)$coefficients[2,"Estimate"])
	seDiff<-abs(c(summary(t.ucl)$coefficients[2,"Std. Error"],
	summary(t.aber)$coefficients[2,"Std. Error"],
	summary(t.twins)$coefficients[2,"Std. Error"],
	summary(t.iop)$coefficients[2,"Std. Error"],
	summary(t.dublin)$coefficients[2,"Std. Error"],
	summary(t.eugei)$coefficients[2,"Std. Error"],
	summary(t.swedish)$coefficients[2,"Std. Error"])
	)
	tab.cellcomp<-rbind(tab.cellcomp,c(summary(t.ucl)$coefficients[2,"Pr(>|t|)"],
	summary(t.aber)$coefficients[2,"Pr(>|t|)"],
	summary(t.twins)$coefficients[2,"Pr(>|t|)"],
	summary(t.iop)$coefficients[2,"Pr(>|t|)"],
	summary(t.dublin)$coefficients[2,"Pr(>|t|)"],
	summary(t.eugei)$coefficients[2,"Pr(>|t|)"],
	summary(t.swedish)$coefficients[2,"Pr(>|t|)"]))
	tab.cellcomp<-rbind(tab.cellcomp,meanDiff/seDiff)
	
	out<-metagen(meanDiff, seDiff)
	tab.cellcomp.meta[each,1:11]<-c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q))
	forest(out, col.square = "blue", col.diamond = "forestgreen",test.overall.fixed = FALSE, print.tau2 = FALSE,  studlab = c("UCL", "Aberdeen", "Twins", "IoPPN", "Dublin", "EUGEI", "Sweden"), xlab = "Mean difference", digits = 4)
grid.text(each, .5, .85, gp=gpar(cex=1.5))
}

for(each in cellTypes[c(7:10)]){

	t.ucl<-lm(pheno.ucl[,each] ~ pheno.ucl$STATUS + pheno.ucl$smokingScore + pheno.ucl$Age + pheno.ucl$Sex)
	t.aber<-lm(pheno.aber[,each] ~ pheno.aber$Status + pheno.aber$DNAmSmokingScore + pheno.aber$Age + pheno.aber$gender.x.chr)
	t.twins<-lm(pheno.twins[,each] ~ pheno.twins$Affection + pheno.twins$scores_combined_A + pheno.twins$Age + pheno.twins$Female)
	t.iop<-lm(pheno.iop[,each] ~ pheno.iop$DiseaseStatus + pheno.iop$smokingScore + pheno.iop$AgeAtBloodCollection + pheno.iop$gender.x.chr)
	t.dublin<-lm(pheno.dublin[,each] ~ factor(pheno.dublin$Status, levels = c("control", "case")) + pheno.dublin$smokingScore + pheno.dublin$Age + pheno.dublin$gender.x.chr)
	
	meanDiff<-c(summary(t.ucl)$coefficients[2,"Estimate"],
	summary(t.aber)$coefficients[2,"Estimate"],
	summary(t.twins)$coefficients[2,"Estimate"],
	summary(t.iop)$coefficients[2,"Estimate"],
	summary(t.dublin)$coefficients[2,"Estimate"])
	seDiff<-abs(c(summary(t.ucl)$coefficients[2,"Std. Error"],
	summary(t.aber)$coefficients[2,"Std. Error"],
	summary(t.twins)$coefficients[2,"Std. Error"],
	summary(t.iop)$coefficients[2,"Std. Error"],
	summary(t.dublin)$coefficients[2,"Std. Error"])
	)	
	out<-metagen(meanDiff, seDiff)
	tab.cellcomp.meta[each,1:11]<-c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q))
	forest(out, col.square = "blue", col.diamond = "forestgreen",test.overall.fixed = FALSE, print.tau2 = FALSE,  studlab = c("UCL", "Aberdeen", "Twins", "IoPPN", "Dublin"), xlab = "Mean difference", digits = 4)
	grid.text(each, .5, .85, gp=gpar(cex=1.5))
	
	tab.cellcomp<-rbind(tab.cellcomp,c(summary(t.ucl)$coefficients[2,"Pr(>|t|)"],
	summary(t.aber)$coefficients[2,"Pr(>|t|)"],
	summary(t.twins)$coefficients[2,"Pr(>|t|)"],
	summary(t.iop)$coefficients[2,"Pr(>|t|)"],
	summary(t.dublin)$coefficients[2,"Pr(>|t|)"]))
	tab.cellcomp<-rbind(tab.cellcomp,meanDiff/seDiff)
}

write.csv(tab.cellcomp.meta, "")
dev.off()

