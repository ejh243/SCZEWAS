## script to test for differences in DNAm derived smoking scores between cases and controls controling for age and sex

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


## smoking analysis controlling for sex and age

t.ucl<-lm(pheno.ucl$smokingScore ~ pheno.ucl$STATUS + pheno.ucl$Sex + pheno.ucl$Age)
t.aber<-lm(pheno.aber$DNAmSmokingScore ~ pheno.aber$Status + pheno.aber$gender.x.chr + pheno.aber$Age)
t.twins<-lm(pheno.twins$scores_combined_A ~ pheno.twins$Affection + pheno.twins$Female + pheno.twins$Age)
t.iop<-lm(pheno.iop$smokingScore ~ pheno.iop$DiseaseStatus + pheno.iop$gender.x.chr + pheno.iop$AgeAtBloodCollection)
t.dublin<-lm(pheno.dublin$smokingScore ~ factor(pheno.dublin$Status, levels = c("control", "case")) + pheno.dublin$gender.x.chr + pheno.dublin$Age)
t.eugei<-lm(pheno.eugei$SmokingScore ~ factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case")) + pheno.eugei$Sex + pheno.eugei$Age)
t.swedish<-lm(pheno.swedish$SmokingScore ~ factor(pheno.swedish$Status, levels = c("control", "case")) + pheno.swedish$sex + pheno.swedish$ageatsampling)

	
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
	summary(t.swedish)$coefficients[2,"Std. Error"]))

	
out<-metagen(meanDiff, seDiff)
