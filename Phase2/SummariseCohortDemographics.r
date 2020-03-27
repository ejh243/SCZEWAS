## script to generate summary tables, plot and test for biases in demogrpahics across cohorts included in meta-analysis

library(meta)
library(grid)
library(wateRmelon)
library(WGCNA)

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

### tabulate sex

tab.sex<-cbind(c(table(pheno.ucl$Sex), table(pheno.ucl$Sex)/nrow(pheno.ucl)*100),
c(table(pheno.aber$gender.x.chr), table(pheno.aber$gender.x.chr)/nrow(pheno.aber)*100),
c(table(1-pheno.twins$Female), table(1-pheno.twins$Female)/nrow(pheno.twins)*100), ## need to flip so females reported first
c(table(pheno.iop$gender.x.chr), table(pheno.iop$gender.x.chr)/nrow(pheno.iop)*100),
c(table(pheno.dublin$gender.x.chr), table(pheno.dublin$gender.x.chr)/nrow(pheno.dublin)*100),
c(table(pheno.eugei$Sex), table(pheno.eugei$Sex)/nrow(pheno.eugei)*100),
c(table(pheno.swedish$sex), table(pheno.swedish$sex)/nrow(pheno.swedish)*100))
colnames(tab.sex)<-cohortNames
rownames(tab.sex)<-c("Females", "Males", "%Females", "%Males")
write.csv(tab.sex, "")

### tabulate gender by case control status

tab.sexbystatus<-cbind(c(table(pheno.ucl$Sex, pheno.ucl$STATUS)),
c(table(pheno.aber$gender.x.chr, pheno.aber$Status)),
c(table((1-pheno.twins$Female), pheno.twins$Affection)),
c(table(pheno.iop$gender.x.chr, pheno.iop$DiseaseStatus)),
c(table(pheno.dublin$gender.x.chr, factor(pheno.dublin$Status, levels = c("control", "case")))),
c(table(pheno.eugei$Sex, factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case")))),
c(table(pheno.swedish$sex, factor(pheno.swedish$Status, levels = c("control", "case")))))

tab.sexbystatus<-rbind(tab.sexbystatus, c(chisq.test(table(pheno.ucl$Sex, pheno.ucl$STATUS))$p.value,
chisq.test(table(pheno.aber$gender.x.chr, pheno.aber$Status))$p.value,
chisq.test(table(pheno.twins$Female, pheno.twins$Affection))$p.value,
chisq.test(table(pheno.iop$gender.x.chr, pheno.iop$DiseaseStatus))$p.value,
chisq.test(table(pheno.dublin$gender.x.chr, pheno.dublin$Status))$p.value,
chisq.test(table(pheno.eugei$Sex, pheno.eugei$Phenotype.EuGEI))$p.value,
chisq.test(table(pheno.swedish$sex, pheno.swedish$Status))$p.value))
rownames(tab.sexbystatus)<-c("Females_controls", "Males_controls", "Females_cases", "Males_cases", "ChiSqTest")

## add total column
tab.sexbystatus<-cbind(tab.sexbystatus, NA)
colnames(tab.sexbystatus)<-c(cohortNames, "Combined")
tab.sexbystatus[,"Combined"]<-rowSums(tab.sexbystatus,na.rm = TRUE)
tab.sexbystatus[5,"Combined"]<-chisq.test(matrix(tab.sexbystatus[1:4,"Combined"], nrow = 2))$p.value

write.csv(tab.sexbystatus, "")
chisq.test(tab.sexbystatus[1:4,])

## summarise age

tab.age<-rbind(c(mean(pheno.ucl$Age, na.rm = TRUE),
mean(pheno.aber$Age, na.rm = TRUE),
mean(pheno.twins$Age, na.rm = TRUE),
mean(pheno.iop$AgeAtBloodCollection, na.rm = TRUE),
mean(pheno.dublin$Age_at_Collection, na.rm = TRUE),
mean(pheno.eugei$Age, na.rm = TRUE),
mean(pheno.swedish$ageatsampling, na.rm = TRUE)),

c(sd(pheno.ucl$Age, na.rm = TRUE),
sd(pheno.aber$Age, na.rm = TRUE),
sd(pheno.twins$Age, na.rm = TRUE),
sd(pheno.iop$AgeAtBloodCollection, na.rm = TRUE),
sd(pheno.dublin$Age_at_Collection, na.rm = TRUE),
sd(pheno.eugei$Age, na.rm = TRUE),
sd(pheno.swedish$ageatsampling, na.rm = TRUE)))

tab.age<-rbind(tab.age, cbind(aggregate(pheno.ucl$Age, by = list(pheno.ucl$STATUS),mean, na.rm = TRUE)[,2],
aggregate(pheno.aber$Age, by = list(pheno.aber$Status),mean, na.rm = TRUE)[,2],
aggregate(pheno.twins$Age, by = list(pheno.twins$Affection),mean, na.rm = TRUE)[,2],
aggregate(pheno.iop$AgeAtBloodCollection, by = list(pheno.iop$DiseaseStatus),mean, na.rm = TRUE)[,2],
aggregate(pheno.dublin$Age_at_Collection, by = list(factor(pheno.dublin$Status, levels = c("control", "case"))),mean, na.rm = TRUE)[,2],
aggregate(pheno.eugei$Age, by = list(factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case"))),mean, na.rm = TRUE)[,2],
aggregate(pheno.swedish$ageatsampling, by = list(factor(pheno.swedish$Status, levels = c("control", "case"))),mean, na.rm = TRUE)[,2]))

tab.age<-rbind(tab.age, c(t.test(pheno.ucl$Age ~ pheno.ucl$STATUS)$p.value,
t.test(pheno.aber$Age ~ pheno.aber$Status)$p.value,
t.test(pheno.twins$Age ~ pheno.twins$Affection)$p.value,
t.test(pheno.iop$AgeAtBloodCollection ~ pheno.iop$DiseaseStatus)$p.value,
t.test(pheno.dublin$Age_at_Collection ~ pheno.dublin$Status)$p.value,
t.test(pheno.eugei$Age ~ pheno.eugei$Phenotype.EuGEI)$p.value,
t.test(pheno.swedish$ageatsampling ~ factor(pheno.swedish$Status, levels = c("control", "case")))$p.value))
tab.age<-cbind(tab.age, NA)
colnames(tab.age)<-c(cohortNames, "Combined")
rownames(tab.age)<-c("MeanAge", "SDAge", "MeanAge_controls", "MeanAge_cases", "T-Test")


## meta-analysis of age

t.ucl<-t.test(pheno.ucl$Age ~ pheno.ucl$STATUS)
t.aber<-t.test(pheno.aber$Age ~ pheno.aber$Status)
t.twins<-t.test(pheno.twins$Age ~ pheno.twins$Affection)
t.iop<-t.test(pheno.iop$AgeAtBloodCollection ~ pheno.iop$DiseaseStatus)
t.dublin<-t.test(pheno.dublin$Age_at_Collection ~ factor(pheno.dublin$Status, levels = c("control", "case")))
t.eugei<-t.test(pheno.eugei$Age ~ factor(pheno.eugei$Phenotype.EuGEI, levels = c("Control", "Case")))
t.swedish<-t.test(pheno.swedish$ageatsampling ~ factor(pheno.swedish$Status, levels = c("control", "case")))
	
meanDiff<-c(diff(t.ucl$estimate),
	diff(t.aber$estimate),
	diff(t.twins$estimate),
	diff(t.iop$estimate),
	diff(t.dublin$estimate),
	diff(t.eugei$estimate),
	diff(t.swedish$estimate))
seDiff<-abs(c(diff(t.ucl$estimate)/t.ucl$statistic,
	diff(t.aber$estimate)/t.aber$statistic,
	diff(t.twins$estimate)/t.twins$statistic,
	diff(t.iop$estimate)/t.iop$statistic,
	diff(t.dublin$estimate)/t.dublin$statistic,
	diff(t.eugei$estimate)/t.eugei$statistic,
	diff(t.swedish$estimate)/t.swedish$statistic
	))
	
out<-metagen(meanDiff, seDiff)
tab.age[5,"Combined"]<-out$pval.random

## calc pooled means
tab.age[1,"Combined"]<-mean(c(pheno.ucl$Age,pheno.aber$Age, pheno.twins$Age, pheno.iop$AgeAtBloodCollection, pheno.dublin$Age_at_Collection, pheno.eugei$Age, pheno.swedish$ageatsampling), na.rm = TRUE)
tab.age[2,"Combined"]<-sd(c(pheno.ucl$Age,pheno.aber$Age, pheno.twins$Age, pheno.iop$AgeAtBloodCollection, pheno.dublin$Age_at_Collection, pheno.eugei$Age, pheno.swedish$ageatsampling), na.rm = TRUE)

totAge<-cbind(aggregate(pheno.ucl$Age, by = list(pheno.ucl$STATUS),sum, na.rm = TRUE)[,2],
aggregate(pheno.aber$Age, by = list(pheno.aber$Status),sum, na.rm = TRUE)[,2],
aggregate(pheno.twins$Age, by = list(pheno.twins$Affection),sum, na.rm = TRUE)[,2],
aggregate(pheno.iop$AgeAtBloodCollection, by = list(pheno.iop$DiseaseStatus),sum, na.rm = TRUE)[,2],
aggregate(pheno.dublin$Age_at_Collection, by = list(factor(pheno.dublin$Status, levels = c("control", "case"))),sum, na.rm = TRUE)[,2],
aggregate(pheno.eugei$Age, by = list(factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case"))),sum, na.rm = TRUE)[,2],
aggregate(pheno.swedish$ageatsampling, by = list(factor(pheno.swedish$Status, levels = c("control", "case"))),sum, na.rm = TRUE)[,2])

totN<-cbind(aggregate(!is.na(pheno.ucl$Age), by = list(pheno.ucl$STATUS),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.aber$Age), by = list(pheno.aber$Status),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.twins$Age), by = list(pheno.twins$Affection),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.iop$AgeAtBloodCollection), by = list(pheno.iop$DiseaseStatus),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.dublin$Age_at_Collection), by = list(factor(pheno.dublin$Status, levels = c("control", "case"))),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.eugei$Age), by = list(factor(pheno.eugei$Phenotype.EuGEI,levels = c("Control", "Case"))),sum, na.rm = TRUE)[,2],
aggregate(!is.na(pheno.swedish$ageatsampling), by = list(factor(pheno.swedish$Status, levels = c("control", "case"))),sum, na.rm = TRUE)[,2])


tab.age[3,"Combined"]<-sum(totAge[1,])/sum(totN[1,])
tab.age[4,"Combined"]<-sum(totAge[2,])/sum(totN[2,])

write.csv(tab.age, "")

### summarise DNAmAge against age 

pdf("", height = 8, width = 14)
par(mfrow = c(2,4))
par(mar = c(5,5,3,0.2))
xlim<-range(c(pheno.ucl$Age, pheno.aber$Age, pheno.iop$AgeAtBloodCollection, pheno.dublin$Age_at_Collection, pheno.twins$Age, pheno.eugei$Age, pheno.swedish$ageatsampling), na.rm = TRUE)
ylim<-range(c(pheno.ucl$DNAmAge, pheno.aber$DNAmAge, pheno.iop$DNAmAge, pheno.dublin$DNAmAge, pheno.twins$DNAmAge, pheno.eugei$PredictedAge, pheno.swedish$PredictedAge), na.rm = TRUE)

plot(pheno.ucl$Age, pheno.ucl$DNAmAge, col = c("red", "blue")[as.factor(pheno.ucl$STATUS)], , main = "UCL", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
abline(a = 0, b = 1)
model<-lm(pheno.ucl$DNAmAge ~ pheno.ucl$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.ucl$Age, pheno.ucl$DNAmAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.ucl$DNAmAge ~ pheno.ucl$Age*as.factor(pheno.ucl$STATUS))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)
legend("topleft", c("SCZ", "CON"), col = c("blue", "red"), pch = 16) 

plot(pheno.aber$Age,pheno.aber$DNAmAge,col = c("red", "blue")[as.factor(pheno.aber$Status)], main = "Aberdeen", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
abline(a = 0, b = 1)
model<-lm(pheno.aber$DNAmAge ~ pheno.aber$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.aber$Age,pheno.aber$DNAmAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.aber$DNAmAge ~ pheno.aber$Age*as.factor(pheno.aber$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.twins$Age, pheno.twins$DNAmAge,col = c("red", "blue")[as.factor(pheno.twins$Affection)], main = "Twins", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
abline(a = 0, b = 1)
model<-lm(pheno.twins$DNAmAge ~ pheno.twins$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.twins$Age, pheno.twins$DNAmAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.twins$DNAmAge ~ pheno.twins$Age*as.factor(pheno.twins$Affection))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.iop$AgeAtBloodCollection, pheno.iop$DNAmAge,col = c("red", "blue")[as.factor(pheno.iop$DiseaseStatus)], main = "IoPPN", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
abline(a = 0, b = 1)
model<-lm(pheno.iop$DNAmAge ~pheno.iop$AgeAtBloodCollection)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.iop$AgeAtBloodCollection, pheno.iop$DNAmAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.iop$DNAmAge ~ pheno.iop$AgeAtBloodCollection*as.factor(pheno.iop$DiseaseStatus))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.dublin$Age_at_Collection, pheno.dublin$DNAmAge,col = c("blue","red")[as.factor(pheno.dublin$Status)], main = "Dublin", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
model<-lm(pheno.dublin$DNAmAge ~ pheno.dublin$Age_at_Collection)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.dublin$Age_at_Collection, pheno.dublin$DNAmAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.dublin$DNAmAge ~ pheno.dublin$Age_at_Collection*as.factor(pheno.dublin$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.eugei$Age, pheno.eugei$PredictedAge,col = c("blue","red")[as.factor(pheno.eugei$Phenotype.EuGEI)], main = "EUGEI", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
model<-lm(pheno.eugei$PredictedAge ~ pheno.eugei$Age)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.eugei$Age, pheno.eugei$PredictedAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.eugei$PredictedAge ~ pheno.eugei$Age*as.factor(pheno.eugei$Phenotype.EuGEI))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.swedish$ageatsampling, pheno.swedish$PredictedAge,col = c("blue","red")[as.factor(pheno.swedish$Status)], main = "Sweden", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "DNAmAge")
model<-lm(pheno.swedish$PredictedAge ~ pheno.swedish$ageatsampling)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.swedish$ageatsampling, pheno.swedish$PredictedAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.swedish$PredictedAge ~ pheno.swedish$ageatsampling*as.factor(pheno.swedish$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)
dev.off()

## try summarise PhenoAge against age

pdf("", height = 8, width = 14)
par(mfrow = c(2,4))
par(mar = c(5,5,3,0.2))
xlim<-range(c(pheno.ucl$Age, pheno.aber$Age, pheno.iop$AgeAtBloodCollection, pheno.dublin$Age_at_Collection, pheno.twins$Age, pheno.eugei$Age, pheno.swedish$ageatsampling), na.rm = TRUE)
ylim<-range(c(pheno.ucl$PhenoAge, pheno.aber$PhenoAge, pheno.iop$PhenoAge, pheno.dublin$PhenoAge, pheno.twins$PhenoAge, pheno.eugei$PhenoAge, pheno.swedish$PhenoAge), na.rm = TRUE)

plot(pheno.ucl$Age, pheno.ucl$PhenoAge, col = c("red", "blue")[as.factor(pheno.ucl$STATUS)], , main = "UCL", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
abline(a = 0, b = 1)
model<-lm(pheno.ucl$PhenoAge ~ pheno.ucl$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.ucl$Age, pheno.ucl$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.ucl$PhenoAge ~ pheno.ucl$Age*as.factor(pheno.ucl$STATUS))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)
legend("topleft", c("SCZ", "CON"), col = c("blue", "red"), pch = 16) 

plot(pheno.aber$Age,pheno.aber$PhenoAge,col = c("red", "blue")[as.factor(pheno.aber$Status)], main = "Aberdeen", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
abline(a = 0, b = 1)
model<-lm(pheno.aber$PhenoAge ~ pheno.aber$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.aber$Age,pheno.aber$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.aber$PhenoAge ~ pheno.aber$Age*as.factor(pheno.aber$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.twins$Age, pheno.twins$PhenoAge,col = c("red", "blue")[as.factor(pheno.twins$Affection)], main = "Twins", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
abline(a = 0, b = 1)
model<-lm(pheno.twins$PhenoAge ~ pheno.twins$Age)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.twins$Age, pheno.twins$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.twins$PhenoAge ~ pheno.twins$Age*as.factor(pheno.twins$Affection))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.iop$AgeAtBloodCollection, pheno.iop$PhenoAge,col = c("red", "blue")[as.factor(pheno.iop$DiseaseStatus)], main = "IoPPN", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
abline(a = 0, b = 1)
model<-lm(pheno.iop$PhenoAge ~pheno.iop$AgeAtBloodCollection)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.iop$AgeAtBloodCollection, pheno.iop$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.iop$PhenoAge ~ pheno.iop$AgeAtBloodCollection*as.factor(pheno.iop$DiseaseStatus))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.dublin$Age_at_Collection, pheno.dublin$PhenoAge,col = c("blue","red")[as.factor(pheno.dublin$Status)], main = "Dublin", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
model<-lm(pheno.dublin$PhenoAge ~ pheno.dublin$Age_at_Collection)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.dublin$Age_at_Collection, pheno.dublin$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.dublin$PhenoAge ~ pheno.dublin$Age_at_Collection*as.factor(pheno.dublin$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)

plot(pheno.eugei$Age, pheno.eugei$PhenoAge,col = c("blue","red")[as.factor(pheno.eugei$Phenotype.EuGEI)], main = "EUGEI", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
model<-lm(pheno.eugei$PhenoAge ~ pheno.eugei$Age)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.eugei$Age, pheno.eugei$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.eugei$PhenoAge ~ pheno.eugei$Age*as.factor(pheno.eugei$Phenotype.EuGEI))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)


plot(pheno.swedish$ageatsampling, pheno.swedish$PhenoAge,col = c("blue","red")[as.factor(pheno.swedish$Status)], main = "Sweden", cex.main = 2, cex.axis = 1.75, cex.lab = 1.75, xlim = xlim, ylim = ylim, pch = 16, xlab = "Age", ylab = "PhenoAge")
model<-lm(pheno.swedish$PhenoAge ~ pheno.swedish$ageatsampling)
abline(a = 0, b = 1)
abline(model, lty = 2)
mtext(paste("r =", signif(cor(pheno.swedish$ageatsampling, pheno.swedish$PhenoAge, use = "pairwise.complete.obs"),3)), side = 3, line = 0.5, adj = 1)
model<-lm(pheno.swedish$PhenoAge ~ pheno.swedish$ageatsampling*as.factor(pheno.swedish$Status))
text(paste("Int P =", signif(summary(model)$coefficients[4,4],3)),x = xlim[2]*0.7, y = ylim[1], adj = 0, cex = 1.45)
dev.off()

