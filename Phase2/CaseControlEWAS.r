## This script is an example of the case control analysis run within each cohort
## File paths have been removed for security reasons

# Setting up parallel processors

library(doParallel)
cl<-makeCluster(20)
registerDoParallel(cl)

load("Normalised.rdat")

testCpG<-function(row, pheno){
	model<-lm(row ~ pheno$Status + factor(pheno$sex) + pheno$ageatsampling + factor(pheno$CHIP.ID) + pheno$Gran + pheno$Mono + pheno$Bcell + pheno$NK + pheno$CD4T + pheno$CD8T + pheno$SmokingScore)
  	return(summary(model)$coefficients["pheno$Statuscontrol",c(1,2,4)])
}



# Run EWAS using foreach() and %dopar% to run analysis is parallel
res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
		testCpG(betas[i,], pheno)
	}

colnames(res)<-c("Status_Beta", "Status_SE", "Status_P")
rownames(res)<-rownames(betas)


# Convert 0-1 proportions to percentages
res[,"Status_Beta"]<-res[,"Status_Beta"]*-100 ## flip effect sizes so controls taken as baseline
res[,"Status_SE"]<-res[,"Status_SE"]*100

# Annotate results
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7) ## manifest downloaded from illumina website
epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])


# Sorting results by P.value
res<-res[order(res$Status_P),]

# Saving results
write.csv(res, "CaseControlEWAS.csv")
