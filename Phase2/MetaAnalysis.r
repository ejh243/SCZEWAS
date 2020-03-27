## Meta-analysis of EWAS results using mean difference and standard error
## this script is the basis for 1) psychosis vs controls; 2) schizophrenia vs controls; 3) treatment resistent schizophrenia cases vs non-treatment resistent schizophrenia cases
## only for sites tested in at least 2 cohorts 

library(bacon)
library(meta)
library(BiocParallel)

setwd("")


## load EWAS results and filter out Y chromosome
res.ucl<-read.csv("", row.names = 1, stringsAsFactors = FALSE)
res.ucl<-res.ucl[which(res.ucl$CHR != "Y"),]

res.aber<-read.csv("", stringsAsFactors = FALSE,row.names = 1,)
res.aber<-res.aber[which(res.aber$CHR != "Y"),]

res.iop<-read.csv("", row.names = 1, stringsAsFactors = FALSE)
res.iop<-res.iop[which(res.iop$CHR != "Y"),]

res.dublin<-read.csv("", stringsAsFactors = FALSE,row.names = 1)
res.dublin<-res.dublin[which(res.dublin$CHR != "Y"),]

res.eugei<-read.csv("", stringsAsFactors = FALSE,row.names = 1)
res.eugei<-res.eugei[which(res.eugei$CHR != "Y"),]

res.twins<-read.csv("", row.names = 1, stringsAsFactors = FALSE)
res.twins<-res.twins[which(res.twins$CHR != "Y"),]

res.swedish<-read.csv", row.names = 1, stringsAsFactors = FALSE)
res.swedish<-res.swedish[which(res.swedish$CHR != "Y"),]


## need to convert mean difference from some cohorts so all relative to controls
## multiply by 100 to convert from proportion to percentage
res.ucl$Beta<-res.ucl$Beta*100
res.aber$Beta<-res.aber$Beta*100
res.twins$ClusteredRobust.Beta<-res.twins$ClusteredRobust.Beta*100
res.dublin$Beta<-res.dublin$Beta*-100
res.iop$Beta<-res.iop$Beta*100

res.ucl$SE<-res.ucl$SE*100
res.aber$SE<-res.aber$SE*100
res.twins$ClusteredRobust.SE<-res.twins$ClusteredRobust.SE*100
res.dublin$SE<-res.dublin$SE*100
res.iop$SE<-res.iop$SE*100

## apply bacon correction
bc.ucl<-bacon(NULL, res.ucl$Beta, res.ucl$SE)
bc.aber<-bacon(NULL, res.aber$Beta, res.aber$SE)
bc.twins<-bacon(NULL, res.twins$ClusteredRobust.Beta, res.twins$ClusteredRobust.SE)
bc.dublin<-bacon(NULL, res.dublin$Beta, res.dublin$SE)
bc.swedish<-bacon(NULL, res.swedish$Status_Beta, res.swedish$Status_SE)
bc.iop<-bacon(NULL, res.iop$Beta, res.iop$SE)
bc.eugei<-bacon(NULL, res.eugei$Beta, res.eugei$SE)

## identify probes that were tested in at least 2 cohorts 
probes<-table(c(rownames(res.ucl), rownames(res.aber), rownames(res.twins), rownames(res.iop), rownames(res.dublin),  rownames(res.swedish),  rownames(res.eugei)))
probes<-probes[which(probes > 1)]


## create matrix to store results
res.meta<-matrix(data = NA, nrow = length(probes), ncol = 11)
rownames(res.meta)<-names(probes)
colnames(res.meta)<-c("N_cohorts", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")

## run meta-analysis
for(i in 1:length(probes)){
	probeID<-names(probes)[i]
	## first run meta analysis with dublin included
	meanDiff<-c(res.ucl$Beta[match(probeID, rownames(res.ucl))], 
		res.aber$Beta[match(probeID, rownames(res.aber))], 
		res.twins$ClusteredRobust.Beta[match(probeID, rownames(res.twins))],  
		res.iop$Beta[match(probeID, rownames(res.iop))], 
		res.dublin$Beta[match(probeID, rownames(res.dublin))],  
		res.swedish$Status_Beta[match(probeID, rownames(res.swedish))],  
		res.eugei$Beta[match(probeID, rownames(res.eugei))])
	
	seDiff<-c(res.ucl$SE[match(probeID, rownames(res.ucl))], 
		res.aber$SE[match(probeID, rownames(res.aber))], 	
		res.twins$ClusteredRobust.SE[match(probeID, rownames(res.twins))], 
		res.iop$SE[match(probeID, rownames(res.iop))], 
		res.dublin$SE[match(probeID, rownames(res.dublin))], 
		res.swedish$Status_SE[match(probeID, rownames(res.swedish))], 
		res.eugei$SE[match(probeID, rownames(res.eugei))])
	
	out<-metagen(meanDiff, seDiff)
	res.meta[i,1:11]<-c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q))
	}

## add in results from cohorts
res.meta<-cbind(res.meta, 
	res.ucl[match(names(probes), rownames(res.ucl)), c("Beta", "SE")], 
	res.aber[match(names(probes), rownames(res.aber)), c("Beta", "SE")], 
	res.twins[match(names(probes), rownames(res.twins)), c("ClusteredRobust.Beta", "ClusteredRobust.SE")], 
	res.iop[match(names(probes), rownames(res.iop)), c("Beta", "SE")], 
	res.dublin[match(names(probes), rownames(res.dublin)), c("Beta", "SE")], 
	res.swedish[match(names(probes), rownames(res.swedish)),c("Status_Beta", "Status_SE")], 
	res.swedish[match(names(probes), rownames(res.eugei)),c("Beta", "SE")])

colnames(res.meta)[-c(1:11)]<-outer(c("Beta_", "SE_"),c("UCL", "Aberdeen", "Twins", "IoPPN", "Dublin", "Sweden", "EUGEI"),  FUN = "paste0")
	
write.csv(res.meta, "")

