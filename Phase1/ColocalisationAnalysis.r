### Script to perform bayesian co-localisation of mQTL and schizophreniz GWAS results
### perform analysi separately for each chromosome


library(coloc)
library(data.table)

## function to reformat SNP names from imputation
removeAlleles<-function(text){
	tmp<-unlist(strsplit(unlist(strsplit(gsub("X", "", text), "\\.")), "_"))
	if(length(tmp) == 3){
		return(tmp)
	} else {
		return(c(tmp[1:2], paste(tmp[3:length(tmp)], collapse = "_")))
	}
}


### load PGC GWAS regions
gwas<-fread("/mnt/data1/Schizophrenia/GWAS/PGC2/daner_PGC_SCZ52_0513a.hq2")
regions<-read.table("/mnt/data1/Schizophrenia/GWAS/PGC2/scz2.anneal.108.txt", header = TRUE, stringsAsFactors = FALSE)


setwd("")
meth<-read.table("", header = TRUE, row.names = 1) ## load methylation data
load("") ## load all illumina 450K probe annotation
probeAnnot<-probeAnnot[rownames(meth),] 



output<-NULL
for(chr in c(1:12, 14:20,22)){
	mQTL<-read.table(paste("", chr, "", sep = ""), stringsAsFactors = FALSE, header = TRUE) ## load mQTL results
	regions.sub<-regions[which(regions$hg19chrc == paste("chr", chr, sep = "")),] ### identify GWAS regions on chromosome in question
	locs<-matrix(unlist(lapply(mQTL$SNP, removeAlleles)), ncol = 3, byrow = TRUE) 
	mQTL<-cbind(mQTL, as.character(paste(locs[,1], locs[,2], sep = ":")))
	
	freq<-read.table(paste("", chr, ".frq", sep = ""), header = TRUE, stringsAsFactors = FALSE) ## load SNP frequency statistics
	for(i in 1:nrow(regions.sub)){ ### for each chromosome take each GWAS region
		gwas.sub<-gwas[which(gwas$CHR == chr & gwas$BP <= regions.sub$anneal2[i]+500000 & gwas$BP >= regions.sub$anneal1[i]-500000),]
		dataset2<-list(beta= log(gwas.sub$OR), varbeta= (gwas.sub$SE^2), type = "cc", snp = paste(gwas.sub$CHR, gwas.sub$BP, sep = ":"), MAF = gwas.sub$FRQ_U_46839, N=82315)
		probes<-as.character(probeAnnot[which(probeAnnot$CHR == chr & probeAnnot$MAPINFO <= regions.sub$anneal2[i]+500000 & probeAnnot$MAPINFO >= regions.sub$anneal1[i]-500000),"NAME"]) ### identify 450K probes within each region
			
		probes<-intersect(probes, rownames(meth))
		for(each in probes){ ## for each region take each 450K probe within this region and compare mQTL results to schizophrenia GWAS results 
			mQTL.sub<-mQTL[which(mQTL$gene == each),]
			if(nrow(mQTL.sub) > 1){
				dataset1<-list(beta=mQTL.sub$beta,varbeta=(mQTL.sub$beta/mQTL.sub$t.stat)^2,type = "quant",snp = mQTL.sub[,7], MAF = freq$MAF[match(mQTL.sub[,7], freq$SNP)], N=166)
				my.res<-coloc.abf(dataset1, dataset2)
				output<-rbind(output, c("region"=regions.sub$anneal.rank[i], "trait2"=each, my.res$summary))
			}
		}
	}
}


### save output
write.csv(output, "")

load("/mnt/data1/450K_reference/AllProbeIlluminaAnno.Rdata")
probeAnnot<-probeAnnot[output[,2],]
output<-cbind(output, probeAnnot[,c("CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP")])

### recreate Fig 1 from http://hmg.oxfordjournals.org/content/early/2015/04/03/hmg.ddv077.full

d0<-density(as.numeric(as.character(output[,4])), from = 0, to = 1, bw = 0.01)
d1<-density(as.numeric(as.character(output[,5]))+as.numeric(as.character(output[,6])), from = 0, to = 1, bw = 0.01)
d2<-density(as.numeric(as.character(output[,7]))+as.numeric(as.character(output[,8])), from = 0, to = 1, bw = 0.01)

pdf("MatrixEQTL/Coloc/DensityPlot_PP_ColocAnalysis.pdf")
plot(d0)
lines(d1, col = "blue")
lines(d2, col = "red")
legend("topright", lty = 1, col = c("black", "blue", "red"), c("PP0", "PP1+PP2", "PP3+PP4"))
dev.off()

### match to brain region coloc output

fetal<-read.csv("/mnt/data1/Eilis/Projects/Helen/MatrixEQTL/Output/Results/AdditiveModel/Imputed/ColocAnalysiswithin500kb.csv", stringsAsFactors = FALSE, row.names = 1)
output<-cbind(output, fetal[match(paste(output$region, output$trait2, sep = "|"), paste(fetal$region, fetal$trait2, sep = "|")),])

write.csv(output, "MatrixEQTL/Coloc/ColocAnalysiswithin500kb_withFetal.csv")

pdf("MatrixEQTL/Coloc/ComparePOsterierProbabilitiesBloodBrain.pdf")
for(i in 4:8){
	plot(output[,i], output[,(i+12)], main = colnames(output)[i], pch = 16, xlab = "Blood", ylab = "Brain")

}
dev.off()

