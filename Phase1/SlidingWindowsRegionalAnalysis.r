### Script to perform sliding window regional analysis on results of EWAS results

library(wateRmelon)

setwd()
load() ## final QCed R object

betas<-betas(mset450k.pf.dasen) ### need beta matrix to caluclate correlations between probes

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

res<-read.csv("", row.names = 1, stringsAsFactors = FALSE) ## load EWAS results

res<-res[which(res$CHR != "Y"),] ## exclude Y chromosome probes
res<-res[which(res$SmokingPCellCompAdj > 1e-7),] ## exclude probes associated with smoking


### order by bp NOTE do not need to sort by CHR as extracting each CHR in turn
res<-res[order(res$MAPINFO),]
betas<-betas[rownames(res),]

calcCovMatrix<-function(x){
   ### x is a vector of correlations
   y<-vector(length = length(x))
   for(i in 1:length(x)){
      if(x[i] <= 1 & x[i] >= 0){
         y[i]<-x[i]*(3.25 + 0.75*x[i])
    } else {
      if(x[i] <= 0 & x[i] >= -0.5){
         y[i]<-x[i]*(3.27 + 0.71*x[i])
      }
    }
   
  }
  return(y)
}

brownsP<-function(covar, pval){
	## covar is vector of covariances between all pairs of tests NOTE not correlations
	## pval is vector of p values
	ntests<-length(pval)
	var_xsq<-sum(covar)*2 + 4*ntests 	# equation (3)
	exp_xsq<-2*ntests	#equation (2)

	## estimate parameters for chi square distribution using Brown's notation
	f = (2*(exp_xsq^2))/var_xsq
	c = var_xsq/(2*exp_xsq)

	##### NOTE: doesn't match Brown's answer but matches my answer by hand
	chi_sq<- -2*sum(log(pval))

	### to obtain p value
	test.stat<-chi_sq/c
	browns.pval<-1-pchisq(test.stat, df = f)
	return(c(browns.pval, test.stat, f))
}


windows<-c(100,200,500,1000,2000,5000)
out<-NULL
for(chr in 1:22){
	print(paste("CHR: ", chr))
	res.tmp<-res[which(res$CHR == chr),]
	
	out.tmp<-matrix(data = NA, nrow = nrow(res.tmp), ncol = 2+length(windows)*5)
	colnames(out.tmp)<-c("Central Probe", "Chr", rep(c("Start Region", "End Region", "nProbes", "Gene Anno", "Browns P"), length(windows)))
	out.tmp[,1]<-rownames(res.tmp)	
	out.tmp[,2]<-chr
	for(i in 1:nrow(res.tmp)){

		chr<-res.tmp$CHR[i]
		bp<-res.tmp$MAPINFO[i]

		for(j in 1:length(windows)){
		sub<-res.tmp[which(res.tmp$CHR == chr & res.tmp$MAPINFO <= (bp+windows[j]) & res.tmp$MAPINFO >= (bp-windows[j])),]
	
			out.tmp[i,(j*6)-3]<-min(sub$MAPINFO)
			out.tmp[i,(j*6)-2]<-max(sub$MAPINFO)
			out.tmp[i,(j*6)-1]<-nrow(sub)
			
			if(nrow(sub) > 1){
				out.tmp[i,(j*6)]<-paste(unique(unlist(strsplit(sub$UCSC_REFGENE_NAME, ";"))), collapse = ";")
				betas.sub<-betas[rownames(sub),]
				
				cor.mat<-c(unique(as.vector(cor(t(betas.sub)))))[-1]
				covar<-calcCovMatrix(cor.mat)

				pval<-sub$P.value

				brownsp<-brownsP(covar, pval)
				out.tmp[i,(j*6)+1]<-brownsp[1]
			}
			
		}
	}
out<-rbind(out,out.tmp)
	
}


### save results
write.csv(out, "")

