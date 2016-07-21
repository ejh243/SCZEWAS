#### calculate Brown's combined P for each SCZ GWAS region using EWAS results

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


library(data.table)
gwas<-fread("")
scz.anno<-read.table("", stringsAsFactors = FALSE, header = TRUE) ### load clumped GWAS results i.e. list of regions


setwd("")
res.ewas<-read.csv("", row.names = 1, stringsAsFactors = FALSE) ## load case control EWAS results
res.poly<-read.csv("", stringsAsFactors = FALSE, row.names = 1) ## load PRS EWAS results


res.ewas<-res.ewas[which(res.ewas$CHR != "X"),] ## filter to autosomes
res.ewas<-res.ewas[which(res.ewas$CHR != "Y"),] ## filter to autosomes
res.ewas<-res.ewas[-which(res.ewas$SmokingPCellCompAdj < 1e-7),] ## remove smoking associated probes


res.poly<-res.poly[rownames(res.ewas),]
res.poly<-res.poly[rownames(res.ewas),]


### for each region calculate BrownsP
library(wateRmelon)

load("") ## load QC methylation data

betas<-betas(mset450k.pf.dasen) ## extract betas matrix to calculate probe correlations

pheno<-pheno[match(colnames(betas), pheno$Basename),]
betas<-betas[rownames(res.ewas),]


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

out<-matrix(data = NA, nrow = nrow(scz.anno), ncol = 6)
colnames(out)<-c("nProbe", "minCaseConP", "Browns_CaseConP", "minPolygenicP", "Browns_PolygenicP", "MeanCor")

for(i in 1:nrow(scz.anno)){
	## identify probes within each region
	res.sub<-res.ewas[which(res.ewas$CHR == gsub("chr", "", scz.anno$hg19chrc[i]) & res.ewas$MAPINFO >= scz.anno$anneal1[i] & res.ewas$MAPINFO <= scz.anno$anneal2[i]),]
	res.poly.sub<-res.poly[which(res.poly$CHR == gsub("chr", "", scz.anno$hg19chrc[i]) & res.poly$MAPINFO >= scz.anno$anneal1[i] & res.poly$MAPINFO <= scz.anno$anneal2[i]),]
	
	out[i,1]<-nrow(res.sub)
	if(nrow(res.sub) > 1){ ### if region has more than 1 probe we can calculate a combined p value
		betas.sub<-betas[rownames(res.sub),]
		
		cor.mat<-c(unique(as.vector(cor(t(betas.sub)))))[-1]
		out[i,6]<-mean(cor.mat)
		covar<-calcCovMatrix(cor.mat)

		pval<-res.sub$P.value
		out[i,2]<-min(res.sub$P.value)

		brownsp<-brownsP(covar, pval)
		out[i,3]<-brownsp[1]
		
		pval<-res.poly.sub$P.value
		out[i,4]<-min(res.poly.sub$P.value)

		brownsp<-brownsP(covar, pval)
		out[i,5]<-brownsp[1]		
	} 
}

out<-cbind(scz.anno, out) ## merge results with GWAS data for each region

### save output
write.csv("")
