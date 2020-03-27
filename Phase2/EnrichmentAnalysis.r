## script to perform enrichment analysis of dmps
## selects maximal background of probes to match distribution of mean and variance in DNA methylation of dmps. 

setwd("")
res.meta<-read.csv("", stringsAsFactors = FALSE, row.names = 1) ## load meta-analysis results

## calc probe means and variances to subsample 

load("EUGEI/Normalised.rdat") ## load largest EPIC dataset 
mu<-apply(betas, 1, mean)
sigma<-apply(betas, 1, sd)
load("Aberdeen/Normalised.rdat") ## load largest 450K dataset 
mu.aber<-apply(betas.aber, 1, mean)
sigma.aber<-apply(betas.aber, 1, sd)

mu.comb<-mu[match(rownames(res.meta), rownames(betas))]
naIndex<-which(is.na(mu.comb))
mu.comb[naIndex]<-mu.aber[match(rownames(res.meta)[naIndex], rownames(betas.aber))]

sigma.comb<-sigma[match(rownames(res.meta), rownames(betas))]
naIndex<-which(is.na(sigma.comb))
sigma.comb[naIndex]<-sigma.aber[match(rownames(res.meta)[naIndex], rownames(betas.aber))]


## sample probes to match variance distribution include maximal possble number of samples
## nb done across all probes in analysis (some enrichment analyses only on 450K)
sigma.comb.tmp<-sigma.comb
probeBins.sigma<-cut(sigma.comb.tmp, quantile(sigma.comb.tmp, seq(0,1,0.1), na.rm = TRUE))
mu.comb.tmp<-mu.comb
probeBins.mu<-cut(mu.comb.tmp, quantile(mu.comb.tmp, seq(0,1,0.1), na.rm = TRUE))
tabProbeBins<-table(probeBins.sigma, probeBins.mu)
nSamples<-table(probeBins.sigma[which(res.meta$All_P_Random < 9e-8)], probeBins.mu[which(res.meta$All_P_Random < 9e-8)])
# exclude 0s from division
nSamples[which(nSamples == 0)]<-NA
## find maximal number of times we can select a subset of probes matched for variance and mean
nTimes<-floor(min(tabProbeBins/nSamples, na.rm = TRUE))
nSelect<-nTimes*nSamples

subSample<-NULL
for(j in 1:ncol(nSamples)){
	for(i in 1:nrow(nSamples)){
		if(!is.na(nSelect[i,j])){
			subSample<-c(subSample, sample(which(probeBins.mu == colnames(nSamples)[j] & probeBins.sigma == rownames(nSamples)[i]), nSelect[i,j]))
		}
	}
}

## are top associations more likely to be genetically mediated?
## load heritability statistics
herit<-read.csv("", stringsAsFactors = FALSE, row.names = 1)

herit.sub<-herit[rownames(res.meta)[which(res.meta$All_P_Random < 9e-8)],]
herit<-herit[rownames(res.meta)[which(res.meta$All_P_Random > 9e-8)],]
d.a<-density(herit.sub$A., from = 0, to = 1, na.rm = TRUE)
d.c<-density(herit.sub$C., from = 0, to = 1, na.rm = TRUE)
d.e<-density(herit.sub$E., from = 0, to = 1, na.rm = TRUE)

d.a.matched<-density(herit$A.[subSample], from = 0, to = 1, na.rm = TRUE)
d.c.matched<-density(herit$C.[subSample], from = 0, to = 1, na.rm = TRUE)
d.e.matched<-density(herit$E.[subSample], from = 0, to = 1, na.rm = TRUE)

pdf("", width = 12, height = 5)
par(mfrow = c(1,3))
plot(c(0,d.a$x, 1), c(0,d.a$y, 0), col = rgb(1,0,0,0.2), main = "", xlab = "A estimate", , cex.axis = 1.5, cex.lab = 1.5, type = "l", ylab = "Density", ylim = c(0,max(c(d.a$y, d.a.matched$y))))
polygon(c(0,d.a$x, 1), c(0,d.a$y, 0), col = rgb(1,0,0,0.2), border = rgb(1,0,0,0.2))
polygon(c(0,d.a.matched$x, 1), c(0,d.a.matched$y, 0), col = rgb(0,1,0,0.2), border = rgb(0,1,0,0.2))
legend("topright", pch = 15, c("Schizophrenia DMPs","Matched background"), col = c(rgb(1,0,0,0.2), rgb(0,1,0,0.2)), cex = 1.4) 

plot(c(0,d.c$x, 1), c(0,d.c$y, 0), col = rgb(1,0,0,0.2), main = "", xlab = "C estimate",  cex.axis = 1.5, cex.lab = 1.5, type = "l", ylab = "Density", ylim = c(0,max(c(d.c$y, d.c.matched$y))))
polygon(c(0,d.c$x, 1), c(0,d.c$y, 0), col = rgb(1,0,0,0.2), border = rgb(1,0,0,0.2))
polygon(c(0,d.c.matched$x, 1), c(0,d.c.matched$y, 0), col = rgb(0,1,0,0.2), border = rgb(0,1,0,0.2))

plot(c(0,d.e$x, 1), c(0,d.e$y, 0), col = rgb(1,0,0,0.2), main = "", xlab = "E estimate",cex.axis = 1.5, cex.lab = 1.5, type = "l", ylab = "Density", ylim = c(0,max(c(d.e$y, d.e.matched$y))))
polygon(c(0,d.e$x, 1), c(0,d.e$y, 0), col = rgb(1,0,0,0.2), border = rgb(1,0,0,0.2))
polygon(c(0,d.e.matched$x, 1), c(0,d.e.matched$y, 0), col = rgb(0,1,0,0.2), border = rgb(0,1,0,0.2))
dev.off()

wilcox.test(herit.sub$A., herit$A.[subSample])
wilcox.test(herit.sub$C., herit$C.[subSample])
wilcox.test(herit.sub$E., herit$E.[subSample])

