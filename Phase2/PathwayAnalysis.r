## script to perform gene ontology enrichment analysis on significant dmps
## uses logistic regression model to control for number of probes per gene

### create unique list of gene anno
uniqueAnno<-function(row){
if(row != ""){
	return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
	} else {
	return(row)
	}
	}
	

apply_tests<-function(each){

	## logistic regres.metasion approach
	pathway<-vector(length = nrow(bg_gene_go))
	pathway[grep(each, bg_gene_go[,2])]<-1
	model<-glm(pathway ~ gene_test + gene_size)
	p1<-summary(model)$coefficients[c("gene_test"),c(4)]
	p2<-summary(model)$coefficients[c("gene_size"),c(4)]
	or<-exp(summary(model)$coefficients[c("gene_test"),c(1)])
	beta<-summary(model)$coefficients[c("gene_size"),c(1)]
	if(or > 1){
		p1<-p1/2
	} else {
		p1<-1-(p1/2)
	}
		
	model<-glm(pathway ~ gene_test_ind + gene_size)
	p1.ind<-summary(model)$coefficients[c("gene_test_ind"),c(4)]
	p2.ind<-summary(model)$coefficients[c("gene_size"),c(4)]
	or.ind<-exp(summary(model)$coefficients[c("gene_test_ind"),c(1)])
	beta.ind<-summary(model)$coefficients[c("gene_size"),c(1)]
	if(or.ind > 1){
		p1.ind<-p1.ind/2
	} else {
		p1.ind<-1-(p1.ind/2)
	}
		
	return(unlist(c(each, names[match(each, names[,1]),2:3], sum(gene_size[which(pathway == 1)]), length(which(pathway == 1)), sum(gene_test[which(pathway == 1  & gene_test >= 1)]), length(which(pathway == 1 & gene_test >= 1)), p1,or,p2,beta,p1.ind,or.ind,p2.ind,beta.ind,  paste(bg_gene_go[which(pathway == 1 & gene_test >= 1),1], collapse = "|"))))
}

library(doParallel)

setwd("")
res.meta<-read.csv("", stringsAsFactors = FALSE, row.names = 1)

## as combination of 450K and EPIC probes need both annotation files
epicManifest<-read.csv("MethylationEPIC_v-1-0_B2.csv", skip = 7) ## load epic manifest file
epicManifest<-epicManifest[match(rownames(res.meta), epicManifest$Name),]
epicManifest$CHR<-as.character(epicManifest$CHR)
epicManifest$MAPINFO<-as.numeric(epicManifest$MAPINFO)
epicManifest$UCSC_RefGene_Name<-as.character(epicManifest$UCSC_RefGene_Name)
epicManifest$UCSC_RefGene_Group<-as.character(epicManifest$UCSC_RefGene_Group)
res.meta<-cbind(res.meta, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])

load("") ## load 450k manifest data
probeAnnot$CHR<-as.character(probeAnnot$CHR)
probeAnnot$MAPINFO<-as.numeric(probeAnnot$MAPINFO)
probeAnnot$UCSC_REFGENE_NAME<-as.character(probeAnnot$UCSC_REFGENE_NAME)
probeAnnot$UCSC_REFGENE_GROUP<-as.character(probeAnnot$UCSC_REFGENE_GROUP)
index<-which(is.na(res.meta$MAPINFO))
res.meta[index,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]<-probeAnnot[rownames(res.meta)[index], c("CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP")]

res.meta$CHR[which(res.meta$CHR == "X")]<-"23"
res.meta$CHR<-as.numeric(res.meta$CHR)
	
res.meta$UCSC_RefGene_Name<-unlist(lapply(res.meta$UCSC_RefGene_Name, uniqueAnno))

gene_test<-table(unlist(strsplit(res.meta$UCSC_RefGene_Name[which(res.meta$All_P_Random <= 9e-8)], "\\;")))
gene_size<-table(unlist(strsplit(res.meta$UCSC_RefGene_Name, "\\;")))

bg_genes<-names(gene_size)

### load gene pathway map 
gene_go<-read.csv("", stringsAsFactors = FALSE,)
# filter out genes not annotated to any pathway
gene_go<-gene_go[which(gene_go[,2] != ""),]

bg_gene_go<-gene_go[match(intersect(bg_genes, gene_go[,1]), gene_go[,1]),]
termCount<-table(unlist(strsplit(as.character(bg_gene_go[,2]), "\\|")))
terms<-names(termCount)[which(termCount > 9 & termCount <= 2000)]
names<-read.csv("", stringsAsFactors = FALSE) ## file to annotate pathway IDs with names

gene_size<-gene_size[bg_gene_go[,1]]
gene_test<-gene_test[bg_gene_go[,1]]
gene_test[is.na(gene_test)]<-0

gene_test_ind<-gene_test
gene_test_ind[which(gene_test_ind > 0)]<-1


cl<-makeCluster(30)
registerDoParallel(cl)

r<-foreach(i=1:length(terms), .combine=rbind, .export = c("gene_size", "gene_test_ind")) %dopar%{
	apply_tests(terms[i])
}
colnames(r)<-c("ID", "Name", "Type", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "P:GenesinTestList", "OR", "P:GeneSize", "Beta:GeneSize","P:GenesinTestList", "OR", "P:GeneSize", "Beta:GeneSize", "GenesinTestListAndPathway")

### filter to terms with between 10 and 2000 genes
r<-as.data.frame(r, stringsAsFactors = FALSE)
r<-r[which(as.numeric(r[,5]) < 2000 & as.numeric(r[,5]) > 9),]

write.csv(r, "")

### identify independent significant terms

## rank analysis by p value and filter to significant terms
r<-r[order(as.numeric(r[,12])),]
r<-r[which(as.numeric(r[,12]) < 0.05/nrow(r)),]

output<-c()
while(!is.null(r)){

	if(class(r) != "character"){

	### for all terms repeat analysis controlling for most significant terms
	best_term<-vector(length =  nrow(bg_gene_go))
	best_term[grep(r[1,1], bg_gene_go[,2])]<-1
	merge.id<-c()
	merge.name<-c()
	remove<-c()
	
	for(j in 2:nrow(r)){
	
		pathway<-vector(length = nrow(bg_gene_go))
		pathway[grep(r[j,1], bg_gene_go[,2])]<-1
		## error if all genes in test list and pathway overlap with best_term
		#if(sum(
		model<-glm(pathway ~ gene_test_ind + gene_size + best_term)
		if(!is.na(summary(model)$coefficients["gene_test_ind", "Pr(>|t|)"]) & summary(model)$coefficients["gene_test_ind", "Pr(>|t|)"] > 0.05){
			merge.id<-append(merge.id, r[j,1])
			merge.name<-append(merge.name, r[j,2])
			remove<-append(remove, j)
		} 
	
	}
	merge.id<-paste(unlist(merge.id), collapse = "|")
	merge.name<-paste(unlist(merge.name), collapse = "|")

	output<-rbind(output, c(r[1,], merge.id, merge.name))
	
	r<-r[-c(1, remove),]
	} else {
	
	output<-rbind(output, c(r, "", ""))
	r<-NULL
}
}
	
write.csv(output, "")	
