### Calculate mQTLs for all probes within schizophrenia GWAS regions identified in the PGC GWAS as input for bayesian co-localistaion analysis
### Analysis is split by chromosomes

setwd("")

library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR; 

covariates_file_name = ""; 

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);



for(chr in 1:22){
	SNP_file_name = paste("", chr, ".txt", sep = "")
	snps_location_file_name = paste("", chr, ".txt", sep = "");
	expression_file_name = paste("", chr, ".txt", sep = "")

	output_file_name = paste("", chr, ".txt", sep = "")

	## threshold of results to save
	pvOutputThreshold = 1; ### need all mQTL results for each probe
	cisDist<-500000
	errorCovariance = numeric();

	 ## load SNP data
	snps = SlicedData$new();
	snps$fileDelimiter = "\t";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
	snps$LoadFile( SNP_file_name );

	## load methylation data
	gene = SlicedData$new();
	gene$fileDelimiter = "\t";      # the TAB character
	gene$fileOmitCharacters = "NA"; # denote missing values;
	gene$fileSkipRows = 1;          # one row of column labels
	gene$fileSkipColumns = 1;       # one column of row labels
	gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
	gene$LoadFile( expression_file_name);

	## Run the analysis
	snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
	snpspos[,1]<-rownames(snps)
	load("") ### load all illumina 450K annotation
	probeAnnot<-probeAnnot[rownames(gene), c("NAME", "CHR", "MAPINFO", "MAPINFO")]
	probeAnnot[,4]<-probeAnnot[,4]+1
	genepos<-probeAnnot

	### run eqtls
	me = Matrix_eQTL_main(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = output_file_name,
		pvOutputThreshold = 0,
		output_file_name.cis = output_file_name,
		pvOutputThreshold.cis = pvOutputThreshold,
		snpspos = snpspos, 
		genepos = genepos,
		cisDist = cisDist,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = "qqplot",
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE)
		
}


