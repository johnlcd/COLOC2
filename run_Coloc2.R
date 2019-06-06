args <- commandArgs(trailingOnly=TRUE)
gwas_file <- args[1]
eqtl_file <- args[2]
outfolder <- args[3]
type <- args[4] #"cc" or "quant"
NT <- as.numeric(args[5])
NC <- NT/2 
p12 <- 1e-6

source("/home/chenjiabin/tools/COLOC2/coloc_functions.R")
source("/home/chenjiabin/tools/COLOC2/functions_coloc_likelihood_summary_integrated.R")
library(data.table)

biom.df <- formatColoc(fname = gwas_file, type=type, N=NT, Ncases=NC, info_filter=0.6, maf_filter=0.05, fread=T, eqtl=FALSE)
eqtl.df <- formatColoc(fname = eqtl_file, type=type, N=NT, Ncases=NA, info_filter=0.6, maf_filter=0.05, fread=T, eqtl=TRUE)
biom.df$type <- type
eqtl.df$type <- type
#eqtl.df$MAF <- 0.5

#names(biom.df) = c("CHR", "POS", "SNP", "N", "A1", "A2", "Z", "BETA", "SE", "N_CONTROLS")
#biom.df$Ncases = biom.df$N-biom.df$N_CONTROLS
#biom.df$PVAL = 2*pnorm(-abs(biom.df$Z))

names(eqtl.df)[which(names(eqtl.df)=="REF")] <- "A1"
names(eqtl.df)[which(names(eqtl.df)=="ALT")] <- "A2"
names(eqtl.df)[which(names(eqtl.df)=="beta")] <- "BETA"
names(eqtl.df)[which(names(eqtl.df)=="se.beta")] <- "SE"


# If this takes too long, submit one script per gene by specifying character vector for list.probes
coloc2 <- coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=TRUE, outfolder="outfolder", prefix="", plot=FALSE, save.coloc.output=FALSE, match_snpid=T)
#coloc.eqtl.biom(merged.data, list.probes=NULL, useBETA=TRUE, outfolder, prefix= prefix, min_snps=50, p1=1e-4, p2=1e-4, p12=1e-6, moloc=FALSE)

