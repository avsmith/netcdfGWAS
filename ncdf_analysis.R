# Run example analysis using netcdf stored genotypes (as created in "create_ncdf.R")
# This example runs linear regression, but the modeling functions can be changed to any
# analysis. The capture values for the results would need to be changed appropriately.
library(ncdf)

# Simulate phenotype data
nsamples <- 500
prevalence <- 0.1
set.seed(42)
pheno <- data.frame(id=1:nsamples,OUTCOME=sample(c(0,1),nsamples,replace=TRUE,prob=c(1-prevalence,prevalence)),AGE=sample(50:70,nsamples,replace=TRUE),SEX=sample(c(1,2),nsamples,replace=TRUE))

# Filtering thresholds for exclusion. Set to zero if all data wanted.
rsq_threshold <- 0.5
maf_threshold <- 0.01

# Seed the initial model
formula<-as.formula("OUTCOME ~ AGE + SEX" )
id <- pheno[["id"]]
m <- glm(formula,family=binomial(), ,data = pheno, na.action =na.exclude)

# Cycle through the chromosomes
for (CHR in 99:99){
  # Open NetCDF file
  genonc<-open.ncdf(sub("@",as.character(CHR),"netcdf_files/chr@.mldose.ncdf"))

  # Get list of SNPs
  SNPs <- get.var.ncdf(genonc,"rsid")
  # Get list of IDs
  IDs <- get.var.ncdf(genonc,"sampleID")
  # Number of SNPs
  NSNP<-dim(SNPs)
  # Number of Samples
  NID<-dim(IDs)

  # Placeholders for results
  betasnp<-rep(NA,NSNP)
  sesnp<-rep(NA,NSNP)
  n<-rep(NA,NSNP)
  pzsnp<-rep(NA,NSNP)
  # While the netcdf holds EAF for the entire sample, this is calculate for those in the analysis
  eaf1 <- rep(NA,NSNP)
  ref<-c("A","C","G","T")[get.var.ncdf(genonc,"ref")]
  alt<-c("A","C","G","T")[get.var.ncdf(genonc,"alt")]
  # ids stored as intergers, paste with "rs" to create rsIDs (as seen below)
  id<-get.var.ncdf(genonc,"rsid")

  pos<-get.var.ncdf(genonc,"pos")
  usedimp<-get.var.ncdf(genonc,"UsedImp")
  EAF<-get.var.ncdf(genonc,"EAF")
  MAF<-get.var.ncdf(genonc,"MAF")
  Rsq<-get.var.ncdf(genonc,"Rsq")

  # Used to compare order in dosages and input samples
  keep <- match(pheno$id, IDs)

  for(i in 1:10){
    dose <- get.var.ncdf(genonc,"dose",start=c(i,1),count=c(1,NID))
    dose <- dose[keep]

    if(Rsq[i]>=rsq_threshold & MAF[i]>=maf_threshold){

      gm  <-  update(m,.~ . + dose,na.action =na.omit)

      betasnp[i] <- coef(summary(gm))["dose","Estimate"]
      sesnp[i] <- coef(summary(gm))["dose","Std. Error"]
      pzsnp[i] <- coef(summary(gm))["dose","Pr(>|z|)"]
      n[i] <- length(resid(gm))
      eaf1[i] <- (mean(dose,na.rm=TRUE)/2)

    }
    if(i%%100 == 0) print(paste("SNP nr: ",i," of ",NSNP,sep=""))
  }

  # Write out the results
  chrom <- rep(CHR,NSNP)
  results <- data.frame(chr=chrom, snp=paste("rs",id,sep=""), pos, ref, alt,
                        freq1 = signif(eaf1,3), strand ="+",N=n,
                        beta_SNP=signif(betasnp,5),
                        se_SNP=signif(sesnp,5),
                        pval_SNP=signif(pzsnp,5), Rsq=signif(Rsq,3),UsedForImputation=usedimp)


  write.table(results,gzfile(paste0("output_dir/Analysis_chr",CHR,".txt.gz")),sep="\t",quote=FALSE,row.names=FALSE)

}