# Script to create netcdf files from mach imputations

# Example data created via mach/minimac2
# Note: mach and minimac installed for Mac OS X via: https://github.com/avsmith/homebrew-genetics 
# Example data source: http://www.sph.umich.edu/csg/cfuchsb/minimac2_example.tgz
# Note: example data includes linux binaries for mach, minimac, and minimac2
# mach command to phase:
# > mach1 --datfile minimac2_example/sample.dat --pedfile minimac2_example/sample.ped --phase --prefix minimac2_example/sample.pp
# minimac2 command to impute:
# > minimac2 --refHaps minimac2_example/hapmap.hap --refSnps minimac2_example/hapmap.snps --haps minimac2_example/sample.pp.gz --snps minimac2_example/sample.snps --prefix minimac2_example/chr99.imp

# NOTE
# The example below requires a position for each SNP. While not required for analysis, convinient
# to store within the netcdf object (making it self-contained with all relevant info)
# The fake postions for this example was create in R as follows:
# > library(dplyr)
# > library(tidyr)
# > chr99info <- read.table("minimac2_example/chr99.imp.info",header=TRUE)
# > mlinfoNames <- names(chr99info)
# > chr99info <- chr99info %>% mutate(POS=extract_numeric(SNP)) %>% select(one_of(c("SNP","POS",mlinfoNames[2:length(mlinfoNames)])))
# > write.table(chr99info, "minimac2_example/chr99.imp.info.pos",quote=FALSE,row.names=FALSE)

options(stringsAsFactors = FALSE)

# NOTE:
# IDs for the individuals have been extracted, and limited to integeters
# The IDs in the netcdf files contains characters.
# netcdf does not work well with characters

idfile<-"minimac2_example/sampleIDs.txt" 
ids<-read.table(idfile, header=TRUE)

# NOTE:
# Following not used in the example. The netcdf object has a field for 
# SNPs used for impuation. Commented out code could be used for that.
# imp <- read.table("snps_used_for_imputation.txt")
# names(imp) <- c("SNP","IMP")
# row.names(imp) <- imp$SNP

# Function call used to create one netcdf file for each chromosome
read.dosage<-function(chromosome){
  
  # NOTE
  # Am using ncdf rather than ncdf4. When installing use version found at:
  # http://cirrus.ucsd.edu/~pierce/ncdf/ 
  # This site has v1.8.6. CRAN has version 1.6.6, which has a bug for reading 
  # Missing values stored as byte
  
  require(ncdf)
  
  # Read in the dosages made as the example above
  imputed<-sub("@",as.character(chromosome),"minimac2_example/chr@.imp.dose")
  
  # NOTE
  # This is the info file from minimac, with a column for SNP position added
  snpfile<-sub("@",as.character(chromosome), "minimac2_example/chr@.imp.info.pos")
  snps<-read.table(snpfile,header=TRUE)
  
  NSNPS<-nrow(snps)
  NID<-nrow(ids)
  
  # Create the dimensions for the netcdf object
  snpdim<-dim.def.ncdf("snp","bases", 1:NSNPS)
  sampledim<-dim.def.ncdf("id","th",1:NID, unlim=TRUE)
  
  # NOTE:
  # IDs are stored as integers!
  varID<-var.def.ncdf("sampleID","person",dim=sampledim, missval=0, prec="integer")
  # Position of the SNP
  varSNP<-var.def.ncdf("pos","bp",dim=snpdim, missval=0, prec="integer")
  
  # Dosages and other numbers stored as "float"
  # R does not support floats, but netcdf does. This makes for much smaller output files.
  # When reading back in, the values are converted to doubles for analysis in R.
  # Nothing is lost with float/double given the precision of the impuations.
  vargenotype<-var.def.ncdf("dose","dose",dim=list(snpdim, sampledim), missval=-1, prec="float")
  
  varrefallele<-var.def.ncdf("ref","base",dim=snpdim,missval=0,prec="byte")
  varaltallele<-var.def.ncdf("alt","base",dim=snpdim,missval=0,prec="byte")
  varsnpid<-var.def.ncdf("rsid","base",dim=snpdim,missval=0,prec="integer")
  vareaf<-var.def.ncdf("EAF","freq",dim=snpdim,missval=-1,prec="float")
  varmaf<-var.def.ncdf("MAF","freq",dim=snpdim,missval=-1,prec="float")
  varrsq<-var.def.ncdf("Rsq","rsq",dim=snpdim,missval=-1,prec="float")
  varimp<-var.def.ncdf("UsedImp","imp",dim=snpdim,missval=-1,prec="byte")
  
  outfname<-sub("@",as.character(chromosome),"netcdf_files/chr@.mldose.ncdf")
  
  
  outfile<-create.ncdf(outfname,
                       list(varID,varSNP,varrefallele, 
                            varaltallele,varsnpid,vareaf,varmaf,varrsq,vargenotype,varimp),
                       verbose=FALSE)
  
  put.var.ncdf(outfile,varSNP,as.integer(snps$POS),start=1, count=NSNPS)
  
  # SNPs are stored as integers, normally as rsID.
  # For this example, also have integers normally prefixed by SNP (from the example)
  put.var.ncdf(outfile,varsnpid,as.integer(sub("^(rs|SNP)","",snps$SNP)),start=1, count=NSNPS)
  put.var.ncdf(outfile, varID,ids$id,start=1, count=NID)
  put.var.ncdf(outfile, varrefallele, match(snps$Al1, c("A","C","G","T"), nomatch=0) ,start=1, count=NSNPS)
  put.var.ncdf(outfile, varaltallele, match(snps$Al2, c("A","C","G","T","-"), nomatch=0) ,start=1, count=NSNPS)
  put.var.ncdf(outfile, vareaf, snps$Freq1,start=1, count=NSNPS)
  put.var.ncdf(outfile, varmaf, snps$MAF,start=1, count=NSNPS)
  put.var.ncdf(outfile, varrsq, snps$Rsq,start=1, count=NSNPS)
  
  # Commented out. This would normally put a flag for imputation.
  # put.var.ncdf(outfile, varimp, imp[snps$SNP,]$IMP,start=1, count=NSNPS)
  
  genotypes<-file(imputed,"r")
  
  header<- 1
  
  print(paste(NSNPS, "SNPs to process"))
  for(i in 1:NID) {
    snpi<-scan(genotypes,what=character(0), nlines=1,quiet=TRUE, skip=0)
    header<-0
    # The key bit to install 
    put.var.ncdf(outfile,vargenotype, vals=as.double(snpi[3:length(snpi)]),start=c(1,i),count=c(NSNPS,1))
    if (!(i%%100)) print(paste("individual #",i," of ",NID,sep=""))
  }
  close(genotypes)
  close(outfile)
}

# Example is a chr99 test. Normally set as 1:22 (or 1:23)
# This iterates over all the files
for(i in 99:99){
  print(paste("Starting: chr",i,sep=""))
  read.dosage(i)
}


