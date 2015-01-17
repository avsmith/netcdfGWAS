Directory contains example genotypes to use as input for creation of netcdf files. The current example creates genotypes with mach/minimac. Documentation and details on these can be found 

* [MACH](http://www.sph.umich.edu/csg/abecasis/MACH/index.html)
* [minimac](http://genome.sph.umich.edu/wiki/Minimac)
* [minimac2](http://genome.sph.umich.edu/wiki/Minimac2)

Example data is taken from [http://www.sph.umich.edu/csg/cfuchsb/minimac_example.tgz]()

Steps to make impute dosages
---

1.  pre-phase with mach  
    `mach1 -d sample.dat -p sample.ped --rounds 20 --states 50 --phase --prefix sample.pp`
2.  create dosages with either minimac or minimac2  
  *  with minimac  
     `minimac --refHaps hapmap.hap --refSnps hapmap.snps --haps sample.pp.gz --snps sample.snps --prefix minimac.chr99.imp`  
  *  with minimac2  
    `minimac2 --refHaps hapmap.hap --refSnps hapmap.snps --haps sample.pp.gz --snps sample.snps --prefix minimac2.chr99.imp`  


