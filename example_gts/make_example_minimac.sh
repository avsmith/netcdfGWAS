#!/bin/sh

mach1 -d sample.dat -p sample.ped --rounds 20 --states 50 --phase --prefix sample.pp
minimac2 --refHaps hapmap.hap --refSnps hapmap.snps --haps sample.pp.gz --snps sample.snps --prefix minimac2.chr99.imp
