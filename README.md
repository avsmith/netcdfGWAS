Genetic Analysis with netcdf
======

The scripts in this directory give examples of how to use netcdf for genetic analyses, notably genetic association analysis based on imputed genotypes. While there are several excellent packages for association analysis, they typically only allow relatively standard regression (linear, logistic, coxph). As dedicated packages they can run relatively fast. However, these fixed programs don't allow for flexible models (GEE, for example.) 

The approach outlined in these example scripts allows for much greater flexibility in approach, though this comes at the expense of speed, as the R based analysis is typically much slower than dedicated programs.