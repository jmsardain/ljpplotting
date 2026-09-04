# Plotting script for Lund jet plane analysis

The code is structured as follows : 

- The main function is in `plotting.cc` 
- Helper functions are in `include/` 
- The error names and which histos to use to determine them are in `errorNames.txt` 

To compile the code you need : 

```
setupATLAS 
lsetup "root 6.36.02-x86_64-el9-gcc14-opt"

source compile.sh 
```

## 1D distributions (differential bin)
To run the code : 
```
./draw1DHistos errorNames.txt  
```

This code will create 2 pdf files :

- `kt_1D_dataMC.pdf` contains 1D histos, ratio to data 
- `kt_1Derror_dataMC.pdf` contains 1D histos, ratio to data, and relative uncertainties for errors defined in `errorNames.txt` 

Improvement (maybe not necessary) : the code expects the errors from `errorNames.txt`, if error names are not the same, it crashes. This becomes important when looking at the breakdown of each systematic. This is an easy fix. 

## Relative uncertainty plot for all bins 

This plot is a bit too crowded, but useful to see behavior across all the bins 

``` 
./drawRelUnc errorNames.txt 
``` 

This code will create one pdf file : 

- `totalUnc_allBins.pdf` contains the relative uncertainty for all of the errors in `errorNames.txt` 