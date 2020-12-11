Welcome.  Thank you for using PSMix program.
This program uses maximum likelihood approach (EM algorithm) to infer 
population structure via mixture model using multilocus genotype data.  


1. Install

For windows:
Download PSMix_1.0-6.zip into your local hard drive.
Open Rgui (R graphical interface), go to menu "Packages\install packages
from local (zip) files ...", then select the provided zip file.  Go to menu 
"Packages\load package..." to load the package into R for use.  

For Linux/Unix:
Download PSMix_1.0-6.tar.gz into your local hard drive.
Install the R package from the source file using the following commands:
first go to the directory where the downloaded file is stored and then use
"R CMD INSTALL PSMix_1.0-6.tar.gz". For more options please see R help on
"INSTALL".


2. Input

Genotype data should be stored in a dataframe or matrix format, with 
columns being marker information and two consecutive rows denoting 
the genotypes of each individual. 

The following is a sample input file format 
234  130  na  227  ...  
232  130  na  217  ...  

where the first column is the allele information for the first marker, 
and the first two rows represent the genotype information for first 
individual. "na" for missing.  


3. Usage

PSMix(K = 2, Geno, itMax = 10000, eps = 1e-06, seed = 1, MarkerVar = FALSE)

Arguments
K: number of underlying populations  
Geno: genotype data  
itMax: maximum number of iterations in EM algorithm  
eps: convergence criterion  
seed: set up the seed for random number generation (reproducibility)  
MarkerVar: allow individual marker admixture.  This feature is not fully 
           implemented in this version.  

There is also a demo illuminating the use of the program using a subset of 
Pima-Surui dataset.  To see the demo, simply type demo(struc).  


4. Output

A list containing the following components:

PIk: estimated probability of belonging to each subpopulation for each individual 
     (at each marker if MarkerVar is set "True" which corresponding the third model
     in the paper).
Zimak: estimated probability of belonging to each population for each allele 
       (estimated origin of allele for ith individual, at mth marker, on ath 
       allele copy). 
Gkmj: estimated allele frequencies for each subpopulation.  


5. Credit 

This program is developed based on the algorithm proposed in 
Nianjun Liu, Baolin Wu and Hongyu Zhao (2004). "Inference of Population 
Structure Using Mixture Model". Technical report, Division of Biostatistics, 
Department of Epidemiology and Public Health, Yale University. Please cite 
this paper if you use this program in your research for publication.
 

6. Support 

All questions and comments should be direct to Baolin Wu or Nianjun Liu.  
Email: baolin@biostat.umn.edu, nianjun.liu@yale.edu.  

