# gghatvr4

## example gghatvr4.inp file (included in repository)

## description:
calculates GRM matrices (genomic relationship matrices) using VanRaden (2008) methods 1 or 2. And inverse of G
input formats: 
. LDMIP
. general marker genotypes (produces re-coded genotype matrix)
. outputfile from Beagle (.phased file)
. PLINK format (.bed file)

##compilation: 
ifort -r8  gghatvr4.f90  -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -o gghatvr4


# gwablup_wgts.jl
Julia script to calculate GWABLUP weights, that can be used in gghatvr4 to perform GWABLUP.   

Example usage:   
julia gwablup_wgts.jl gwas.lrt  0.001  5   
gwas.lrt = input file with a single column of log-likelood-ratio-values from GWAS (SNPs are in map order; 1 value per SNP)   
0.001 = PI-value to be used   
5 = number of SNPs used in moving average   
output file = gwablupwgts.out (2 columns : SNPweight & smoothedLRTvalue)
