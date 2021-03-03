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

