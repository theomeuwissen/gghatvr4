# example of usage: 
# julia gwablup_wgts.jl gwas.lrt  0.01  5
# gwas.lrt = input file with a single column of log-likelood-ratio-values from GWAS (SNPs are in map order; 1 value per SNP)
# 0.01 = PI-value to be used
# 5 = number of SNPs used 
# output file = gwablup_wgts.out (2 columns : SNPweight & smoothed_LRT_value) 
using DelimitedFiles, Statistics
dat=readdlm(ARGS[1])
PI=parse(Float64,ARGS[2])
Navg=parse(Int,ARGS[3])
Nleft=round(Int,(Navg-1)/2)
if(Nleft<0)
    Nleft=0
end    
println("no of SNPs ",size(dat,1))
println("PI= ",PI)
println("No SNPs to left (and right) included in smoothed average ",Nleft)        
avg=zeros(size(dat,1),2);
for i=1:size(dat,1)
	   avg[i,2]=mean(dat[max(1,i-Nleft):min(i+Nleft,size(dat,1))])
	   avg[i,1]=PI*exp(avg[i,2])/(PI*exp(avg[i,2])+1-PI)
end
writedlm("gwablup_wgts.out",avg)
exit()
