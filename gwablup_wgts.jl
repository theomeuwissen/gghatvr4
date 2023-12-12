# example of usage: 
# julia gwablup_wgts.jl gwas.lrt  0.001  5
# gwas.bhat.se = input GWAS results file: 2 columns: SNPeffect & stand.error (SNPs are in map order)
# 0.001 = PI-value to be used
# 5 = number of SNPs used 
# output file = gwablup_wgts.out (2 columns : SNPweight & smoothed_LR_value) 
using DelimitedFiles, Statistics
dat=readdlm(ARGS[1])
lr=0.5*(dat[:,1]./dat[:,2]).^2
PI=parse(Float64,ARGS[2])
Navg=parse(Int,ARGS[3])
Nleft=round(Int,(Navg-1)/2)
if(Nleft<0)
    Nleft=0
end    
println("no of SNPs ",size(lr,1))
println("PI= ",PI)
println("No SNPs to left (and right) included in smoothed average ",Nleft)        
avg=zeros(size(lr,1),2);
for i=1:size(lr,1)
	   avg[i,2]=mean(lr[max(1,i-Nleft):min(i+Nleft,size(lr,1))])
	   avg[i,1]=PI*exp(avg[i,2])/(PI*exp(avg[i,2])+1-PI)
end
writedlm("gwablup_wgts.out",avg)
exit()
