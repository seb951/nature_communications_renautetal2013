#!/usr/bin/Rscript --verbose

args = commandArgs(TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

library(ape)

for(a in from:to)
{

if(a == 1) fst_map_match = read.delim("results/ann_pet_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(a == 2) fst_map_match = read.delim("results/ann_deb_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(a == 3) fst_map_match = read.delim("results/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst
if(a == 4) fst_map_match = read.delim("results/pet_deb_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(a == 5) fst_map_match = read.delim("results/pet_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(a == 6) fst_map_match = read.delim("results/deb_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst


###calculate Moran's I within and between chromosomes###
spatial_res_17_chromo = matrix(0,nrow = 17, ncol = 21)
colnames(spatial_res_17_chromo) = c("chromosome","fix_observed", "fix_expected","fix_sd","fix_pval", "top1_observed", "top1_expected","top1_sd","top1_pval", "top3_observed", "top3_expected","top3_sd","top3_pval", "top5_observed", "top5_expected","top5_sd","top5_pval", "prop_obs_095","prop_obs_top1","prop_obs_top3","prop_obs_top5")
spatial_res_17_chromo[,1] = c(1:17)
col = ncol(fst_map_match)
fst_map_match[is.na(fst_map_match[,(col-8)]),(col-8):(col-6)] = 0
fst_map_match[,(col-8)] = jitter(fst_map_match[,(col-8)],0.1) #add/remove between -/+ 2e-05. to all fst values 
top1 = quantile(fst_map_match[,(col-8)],0.99, na.rm =T)
top3 = quantile(fst_map_match[,(col-8)],0.97, na.rm =T)
top5 = quantile(fst_map_match[,(col-8)],0.95, na.rm =T)
aaa = Sys.time()


for(chr in 1:17)
#for(chr in 1:1)
	{
		dataset_mapped = cbind(0,fst_map_match[as.numeric(fst_map_match[,col-4]) == chr,col-8],fst_map_match[as.numeric(fst_map_match[,col-4]) == chr,col-3],0,0,0)

		dataset_mapped[dataset_mapped[,2] >= 0.95,1] = 1 # fst fixed
		dataset_mapped[dataset_mapped[,2] >= top1,4] = 1 # fst top1
		dataset_mapped[dataset_mapped[,2] >= top3,5] = 1 # fst top3
		dataset_mapped[dataset_mapped[,2] >= top5,6] = 1 # fst top5
		
		dists = as.matrix(dist(cbind(dataset_mapped[,2], dataset_mapped[,3]))) # construct a distance matrix. EUCLIDEAN
		
		dists_round_3decimals = round(dists,2) # round up so any two markers at the same position are at 0 distance.

		#diag(dists) = 0.00000001

		#for(j in 1:nrow(dists))
		#	{
		#	x = c(1:nrow(dists))[dists[,j] == 0]
		#	dists[x,j] = 0.0000001
		#	}

		dists.inv = 1 / dists_round_3decimals
		diag(dists.inv) = 0
		dists.inv.2 = gsub(Inf,100,dists.inv) #anything above 100 is essentially one hundred. (because at most markers are 0.001 cM apart, thus max(inverse of distance) = 100)
		dists.inv.final = matrix(as.numeric(dists.inv.2), nrow = nrow(dists.inv), ncol = nrow(dists.inv))
		x = Moran.I(dataset_mapped[,1], dists.inv.final , na.rm = T, alternative = "greater")
		x_top1 = Moran.I(dataset_mapped[,4], dists.inv.final , na.rm = T, alternative = "greater")
		x_top3 = Moran.I(dataset_mapped[,5], dists.inv.final , na.rm = T, alternative = "greater")
		x_top5 = Moran.I(dataset_mapped[,6], dists.inv.final , na.rm = T, alternative = "greater")
		spatial_res_17_chromo[chr,2:5] = cbind(signif(x[[1]],4),signif(x[[2]],4),signif(x[[3]],4),signif(x[[4]],4))
		spatial_res_17_chromo[chr,6:9] = cbind(signif(x_top1[[1]],4),signif(x_top1[[2]],4),signif(x_top1[[3]],4),signif(x_top1[[4]],4))
		spatial_res_17_chromo[chr,10:13] = cbind(signif(x_top3[[1]],4),signif(x_top3[[2]],4),signif(x_top3[[3]],4),signif(x_top3[[4]],4))
		spatial_res_17_chromo[chr,14:17] = cbind(signif(x_top5[[1]],4),signif(x_top5[[2]],4),signif(x_top5[[3]],4),signif(x_top5[[4]],4))
		spatial_res_17_chromo[chr,18:21] = cbind(signif(sum(dataset_mapped[,1]) / nrow(dataset_mapped),4),signif(sum(dataset_mapped[,4]) / nrow(dataset_mapped),4),signif(sum(dataset_mapped[,5]) / nrow(dataset_mapped),4),signif(sum(dataset_mapped[,6]) / nrow(dataset_mapped),4))
		
		print(paste("comparison",chr,"___ of 17 chromosome done", Sys.time(),"___",nrow(dataset_mapped)))
	}

if(a == 1) write.table(spatial_res_17_chromo,"results/ann_pet_moran", row.names = F)
if(a == 2) write.table(spatial_res_17_chromo,"results/ann_deb_moran", row.names = F)
if(a == 3) write.table(spatial_res_17_chromo,"results/ann_arg_moran", row.names = F)
if(a == 4) write.table(spatial_res_17_chromo,"results/pet_deb_moran", row.names = F)
if(a == 5) write.table(spatial_res_17_chromo,"results/pet_arg_moran", row.names = F)
if(a == 6) write.table(spatial_res_17_chromo,"results/deb_arg_moran", row.names = F)

}
#
#top 3 or 5 % of the genome!
#graphs
library(ggplot2)


m1 = read.delim("results/ann_pet_moran", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
m2 = read.delim("results/ann_deb_moran", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
m3 = read.delim("results/ann_arg_moran", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst
m4 = read.delim("results/pet_deb_moran", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
m5 = read.delim("results/pet_arg_moran", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
m6 = read.delim("results/deb_arg_moran", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst


###
all = cbind(c(rep("ann_pet",17),rep("ann_deb",17),rep("ann_arg",17),rep("pet_deb",17),rep("pet_arg",17),rep("deb_arg",17)),rbind(m1,m2,m3,m4,m5,m6))
all4 = cbind(c(rep("ann_pet",17),rep("ann_deb",17),rep("deb_arg",17),rep("pet_arg",17)),rbind(m1,m2,m6,m5))

fit_all4 <- lm(all4[,11]~all4[,1], data=all4)
anova(fit_all4)

all = cbind(c(rep(1,17),rep(2,17),rep(3,17),rep(4,17),rep(5,17),rep(6,17)),rbind(m1,m2,m3,m4,m5,m6))
all4 = cbind(c(rep(1,17),rep(2,17),rep(3,17),rep(4,17)),rbind(m1,m2,m6,m5))


#p = ggplot(all[c(1,11)],aes(all[,1],all[,11])); p +  geom_boxplot(aes(colour = factor(all[,1]))) + geom_point(aes(colour = factor(all[,1])))
p = ggplot(all4[,c(1,11)],aes(all4[,1],all4[,11])); p +  geom_boxplot(aes(colour = factor(all4[,1]))) + geom_point(aes(colour = factor(all4[,1])))

boxplot(all4[,11] ~ all4[,1], data = all4, col = "lightgray")
dev.print(device=svg, "~/Desktop/Moran_top3_4comp.svg", onefile=FALSE)
dev.off()



