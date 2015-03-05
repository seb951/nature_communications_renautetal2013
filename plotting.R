#!/usr/bin/Rscript --verbose


#args = commandArgs(TRUE)
#a = as.numeric(args[1])


setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

####################
######
### MAPPING FST VALUES ###
##########################
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",]
map = cbind(map,0,0)
colnames(map) = c("name","LG","map_centiMorgan","unique_position","nbnuc")

ref = read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F, stringsAsFactors = F) # new reference assembly
len = cbind(gsub(">","",ref[seq(from = 1, to = nrow(ref), by = 2),]),apply(ref, 1,nchar)[seq(from = 2, to = nrow(ref), by = 2)])

#unique position for the markers#
max_position = cbind(c(1:17), 0,0) # the size of each chromosome and the cumulative number of centimorgans.

for(i in 1:17)
	{
	max_position[i,2] = max(as.numeric(map[as.numeric(map[,2]) == i,3]))
	max_position[i,3] = sum(max_position[1:i,2])-max_position[i,2]
	}


for(i in 1:nrow(map))
#for(i in 1:200)
	{
	for(j in 1:17)
		{
		if(as.numeric(map[i,2]) == max_position[j,1]) map[i,4] = (as.numeric(map[i,3]) + sum(max_position[1:j,2]) - sum(max_position[j,2]) + 1)
		}
		map[i,5] = len[map[i,1] == len[,1],2]
	if(i %% 1000 == 0) print(paste(i,"of",nrow(map),Sys.time()))
}
max_position[,3] = (max_position[,3] +1)


par(mfrow = c(4,1))
for(a in c(1,2,6,5))
{
################
### PLOTTING ###
################
#if(a == 7) c1 =  read.table("results_6species/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

###fst outlier regions
for(i in 1:17) #17 linkage groups
{
	cc = c1[c1[,2] == i, ]

	if(i == 1) plot(cc[cc[,3]>0,1],cc[cc[,3] > 0 ,3], type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), xaxt = "n", lwd = 4,col = ifelse(i %%2 == 0,colors()[153],colors()[230] ), main = main, xlab = "", ylab = "")  else 
	points(cc[cc[,3]>0,1],cc[cc[,3] > 0 ,3], type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), lwd  = 4, col =   ifelse(i %%2 == 0,"darkred","darkblue"))
	#gray coloring ### colors()[153],colors()[230] 
	#colours:  ifelse(i %%2 == 0,"darkred","darkblue")
	points(cc[cc[,16]<0.001,1],rep(1,length(cc[cc[,16]<0.001,1])), col = "black", pch = 19)

	axis(1, at = max_position[,3], label = c(1:17))
}
}
dev.print(device=svg, "~/Desktop/top3_4comp.svg", onefile=FALSE)
dev.off()




par(mfrow = c(3,2))
for(a in 1:6)
{
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

###fixed non fixed regions
for(i in 1:17) #17 linkage groups
{
	cc = c1[c1[,2] == i, ]
	#cc = cc[(cc[,4] != 0) & (cc[,5] != 0),]
	if(i == 1) plot(cc[,1],cc[,4]/ (cc[,4]+ cc[,5]), type = "h", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), xaxt = "n", lwd = 4,   col = ifelse(i %%2 == 0,"darkred", "darkblue"), main = main) else 
	points(cc[,1],cc[,4]/ (cc[,4]+ cc[,5]), type = "h", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), lwd = 4, col = ifelse(i %%2 == 0,"darkred", "darkblue"))
	
	points(cc[cc[,10]>0.999,1], rep(0.6,length(cc[cc[,10]>0.999,1])), col = "black", pch = 19)

	axis(1, at = max_position[,3], label = c(1:17))
}
}
 dev.print(device=svg, "graph_2/fixed.svg", onefile=FALSE)
dev.off()


par(mfrow = c(3,2))
for(a in 1:6)
{
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

###top1 regions
for(i in 1:17) #17 linkage groups
{
	cc = c1[c1[,2] == i, ]
	#cc = cc[(cc[,4] != 0) & (cc[,5] != 0),]
	if(i == 1) plot(cc[,1],cc[,6]/ (cc[,6]+ cc[,7]), type = "h", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), xaxt = "n", lwd = 4,   col = ifelse(i %%2 == 0,"darkred", "darkblue"), main = main) else 
points(cc[,1],cc[,6]/ (cc[,6]+ cc[,7]), type = "h", xlim = c(0,max(c1[,1])+1), ylim = c(-0.1,1), lwd = 4, col = ifelse(i %%2 == 0,"darkred", "darkblue"))
	
	points(cc[cc[,11]>0.999,1], rep(0.6,length(cc[cc[,11]>0.999,1])), col = "black", pch = 19)

	axis(1, at = max_position[,3], label = c(1:17))
}
}
dev.print(device=svg, "graph_2/top1.svg", onefile=FALSE)
dev.off()


################
###size based on FST###
################
c_ALL = NULL
size_all_all = NULL
ka = NULL
ks = NULL
dnds = NULL
#for(a in c(1,2,5,6))
for(a in c(1,2,3,4,5,6))
{
unique_map_transcript = read.delim("reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET";aa =1}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB";aa =2}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG";aa = 3}
if(a == 4) {c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB";aa =4}
if(a == 5) {c1 = read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG";aa =5}
if(a == 6) {c1 =  read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG";aa =6}

#for(z_c1 in 1:(nrow(c1)-1)) { if(round(c1[z_c1,1],4) == round(c1[(z_c1+1),1],4)) c1[(z_c1+1),2:22] = 0}


c_ALL = rbind(c_ALL,cbind(aa,c1,unique_map_transcript[,6]))
size = c(0,0,0,0,0,0,0) #LG, start, stop, how many in window
size_all = NULL 
what = 16
pvalue = 0.001
for(chr  in 1:17)
{
c1_temp = c1[c1[,2] == chr,]
c1_temp = rbind(c1_temp[1,],c1_temp)
for(i in 2:nrow(c1_temp))
{
if((i == 2) & (c1_temp[i,what] <= pvalue)) size[c(1,2,4)] = c1_temp[i,c(2,1,(what- 8))]#start of island 
if((c1_temp[i,what] <= pvalue) & (c1_temp[(i-1),what] > pvalue)) size[c(1,2,4)] = c1_temp[i,c(2,1,(what- 8))]#start of island 
if((c1_temp[i,what] >  pvalue) & (c1_temp[(i-1),what] <= pvalue)) {size[3] = c1_temp[(i-1),1]; size_all = rbind(size_all,as.numeric(size))} #end of island 
if( (nrow(c1_temp) == i) &  (c1_temp[i,what] <=  pvalue) & (c1_temp[(i-1),what] <=  0.001)) {size[3] = c1_temp[i,1]; size_all = rbind(size_all,as.numeric(size))} #end of island 
}
size = c(0,0,0,0,0,0,0) #LG, start, stop,, how many in window
}
size_all[,5] = (size_all[,3] - size_all[,2])

for(zi in 1:nrow(size_all)) {size_all[zi,7] = sum(c1[(c1[,1] >= size_all[zi,2]) & (c1[,1] <= size_all[zi,3]),12])}

#hist(size_all[,4], breaks = 50)
print(paste(nrow(size_all),round(mean(size_all[,5]),2), round(sd(size_all[,5]),2), round(max(size_all[,5]),2), main, "product", round(nrow(size_all) *mean(size_all[,4]),4),round(sum(size_all[,4])/   sum(c1[,(what-8)]),4)  ))
size_all_all = rbind(size_all_all, cbind(aa,size_all) )


dn = (c1[c1[,16] < pvalue,21])
dn_a = (c1[c1[,16] > pvalue,21])
ds = (c1[c1[,16] < pvalue,22])
ds_a = (c1[c1[,16] > pvalue,22])

dnds = rbind(dnds,c(signif(mean(dn_a, na.rm = T)/ mean(ds_a, na.rm = T),3),  signif( mean(dn, na.rm = T)/mean(ds, na.rm = T),3)))
ka = rbind(ka,c(signif(mean(dn_a,na.rm = T),3),  signif( mean(dn, na.rm = T),3)))
ks= rbind(ks,c(signif(mean(ds_a,na.rm = T),3),   signif(mean(ds, na.rm = T),3)))

colnames(size_all_all) = c("comp","LG","start","stop","SNP","size","recom_rate","number_genes")
print(sum(size_all_all[size_all_all[,1] == aa,5]) / sum(c1[,8]))

###recombination rates###
###recombination rates###
###recombination rates###

unique_map_transcript = read.delim("reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
for(i in 1:nrow(size_all_all))
			{
			size_all_all[i,7] = mean(unique_map_transcript[(unique_map_transcript[,4] >= size_all_all[i,3]) & (unique_map_transcript[,4] <= size_all_all[i,4]),6])
			}

#color = c("#99000025","#00009925","#00990025","#00000025","#00990025","#00000025")
color2 = c("lightcoral","limegreen ","aquamarine3 ","cyan2 ","dodgerblue1","orchid1 ")
color = c("#f8756b25","#7aae0025","#7aae0025","#00bdc225","#00bdc225","#c77aff25")

if(a == 1) {plot(c1[,3]~log(unique_map_transcript[,6]), col = color[a], xlim = c(-8,5.5), pch = 20, xlab = "recom_rate (cM/MB)",ylab = "Fst",xaxt = "n");
axis(1, at = log(c(0.001,0.01,0.1,1,10,100)),c(0.001,0.01,0.1,1,10,100))}
if((a == 2) | (a == 5) | (a == 6)) points(log(unique_map_transcript[,6]), c1[,3], col = color[a], xlim = c(-8,5.5),pch = 20)
 text(cex = 1.2,x = log(0.002), y =0.95-(aa*0.0375),font = 2, labels = paste("r^2 = ",signif(cor.test(c1[,3],log(unique_map_transcript[,6]))$estimate^2,4)," (",main,")",sep = ""),col = color2[aa])
if(a == 1) lm_regression_1 = lm(c1[,3]~log(unique_map_transcript[,6]))
if(a == 2) lm_regression_2 = lm(c1[,3]~log(unique_map_transcript[,6]))
if(a == 3) lm_regression_3 = lm(c1[,3]~log(unique_map_transcript[,6]))
if(a == 4) lm_regression_4 = lm(c1[,3]~log(unique_map_transcript[,6]))
if(a == 5) lm_regression_5 = lm(c1[,3]~log(unique_map_transcript[,6]))
if(a == 6) lm_regression_6 = lm(c1[,3]~log(unique_map_transcript[,6]))
#abline(lm(c1[,3]~log(unique_map_transcript[,6])), lwd = 6,col = color2[a])
}


abline(lm_regression_1, lwd = 6,col = color2[1])
abline(lm_regression_2, lwd = 6,col = color2[2])
#abline(lm_regression_3, lwd = 6,col = color2[3])
#abline(lm_regression_4, lwd = 6,col = color2[4])
abline(lm_regression_5, lwd = 6,col = color2[5])
abline(lm_regression_6, lwd = 6,col = color2[6])

recom_rate_summary = matrix(0,nrow = 13, ncol = 3)
rownames(recom_rate_summary) = c("isANN_PET","isANN_DEB","isANN_ARG","isPET_DEB","isPET_ARG","isDEB_ARG","null","ANN_PET","ANN_DEB","ANN_ARG","PET_DEB","PET_ARG","DEB_ARG")

colnames(recom_rate_summary) = c("median","CI_plus","CI_minus")

for(aa in c(1,2,3,4,5,6)){
#recom_rate_summary[2,aa] = signif(mean([(c_ALL[,1] == aa) & (ifelse(c_ALL[,16]<0.001,T,F)),24],na.rm = T),4) #mean recombination rate ISLANDS
recom_rate_summary[aa,1] = sign- if(median(c_ALL[(c_ALL[,1] == aa) & (ifelse(c_ALL[,17]<0.001,T,F)),24],na.rm = T),4) #median recombination rate ISLANDS
#recom_rate_summary[1,aa] = signif(mean(c_ALL[(c_ALL[,1] == aa) & (ifelse(c_ALL[,16]<0.001,F,T)),24],na.rm = T),4) #mean recombination rate ISLANDS
recom_rate_summary[aa+7,1] = signif(median(c_ALL[(c_ALL[,1] == aa) & (ifelse(c_ALL[,17]<0.001,F,T)),24],na.rm = T),4) #median recombination rate ISLANDS

r_ALL_is = (c_ALL[(c_ALL[,1] == aa) & (ifelse(c_ALL[,17]<0.001,T,F)),24]); r_ALL_is = r_ALL_is[!is.na(r_ALL_is)]
r_ALL_non = (c_ALL[(c_ALL[,1] == aa) & (ifelse(c_ALL[,17]<0.001,F,T)),24]); r_ALL_non = r_ALL_non[!is.na(r_ALL_non)]
recom_rate_summary[aa,2] = sort(r_ALL_is)[floor(length(r_ALL_is) * 0.5 + 1.96 * sqrt(length(r_ALL_is) * 0.5 * (1-0.5)))] #see https://epilab.ich.ucl.ac.uk/coursematerial/statistics/non_parametric/confidence_interval.html for median CI calculations
recom_rate_summary[aa+7,2] = sort(r_ALL_non)[floor(length(r_ALL_non) * 0.5 + 1.96 * sqrt(length(r_ALL_non) * 0.5 * (1-0.5)))]
recom_rate_summary[aa,3] = sort(r_ALL_is)[floor(length(r_ALL_is) * 0.5 - 1.96 * sqrt(length(r_ALL_is) * 0.5 * (1-0.5)))]
recom_rate_summary[aa+7,3] = sort(r_ALL_non)[floor(length(r_ALL_non) * 0.5 - 1.96 * sqrt(length(r_ALL_non) * 0.5 * (1-0.5)))]

}

###BARPLOTS WITH GGPLOT###
#confidence interval#
library(ggplot2)
recom_rate_summary = as.data.frame(recom_rate_summary[c(1:2,5:9,12:13),])
limits = aes(ymax = CI_plus, ymin = CI_minus)
p = ggplot(recom_rate_summary,aes(x= c(6,7,9,8,5,1,2,4,3),y = median) )
p + geom_bar(stat = "identity", aes(fill = factor(c(1,2,3,4,1,1,2,3,4) ))) + geom_errorbar(limits,width=0.25)
#

dev.print(device=svg, "~/Desktop/recombination_rate.svg", onefile=FALSE)
dev.off()

write.table(size_all_all,"results/size_all_all", row.names = F, col.names = T, quote = F)


###
###what is the recombination rate for the islands###
###
###ANOVA
size_all_all = size_all_all[size_all_all[,1] != 3,]; size_all_all = size_all_all[size_all_all[,1] != 4,] #remove comparison 3 and 4. 
size_all_all = as.data.frame(size_all_all); size_all_all[,1] = as.factor(size_all_all[,1])
fit_size_all_all <- lm(size_all_all[,6]~size_all_all[,1], data=size_all_all)
anova(fit_size_all_all)
kruskal.test(x = size_all_all[,6],g = as.factor(size_all_all[,1])) #island size
chisq.test(c(59,53,44,54)) #island number

#boxplot(size_all_all[,5]~size_all_all[,1], ylim = c(-2,5))

library(ggplot2)
size_all_all = as.data.frame(size_all_all)
p = ggplot(cbind(size_all_all[1],(size_all_all[,6])),aes(size_all_all[,1],(size_all_all[,6])))
p +  geom_boxplot(aes(colour = factor(size_all_all[,1])),outlier.size = 0) + geom_jitter(aes(colour = factor(size_all_all[,1])),position=position_jitter(width=0.1))


dev.print(device=svg, "~/Desktop/number_islands_top3_4comp.svg", onefile=FALSE)
dev.off()





###############
###HISTOGRAM###
###############
par(mfrow= c(2,2))
main = list(c("annuus-petiolaris (sympatric-ancient)"),c("annuus-debilis (parapatric-ancient)"),c("annuus-argophyllus (allopatric-recent)"),c("petiolaris-debilis (allopatric-recent)"), c("petiolaris-argophyllus (allopatric-ancient)"),c("debilis-argophyllus (allopatric-ancient)") )

for(a in c(1,2,6,5))
{
if(a == 1) {fst_map_match = read.delim("results/ann_pet_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "ANN_PET"} # ann_pet_fst
if(a == 2) {fst_map_match = read.delim("results/ann_deb_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "ANN_DEB"} # pet_deb_fst
if(a == 3) {fst_map_match = read.delim("results/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "ANN_ARG"} # ann_deb_fst
if(a == 4) {fst_map_match = read.delim("results/pet_deb_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "PET_DEB"} # ann_pet_fst
if(a == 5) {fst_map_match = read.delim("results/pet_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "PET_ARG"}# pet_deb_fst
if(a == 6) {fst_map_match = read.delim("results/deb_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); main = "DEB_ARG"} # ann_deb_fst

col = ncol(fst_map_match); m = mean(fst_map_match[,col-8], na.rm = T); print(m)

#x1 = hist(fst_map_match[,ncol(fst_map_match)-8], breaks = 100, plot = T, xlim = c(-0.1,1), main = main[1],col = "black") # ,freq = F)

x = fst_map_match[fst_map_match[,ncol(fst_map_match)-8] > 0.9,ncol(fst_map_match)-8] 

print(mean(fst_map_match[,ncol(fst_map_match)-8], na.rm = T))

xx = jitter(x,2000)
xx = xx[!is.na(xx)]
xx[xx > 1] = (1 / xx[xx > 1])
y = fst_map_match[fst_map_match[,ncol(fst_map_match)-8] <= 0.9,ncol(fst_map_match)-8] 

x1 = hist(c(xx,y), breaks =  seq(-0.6,1.1,by = 0.045), plot = T,freq = T, ylim = c(0,25000) , xlim= c(-0.1,1.1))
#x2 = hist(c(x,y), breaks = 100, plot =T,freq = F)

#plot(lowess(x1$counts~x1$mids , f=.08), col = "darkgrey", lwd = 5, type = "l",xlim = c(-0.05,1) , ylim = c(0,max(x2$counts)))
#points(x2$counts~x2$mids)

#plot(x1$counts/sum(x1$counts)~x1$mids ,type = "l", col = "darkblue", lwd = 5, xlab = list(expression(~italic(F)[ST]), cex = 1.8, font = 2), ylab = list("Frequency", cex = 1.8, font = 2), main = main[1], ylim = c(0,max(x1$counts/sum(x1$counts))), yaxt = "n", xlim = c(-0.05,1))
#axis(2,at =  c(0,2000,4000,6000,8000) /sum(x1$counts), labels =   c(0,2000,4000,6000,8000), font = 2)

}

dev.print(device=svg, "~/Desktop/histogram_4comp.svg", onefile=FALSE)
dev.off()


###PLOT GENETIC MAP###
###PLOT GENETIC MAP###
###PLOT GENETIC MAP###
library(qtl)

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new m
map = map[map[,3] != "no",]
write.table(cbind(c(1,"",""),t(cbind(map[,1],map[,2],map[,3]))), "reference/annuus_genetic_map.csv",row.names = F,col.names = F, quote = F, sep = ",")

genetic_map = read.cross("csv", file ="reference/annuus_genetic_map.csv" )
plot.map(genetic_map)

dev.print(device=svg, "~/Desktop/map.svg", onefile=FALSE)
dev.off()

###
### PLOT BOTH GENETIC & PHYSICAL MAP###
### 

library(qtl)
unique_map_transcript_all = read.delim("reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
unique_map_transcript_all[is.na(unique_map_transcript_all[,6]),6] =  1
dup = duplicated(unique_map_transcript_all[,4]); unique_map_transcript_all = unique_map_transcript_all[dup == F, ]

write.table(cbind(c(1,"",""),t(cbind(unique_map_transcript_all[,1],unique_map_transcript_all[,2],unique_map_transcript_all[,3]))), "results/annuus_genetic_map.csv",row.names = F,col.names = F, quote = F, sep = ",")
write.table(cbind(c(1,"",""),t(cbind(unique_map_transcript_all[,1],unique_map_transcript_all[,2],unique_map_transcript_all[,7]))), "results/annuus_physical_map.csv",row.names = F,col.names = F, quote = F, sep = ",")
genetic_map = read.cross("csv", file ="results/annuus_genetic_map.csv" )
physical_map = read.cross("csv", file ="results/annuus_physical_map.csv" )

plot.map(genetic_map,physical_map) # ,  chr = c(1,2))
plot.map(genetic_map, chr = c(1,2))
plot.map(physical_map, chr = c(1,2))

dev.print(device=png, "~/Desktop/physical_map.png", onefile=FALSE)

dev.off()


###
### PLOT BOTH GENETIC  MAPs Grassa and Bowers' map###
### 
setwd("~/Documents/speciation_islands_dec2011") #set up working directory 
map_g2 = read.delim("results/genetic_map_comparisons", header = T, stringsAsFactors = F,sep = "\t") ###recombination rates###

map_g2 = map_g2[map_+-g2[,2] != 0, ] 

plot(map_g2,xaxt = "n", yaxt = "n", pch = 20, lwd = 0.5, cex = 0.5) # , ylim = c(1.5,18.1), xlim = c(1.5,18.1))
axis(1, c(1:17));axis(2, c(1:17))
for(i in 1:18)
{
lines(y = c(i,i),x = c(1,18), col = "#00000050" )
lines(y = c(1,18),x = c(i,i),col = "#00000050" )
}
dev.print(device=svg, "~/Desktop/genetic_map_comparisons.svg", onefile=FALSE)
dev.off()

###SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX###
###SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX###
###SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX######SANDBOX###

x = fst_map_match[fst_map_match[,ncol(fst_map_match)-8] > 0.9,ncol(fst_map_match)-8] 
xx = jitter(x,1000)
xx[xx > 1] = 1
y = fst_map_match[fst_map_match[,ncol(fst_map_match)-8] <= 0.9,ncol(fst_map_match)-8] 

x1 = hist(c(xx,y), breaks = 100)
x2 = hist(c(x,y), breaks = 100)

plot(x2$counts/sum(x2$counts)~x2$mids )

lines(lowess(x1$counts/sum(x1$counts)~x1$mids , f=.05), col = 2, lwd = 5)

#for(i in 1:nrow(fst_map_match))
for(i in 1:10000)
{
if(fst_map_match[i,   ncol(fst_map_match)-8] > 0.9]     ) fst_map_match[i,ncol(fst_map_match)-8] =  jitter(fst_map_match[i,ncol(fst_map_match)-8], 10)

if(i %% 1000 == 0) print(i)
}


x1 = hist(fst_map_match[,ncol(fst_map_match)-8], breaks = 100, plot = T, xlim = c(-0.1,1), main = main[1],col = "black") # ,freq = F)
plot(x1$counts/sum(x1$counts)~x1$mids)

fst_map_match[,(col-8)] = jitter(fst_map_match[,(col-8)],0.0000001)


c_4comp = list(c1,c2,c3)

main = list(c("annuus-petiolaris(sympatric-ancient)"),c("annuus-debilis(parapatric-ancient)"),c("annuus-argophyllus(allopatric-recent)"),c("petiolaris-debilis (allopatric-recent)"),c("petiolaris-argophyllus(allopatric-ancient)"),c("debilis-argophyllus (allopatric-ancient)") )

par(mfrow = c(3,2))
for(b in 1:6)
       {
               ####CLUSTERING####
               #plot qvalues of chi square test
               plot(x = c_4comp[[b]][,1], y = c_4comp[[b]][,9], main = main[[b]],
lwd = 3, xlab = "Linkage Group", ylab ="Fisher Exact Test (-log(p-value))" ,
               col = "darkblue", type = "h", xaxt = "n", ylim = c(0,10))

               axis(1, at = max_position[,3], label = c(1:17))
               abline(h = -log(0.001,10), lwd = 2, lty = 2)
       #       if(b != 4) points(temp_4[,1],rep(max(c_4comp[[b]][,9]),nrow(temp_4)),type = "h", lwd = 7, col = "#99000050")
       }
dev.print(device=svg, "graphs/clustering_top1percent.svg", onefile=FALSE)
dev.off()

###########################################
### EFFECT OF GENE DENSITY - XY SCATTER ###
###########################################

par(mfrow = c(3,2))

for(b in 1:6)
{		
		z =  data.frame(cbind(c_4comp[[b]][c_4comp[[b]][,9] != 0 ,12], c_4comp[[b]][c_4comp[[b]][,9] != 0,9]))
		
		plot(z, pch = 10, lwd = 3, col = densCols(z[,1],colramp = colorRampPalette(c("red","black"))),
		ylab = "Fisher Exact Test (-log(Qvalue))", xlab = "Number of Genes in 1cM Sliding Window", main = main[b]) # number of genes in SW. 
		
		lm.r = lm(z[,2]~ z[,1])
		print(summary(lm.r))
		y_text = c(100,50,30,20,20,200)
		x_text = c(400,400,400,400,400,400)
		text = c("R-squared = 0.61","R-squared = 0.54", "R-squared = 0.28","R-squared = 0.47","R-squared = 073", "R-squared = 0.62")
	#	text_pval = c("P-value < 0.0001","P-value < 0.0001","P-value < 0.0001","n.s.")
		abline(lm.r, col = "black", lwd = 5)
		text(x_text[b], y_text[b], text[b])
	#	text(x_text[b], y_pval[b], text_pval[b])
	}	
dev.print(device=svg, "graphs/speciation_genedensity_effect1.svg", onefile=FALSE)
dev.off()

#########################################
### EFFECT OF GENE DENSITY - ALONG LG ###
#########################################

par(mfrow = c(3,2))
for(b in 1:6)
	{
		####CLUSTERING####
		#plot qvalues of chi square test
		plot(x = c_4comp[[b]][,1], y = (c_4comp[[b]][,5]/(c_4comp[[b]][,5]+c_4comp[[b]][,6])), main = main[[b]], lwd = 3, xlab = "Linkage Group", ylab ="gene density | proportion fixed",
		ylim = c(-1,1), col = "darkred", type = "h",xaxt = "n", yaxt = "n")
		lines(x = c_4comp[[b]][,1], y = -(c_4comp[[b]][,12]/max(c_4comp[[b]][,12]) *1), main = main[[b]], lwd = 3, col = "black")
		#lines(x = c_4comp[[b]][,1], y = c_4comp[[b]][,9], main = main[[b]], lwd = 1, xlab = "Linkage Group", ylab ="Fisher Exact Test (-log(Qvalue))" , col = "darkblue", type = "h", xaxt = "n")
		
		axis(2, at = c(-(809/max(c_4comp[[b]][,12]) *1),-(600/max(c_4comp[[b]][,12]) *1),-(400/max(c_4comp[[b]][,12]) *1),-(200/max(c_4comp[[b]][,12]) *1),0), label = c(800,600,400,200,0))
		
		max = max(c_4comp[[b]][,9])
		axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels =   c(0,0.2,0.4,0.6,0.8,1))
		axis(1, at = max_position[,3], label = c(1:17))
		abline(h = -log(0.0001,10), lwd = 3, lty = 2)
		points(c_4comp[[b]][(c_4comp[[b]][,9] > 3),1],rep(1,length(c_4comp[[b]][(c_4comp[[b]][,9] > 3),1])), col = "darkblue", pch = 20, lwd = 3)

#		points(as.numeric(comp5[[b]][as.numeric(comp5[[b]][,ncol(comp5[[b]])-6]) > 0.9,ncol(comp5[[b]])]), as.numeric(comp5[[b]][as.numeric(comp5[[b]][,ncol(comp5[[b]])-6]) > 0.9,ncol(comp5[[b]])-6])* 50, col = "#99000050")
	}	
dev.print(device=svg, "graphs/speciation_genedensity_effect_v3.svg", onefile=FALSE)
dev.off()


p = ggplot(size_all_all[c(5,1)],aes(size_all_all[,1],size_all_all[,5])); p +  geom_boxplot(aes(colour = factor(size_all_all[,1]))) + geom_point(aes(colour = factor(size_all_all[,1])))


size_all_all = as.data.frame(size_all_all)

p = ggplot(cbind(size_all_all[1],log(size_all_all[,5])),aes(size_all_all[,1],log(size_all_all[,5])))
p +  geom_boxplot(aes(colour = factor(size_all_all[,1]))) + geom_point(aes(colour = factor(size_all_all[,1])))


boxplot()


####reproducibility LOREN MARCH 2012

par(mfrow = c(2,1))
for(a in 3:4)
{
################
### PLOTTING ###
################
#if(a == 7) c1 =  read.table("results_6species/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {c2 =  read.table("results/pet_deb_fst_map_cluster",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {c1 =  read.table("results/pet_arg_fst_map_cluster",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {c1 = read.table("results/deb_arg_fst_map_cluster",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

###fst outlier regions
#for(i in 1:17) #17 linkage groups
for(i in 3:3)
{
	cc = c1[c1[,2] == i, ]

	cc[,1] = cc[,1] - cc[1,1]


	plot(cc[cc[,3]>0,1],cc[cc[,3] > 0 ,3], type = "l",ylim = c(0,0.8), lwd = 4,col = ifelse(i %%2 == 0,"darkred", "darkblue"), main = main, xlab = "", ylab = "")  
	axis(1, at = cc[,1], label = rep("",nrow(cc)))
}
}
cc[,1] = cc[,1] - cc[1,1]
i = 1
cc = c3[c3[,2] == i, ]
cc2 = c4[c4[,2] == i, ]

plot(cc[,3],cc2[,3])

axis(1, at = cc[,1], label = rep("",nrow(cc)))
axis(2, at = cc[,1], label = rep("",nrow(cc)))

###plotting the physical and genetic  map (OLD VERSION DOESNT REALLY WORK)
dup = duplicated(map[,4]);length(dup[dup == F])

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F);map = map[map[,3] != "no",] # new m
map = map[dup == F,] # remove all duplicate positions

map_temp = cbind(map,unique_map_transcript[,6],0,0,0)
map_conversion =  NULL

for(j in 1:17)
{
	map_chr = map_temp[map_temp[,2] == j, ]
	
		for(i in 1:nrow(map_chr))
		{
			if(i == 1) map_chr[i,5] = as.numeric(map_chr[(i+1),3]) - as.numeric(map_chr[i,3])
			if(i == nrow(map_chr)) map_chr[i,5] = 0 # convert to distance between 2 positions
		}
		map_chr[,6] = as.numeric(map_chr[,5]) * (1/map_chr[,4])  #convert to megabases
		
		for(i in 1:nrow(map_chr))#convert back to cumulative distances, this time in megabases.
		{
			if(i == 1) map_chr[i,7] =  map_chr[i,3] 
			if(i > 1) map_chr[i,7] = as.numeric(map_chr[(i-1),7]) + as.numeric(map_chr[i,6])  #convert back to cumulative distances, this time in megabases.
		}
		map_conversion = rbind(map_conversion,map_chr)
}


write.table(cbind(c(1,"",""),t(cbind(map[,1],map[,2],map[,3]))), "~/Desktop/annuus_genetic_map.csv",row.names = F,col.names = F, quote = F, sep = ",")

write.table(cbind(c(1,"",""),t(cbind(map[,1],map[,2],as.numeric(map[,3])*unique_map_transcript[,6] / mean(unique_map_transcript[,6])))), "~/Desktop/annuus_physical_map.csv",row.names = F,col.names = F, quote = F, sep = ",") ###original converter###

write.table(cbind(c(1,"",""),t(cbind(map_conversion[,1],map_conversion[,2],as.numeric(map_conversion[,7]) ))), "~/Desktop/annuus_physical_map.csv",row.names = F,col.names = F, quote = F, sep = ",")

write.table(cbind(c(1,"",""),t(cbind(map[,1],map[,2],as.numeric(map[,3])* unique_map_transcript[,6]))), "~/Desktop/annuus_physical_map.csv",row.names = F,col.names = F, quote = F, sep = ",") ###original converter###

genetic_map = read.cross("csv", file ="~/Desktop/annuus_genetic_map.csv" )
#physical_map = read.cross("csv", file ="results/annuus_physical_map.csv" )

physical_map = read.cross("csv", file ="~/Desktop/annuus_physical_map.csv" )

plot.map(genetic_map,physical_map)
plot.map(genetic_map)
plot.map(physical_map)

