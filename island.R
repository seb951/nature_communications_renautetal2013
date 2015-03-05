#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

for(a in from:to)
{
###clustering
if(a == 1)  {comp4 = read.table("results/ann_pet_comp4", stringsAsFactors =  F, header = T); bayes_limit = 0.9271889; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/ann_pet/dnds_ann_pet", sep = " ", stringsAsFactor = F)}
if(a == 2)  {comp4 =  read.table("results/ann_deb_comp4", stringsAsFactors =  F, header = T); bayes_limit = 0.980655; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/ann_deb/dnds_ann_deb", sep = " ", stringsAsFactor = F)}
if(a == 3)  {comp4 =  read.table("results/ann_arg_comp4", stringsAsFactors =  F, header = T); bayes_limit = 0.91788; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/ann_arg/dnds_ann_arg", sep = " ", stringsAsFactor = F)}
if(a == 4)  {comp4 = read.table("results/pet_deb_comp4", stringsAsFactors =  F, header = T); bayes_limit = 0.9056857; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/pet_deb/dnds_pet_deb", sep = " ", stringsAsFactor = F)}
if(a == 5)  {comp4 =  read.table("results/pet_arg_comp4", stringsAsFactors =  F, header = T); bayes_limit = 1; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/pet_arg/dnds_pet_arg", sep = " ", stringsAsFactor = F)}
if(a == 6)  {comp4 =  read.table("results/deb_arg_comp4", stringsAsFactors =  F, header = T); bayes_limit = 1; dnds = read.delim("~/PAML/ann_deb_pet_arg_RESULTS/deb_arg/dnds_deb_arg", sep = " ", stringsAsFactor = F)}


####################
### MAPPING FST VALUES ###
####################

map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",]
map = cbind(map,0,0,0,0)
colnames(map) = c("name","LG","map_centiMorgan","unique_position","nbnuc","dn","ds")

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
	if(i %% 1000 == 0) print(paste(i,Sys.time()))
}
max_position[,3] = (max_position[,3] +1)

#add dnds
for(i in 1:nrow(map))
#for(i in 1:200)
{
x = dnds[,1] %in% map[i,1]

if(nrow(dnds[x == T,]) == 1)    map[i,6:7] = dnds[x == T,2:3]

if(i %% 1000 == 0) print(i)
}

###
###
###
#comp4 = comp4[1:2012,]
###To keep information seperate for each SNP###
	comp5 =  cbind(comp4,0,0,0,0,0)
	colnames(comp5) = c(colnames(comp4),colnames(map)[1:5])
	map_temp = map[1:100,1:5]; map_temp[,2:5] = 0 
	ceiling = ceiling(nrow(comp5)/100)
	map_all = NULL
	for(j in 1:ceiling)
		{
			if(j != ceiling) {map_temp = map[1:100,1:5]; map_temp[,2:5] = 0} else {map_temp = map[1:(nrow(comp5) %% 100),1:5]; map_temp[,2:5] = 0}
			if(length(map_temp) == 0)map_temp = map[1:100,1:5]; map_temp[,2:5] = 0 
			for(z in 1:nrow(map_temp))
			{
			x = map[,1] %in% comp5[(((j -1)*100)+ z),1]
			if(length(x[x == T]) == 1)  map_temp[z,] = map[(x == T),1:5]
			if(length(x[x == T]) > 1) print(map[(x == T),1:5])
			}			
			map_all = rbind(map_all,map_temp)
			if(j %% 100 == 0) print(paste(j,"of ceiling",ceiling,Sys.time()))		
			}		
			comp5[,(ncol(comp4)+1):(ncol(comp4)+5)] = map_all
	

###add a column to flag whether the SNP is in the top 1% quantile.###

	col = ncol(comp5) 
	top_quantile = quantile(as.numeric(comp5[,col-7]), 0.99,na.rm = T)
	comp5 = cbind(comp5,0)
	colnames(comp5)[col+1] = "top-quant"
	for(i in 1:nrow(comp5))
		{if((as.numeric(comp5[i,col-7]) > top_quantile)   & !is.na(comp5[i,col-7])  ) comp5[i,col +1] = 1}

if(a == 1) write.table(comp5,"results/ann_pet_fst_map", row.names = F, col.names = T)
if(a == 2) write.table(comp5,"results/ann_deb_fst_map", row.names = F, col.names = T)
if(a == 3) write.table(comp5,"results/ann_arg_fst_map", row.names = F, col.names = T)
if(a == 4) write.table(comp5,"results/pet_deb_fst_map", row.names = F, col.names = T)
if(a == 5) write.table(comp5,"results/pet_arg_fst_map", row.names = F, col.names = T)
if(a == 6) write.table(comp5,"results/deb_arg_fst_map", row.names = F, col.names = T)


###################
#### CLUSTERING ###
###################
gmean=function(y=y,x=x, q=q,d=d){  # y=values, x=snp pos vector, q=query pos, d=dist:    Rose's gaussian average
	weights=rep(NA,length(x))
	weights[abs(x-q)<3*d]=exp(-(x[abs(x-q)<3*d]-q)^2/(2*d^2))
	weighted=weights*y/sum(weights,na.rm=T)
	return(sum(weighted,na.rm=T))
}

if(a == 1) fst_map_match = read.delim("results/ann_pet_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(a == 2) fst_map_match = read.delim("results/ann_deb_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(a == 3) fst_map_match = read.delim("results/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst
if(a == 4) fst_map_match = read.delim("results/pet_deb_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(a == 5) fst_map_match = read.delim("results/pet_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(a == 6) fst_map_match = read.delim("results/deb_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst
	
	fst_map_match = cbind(comp4[,6:7],fst_map_match)
	col = ncol(fst_map_match) 
	fst_map_match[is.na(fst_map_match[,(col-8)]),(col-8):(col-6)] = 0
	fst_map_match[,(col-8)] = jitter(fst_map_match[,(col-8)],0.1)  #add/remove between -/+ 2e-05. to all fst values 

	top1 = quantile(fst_map_match[,(col-8)],0.99,na.rm = T)
	top3 = quantile(fst_map_match[,(col-8)],0.97,na.rm = T)
	top5 = quantile(fst_map_match[,(col-8)],0.95,na.rm = T)
	####running window - distance based###
 	d = 0.5 #distance cutoff in cM
	
	c = cbind(unique(as.numeric(map[,4])),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #position you are interrogating along with the result of the distance based running mean. 
	dnds_temp = cbind(unique(as.numeric(map[,4])),0,0) #temporary file with dnds
	for(i in 1:nrow(c))
		{
				temp = map[as.numeric(map[,4]) == c[i,1],]
				c[i,2] = temp[,2][1]
				c[i,18] = sum(as.numeric(temp[,5])) #number nuc
				if(!is.na(mean(as.numeric(temp[as.numeric(temp[,6]) > 0 & as.numeric(temp[,7] > 0),6]),na.rm = T))) dnds_temp[i,2] = signif(mean(as.numeric(temp[as.numeric(temp[,6]) > 0 & as.numeric(temp[,7] > 0),6]),na.rm = T),4) #dn
				if(!is.na(mean(as.numeric(temp[as.numeric(temp[,6]) > 0 & as.numeric(temp[,7] > 0),7]),na.rm = T))) dnds_temp[i,3] = signif(mean(as.numeric(temp[as.numeric(temp[,7]) > 0 & as.numeric(temp[,7] > 0),7]),na.rm = T),4) #ds
					dnds_temp
		}

	c = matrix(as.numeric(c), ncol = 22)
	colnames(c) = c("unique_position","LG","mean_sliding_fst", "fixed_window","nonfixed_window", "bayes_window","nonbayes_window", "top3_window","nontop3_window","top5_window","nontop5_window","number_genes_unique_position","resampling_mean","resampling_fixed","resampling_bayes","resampling_top3","resampling_top5","nbnuc","dn","ds","ns","s") #,"fst","resampling_fst")

	for(chr in 1:17)
#	for(chr in 1:1)
		{
		chromo_temp = fst_map_match[(as.numeric(fst_map_match[,col-4]) == chr),] # subset of chromosome a
		
		for(i in c(1:nrow(c))[c[,2] == chr]) 
			{
			temp = chromo_temp[chromo_temp[,(col-2)] > (c[i,1]-d) & chromo_temp[,(col-2)] < (c[i,1]+d),] # mean fst at unique position in a sliding window of X cM
			c_temp = dnds_temp[c[,1] > (c[i,1]-d) & c[,1] < (c[i,1]+d),2:3] # mean dnds in a sliding window approach
			if(nrow(temp) >0) c[i,3] = gmean(temp[,(col-8)], temp[,(col-2)],c[i,1],d)	 #gaussian mean
			if(nrow(temp) >0) c[i,4] = length(temp[temp[,(col-8)] >= 0.95,(col-8)]) # observed fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,5] = length(temp[temp[,(col-8)] < 0.95,(col-8)])  # observed non fixed SNP per sliding window...

			if(nrow(temp) >0) c[i,6] = length(temp[temp[,(col-8)] >= bayes_limit,(col-8)]) # observed fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,7] = length(temp[temp[,(col-8)] < bayes_limit,(col-8)])  # observed non fixed SNP per sliding window...

			if(nrow(temp) >0) c[i,8] = length(temp[temp[,(col-8)] >= top3,(col-8)]) # observed fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,9] = length(temp[temp[,(col-8)] < top3,(col-8)])  # observed non fixed SNP per sliding window...

			if(nrow(temp) >0) c[i,10] = length(temp[temp[,(col-8)] >= top5,(col-8)]) # observed fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,11] = length(temp[temp[,(col-8)] < top5,(col-8)])  # observed non fixed SNP per sliding window...

			if(nrow(temp) >0) c[i,21] = length(temp[temp[,2] == 1,(col-8)]) # observed ns mutations per sliding window...
			if(nrow(temp) >0) c[i,22] = length(temp[temp[,1] == 1,(col-8)])  # observed s mutations per sliding window...
		
			if(nrow(temp) >0) c[i,12] = length(unique(sort(fst_map_match[fst_map_match[,col-2] == c[i,1],col-5]))) #genes per unique position (not sliding)# 
			
			if(length(c_temp) >2) c[i,19]  = mean(c_temp[c_temp[,1] > 0,1]) else c[i,19]  = c_temp[1] 
			if(length(c_temp) >2) c[i,20]  = mean(c_temp[c_temp[,2] > 0,2]) else c[i,20]  = c_temp[2] 

			#chromo_temp_fst =  chromo_temp[chromo_temp[,(col-2)] == round(c[i,1],5),(col-8)] 
			#if(length(chromo_temp_fst) > 0) c[i,13] = mean(chromo_temp_fst, na.rm = T)	 #mean at position c[i,1] #not sliding#
	
			perm_100 = c(0,0,0,0,0) # mean, fixed, bayes, top3, top5
			for(s in 1:1000) # 1000 resampling
				{			
				temp_perm = sample(fst_map_match[,(col-8)],nrow(temp),replace = T)
				if((nrow(temp) >0) & (c[i,3] > 0) & (c[i,3] > mean(temp_perm) ))	 perm_100[1] = perm_100[1] + 1
				if((nrow(temp) >0) &  (c[i,4] > 0) & (c[i,4] >= length(temp_perm[temp_perm > 0.95 ])))	 perm_100[2] = perm_100[2] + 1
				if((nrow(temp) >0) &  (c[i,6] > 0) & (c[i,6] >= length(temp_perm[temp_perm >= bayes_limit ])))	 perm_100[3]= perm_100[3] + 1
				if((nrow(temp) >0) & (c[i,8] > 0) &  (c[i,8] >= length(temp_perm[temp_perm > top3 ])))	 perm_100[4]= perm_100[4] + 1
				if((nrow(temp) >0) &  (c[i,10] > 0) & (c[i,10] >= length(temp_perm[temp_perm > top5 ])))	 perm_100[5]= perm_100[5] + 1
			#	if((length(chromo_temp_fst) > 0) & (c[i,13] >= mean(temp_perm_fst, na.rm = T)))	 perm_100[4]= perm_100[4] + 1
				}

				if(perm_100[1] > 900 | perm_100[2] > 900 | perm_100[3] > 900 | perm_100[4] > 900  | perm_100[5] > 900) 	for(s in 1:99000) # 10000 resampling
				{			
				temp_perm = sample(fst_map_match[,(col-8)],nrow(temp),replace = T)
				#temp_perm_fst = sample(fst_map_match[,(col-8)],length(chromo_temp_fst)) #mean not sliding
				if((perm_100[1] > 900) & (c[i,3] > 0) &  (c[i,3] >	mean(temp_perm)))	 perm_100[1] = perm_100[1] + 1
				if((perm_100[2] > 900) &  (c[i,4] > 0) & (c[i,4] >= length(temp_perm[temp_perm > 0.95 ]))     )	 perm_100[2] = perm_100[2] + 1
				if((perm_100[3] > 900) & (c[i,6] > 0) &  (c[i,6] >= length(temp_perm[temp_perm >= bayes_limit ])))	 perm_100[3] = perm_100[3] + 1
				if((perm_100[4] > 900) & (c[i,8] > 0) &  (c[i,8] >= length(temp_perm[temp_perm > top3 ])))	 perm_100[4] = perm_100[4] + 1
				if((perm_100[5] > 900) & (c[i,10] > 0) &  (c[i,10] >= length(temp_perm[temp_perm > top5 ])))	 perm_100[5] = perm_100[5] + 1
			#	if((length(chromo_temp_fst) > 0) & (c[i,13] >=  mean(temp_perm_fst, na.rm = T)))	 perm_100[4]= perm_100[4] + 1
				} 

			if(perm_100[1] < 1000)   c[i,13] = 1 - (perm_100[1] / 1000) else c[i,13] = 1 - (perm_100[1] / 100000)
			if(perm_100[2] < 1000)   c[i,14] = 1 - (perm_100[2] / 1000) else c[i,14] = 1 - (perm_100[2] / 100000)
			if(perm_100[3] < 1000)   c[i,15] = 1 - (perm_100[3] / 1000) else c[i,15] = 1 - (perm_100[3] / 100000)
			if(perm_100[4] < 1000)   c[i,16] = 1 - (perm_100[4] / 1000) else c[i,16] = 1 - (perm_100[4] / 100000)
			if(perm_100[5] < 1000)   c[i,17] = 1 - (perm_100[5] / 1000) else c[i,17] = 1 - (perm_100[5] / 100000)
		#	if(perm_100[4] < 1000)   c[i,14] = perm_100[4] / 1000 else c[i,14] = perm_100[4] / 10000
			perm_100 = NULL
			}

		print(paste(chr, Sys.time()))
}

if(a == 1) write.table(c,"results/ann_pet_fst_map_cluster2", row.names = F)
if(a == 2) write.table(c,"results/ann_deb_fst_map_cluster2", row.names = F)
if(a == 3) write.table(c,"results/ann_arg_fst_map_cluster2_nobtm19", row.names = F)
if(a == 4) write.table(c,"results/pet_deb_fst_map_cluster2", row.names = F)
if(a == 5) write.table(c,"results/pet_arg_fst_map_cluster2_nobtm19", row.names = F)
if(a == 6) write.table(c,"results/deb_arg_fst_map_cluster2_nobtm19", row.names = F)


  }

