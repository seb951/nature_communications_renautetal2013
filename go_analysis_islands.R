#!/usr/bin/Rscript --verbose


#args = commandArgs(TRUE)
#a = as.numeric(args[1])

library(qvalue)
#############
### GO analysis###	
#############
setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

size_all_all = read.table("results/size_all_all",  header = T, stringsAsFactors = F)
for(a in 1:6)
{
if(a == 1) {fst_map_match = read.delim("results/ann_pet_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {fst_map_match = read.delim("results/ann_deb_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {fst_map_match = read.delim("results/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {fst_map_match = read.delim("results/pet_deb_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {fst_map_match = read.delim("results/pet_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {fst_map_match = read.delim("results/deb_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

col = ncol(fst_map_match)
	
###reference dataset
description_go = as.matrix(read.delim("reference/go_ID_description.txt", header = T))
go_annotation = as.matrix(read.delim("~/Documents/speciation_islands_dec2011/reference/go_annotations/go_annot", header = T,sep = ";"))
go_annotation[,4] = gsub("c_,","",go_annotation[,4], fixed = T); go_annotation[,4] = gsub("p_,","",go_annotation[,4], fixed = T);go_annotation[,4] = gsub(",f_$","",go_annotation[,4])
go_annotation[,4] = gsub("c_NA,","",go_annotation[,4], fixed = T); go_annotation[,4] = gsub("p_NA,","",go_annotation[,4], fixed = T);go_annotation[,4] = gsub("f_NA$","",go_annotation[,4])

go_annotation[,1] =  gsub(">","",go_annotation[,1], fixed = T)
go_representation = list("cc","bp","mf")

### only genes in dataset###	
ExpressedGenes = unique(fst_map_match[,1])
go_expressed_temp = rep(0,nrow(go_annotation))
	for(i in 1:nrow(go_annotation))
		{
		x = ExpressedGenes %in% go_annotation[i,1]
		if(length(x[x == T]) == 1) go_expressed_temp[i] = 1
			}
	go_annotation_expr = go_annotation[go_expressed_temp == 1,]

###GENOME WIDE DISTRIBUTION###		
	cc_go = strsplit(paste(go_annotation_expr[,4], collapse = ","),",")[[1]];cc_go = cc_go[regexpr("c_",cc_go) > 0];cc_go = cbind(rle(sort(cc_go))$lengths, rle(sort(cc_go))$values)
	mf_go = strsplit(paste(go_annotation_expr[,4], collapse = ","),",")[[1]];mf_go = mf_go[regexpr("f_",mf_go ) > 0];mf_go = cbind(rle(sort(mf_go))$lengths, rle(sort(mf_go))$values)
	bp_go = strsplit(paste(go_annotation_expr[,4], collapse = ","),",")[[1]];bp_go = bp_go[regexpr("p_",bp_go ) > 0];bp_go = cbind(rle(sort(bp_go))$lengths, rle(sort(bp_go))$values)
	
	go_all= list(bp_go,mf_go,cc_go)
	go_category = c("p_","f_","c_")
	
###GENOMIC POSITION PER GENOMIC POSITION###	
c1_go = data.frame(c1[,1],c1[,2],"A","A",0,"A","A",0,"A","A",0, stringsAsFactors = F)
c1_go[,1:2] = c1[,1:2]
colnames(c1_go) = c("unique_position","LG","GO_bp","GO_name_bp","GO_number_bp","GO_mf","GO_name_mf","GO_number_mf","GO_cc","GO_name_cc","GO_number_cc")

for(c in 1:nrow(c1))
#for(c in 1:500)
{
#for island one find out what genes are in there. 
w = 0.5 #window size

island_genes = fst_map_match[fst_map_match[,(col - 2)] >= (c1[c,1] -w),]; island_genes = island_genes[island_genes[,(col-2)] <= (c1[c,1]+w),]; island_genes = unique(island_genes[,1])

island_genes_GO = NULL

for(j in 1:length(island_genes))
{
x = go_annotation_expr[,1] %in% island_genes[j]
if(length(x[x == T]) ==1) island_genes_GO = paste(island_genes_GO,go_annotation_expr[x == T,4], collapse = ",", sep = ",")
}

	for(go in 1:3) #1 is bp, 2 is mf and three is cc
	{
	go_islands = strsplit(paste(island_genes_GO, collapse = ","),",")[[1]];go_islands = go_islands[regexpr(go_category[go],go_islands) > 0];go_islands = cbind(rle(sort(go_islands))$lengths, rle(sort(go_islands))$values,0,0,0,0,0)
	if(length(go_islands) == 5) go_islands = matrix(0,nrow= 2, ncol = 7)
	colnames(go_islands) = c("number","GO","proportion_observed","expected","pval","adjust_pvalBH","GO_name")

	go_islands[,3] = as.numeric(go_islands[,1]) / sum(as.numeric(go_islands[,1]))

	###FET adjusted
	for(go_is in 1:nrow(go_islands))
	{
	a = as.numeric(go_islands[go_is,1])
	b = sum(as.numeric(go_islands[,1]))
	cc = as.numeric(go_all[[go]][go_all[[go]][,2] == go_islands[go_is,2],1])
	d = sum(as.numeric(go_all[[go]][,1]))
	go_islands[go_is,5] = signif(fisher.test(cbind(c(a,b),c(cc,d)), alternative = "g")$p.value,2)
	if((length(cc) != 0) & length((d) != 0)) go_islands[go_is,4] = cc / d
	}
	go_islands[,6] =  p.adjust(as.numeric(go_islands[,5]),"fdr")

	###get the description
	for(des in 1:nrow(go_islands))
		{
		x = description_go[,4] %in% gsub("^..","",go_islands[des,2])
		if(length(x[x == T]) ==1) go_islands[des,7] = description_go[x == T,2]
		}

threshold = 0.01

if((go == 1) & (length(go_islands) > 7)) {c1_go[c,3] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,2],sep = "", collapse = ","); c1_go[c,4] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,7],sep = "", collapse = ",");c1_go[c,5] = length(go_islands[as.numeric(go_islands[,6])< threshold,2])}
if((go == 1) & (length(go_islands) == 7) & (as.numeric(go_islands[6])< threshold)) {c1_go[c,3:4] = go_islands[c(2,7)]; c1_go[c,5] = 1}

if((go == 2) & (length(go_islands) > 7)) {c1_go[c,6] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,2],sep = "", collapse = ","); c1_go[c,7] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,7],sep = "", collapse = ",");c1_go[c,8] = length(go_islands[as.numeric(go_islands[,6])< threshold,2])}
if((go == 2) & (length(go_islands) == 7) & (as.numeric(go_islands[6])<threshold)) {c1_go[c,6:7] = go_islands[c(2,7)]; c1_go[c,8] = 1}

if((go == 3) & (length(go_islands) > 7)) {c1_go[c,9] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,2],sep = "", collapse = ","); c1_go[c,10] =  paste(go_islands[as.numeric(go_islands[,6])< threshold,7],sep = "", collapse = ",");c1_go[c,11] = length(go_islands[as.numeric(go_islands[,6])< threshold,2])}
if((go == 3) & (length(go_islands) == 7) & (as.numeric(go_islands[6])< threshold)) {c1_go[c,9:10] = go_islands[c(2,7)]; c1_go[c,11] = 1}

if(c %% 1000 == 0) print(paste(c, "of 3047",Sys.time()))
}
}
print(rbind(c(mean(c1_go[c1[,16] <0.001,5]),mean(c1_go[c1[,16] >0.001,5])),c(mean(c1_go[c1[,16] <0.001,8]),mean(c1_go[c1[,16] >0.001,8])),c(mean(c1_go[c1[,16] <0.001,11]),mean(c1_go[c1[,16] >0.001,11]))))
}



###Look at over-representation island per island###
###Look at over-representation island per island###
###Look at over-representation island per island###
for(a in 1:1)
{

islands_go = cbind(size_all_all[size_all_all[,1] == a,][,c(2,3,4)],0,0,0,0,0,0)
for(c in 1:nrow(islands_go))
#for(c in 1:500)
{
#for island one find out what genes are in there. 

w = 0.001 #window size
island_genes = fst_map_match[fst_map_match[,(col - 2)] >= (islands_go[c,2]-w),];island_genes = island_genes[island_genes[,(col-2)] <= (islands_go[c,3]+w),]; island_genes = unique(island_genes[,1])

island_genes_GO = NULL
islands_go[c,4] = length(island_genes)
for(j in 1:length(island_genes))
{
x = go_annotation_expr[,1] %in% island_genes[j]
if(length(x[x == T]) ==1) island_genes_GO = paste(island_genes_GO,go_annotation_expr[x == T,4], collapse = ",", sep = ",")
}

	for(go in 1) #1 is bp, 2 is mf and three is cc
	{
	go_islands = strsplit(paste(island_genes_GO, collapse = ","),",")[[1]];go_islands = go_islands[regexpr(go_category[go],go_islands) > 0];go_islands = cbind(rle(sort(go_islands))$lengths, rle(sort(go_islands))$values,0,0,0,0,0)
	if(length(go_islands) == 5) go_islands = matrix(0,nrow= 2, ncol = 7)
	colnames(go_islands) = c("number","GO","proportion_observed","expected","pval","adjust_pvalBH","GO_name")

	go_islands[,3] = as.numeric(go_islands[,1]) / sum(as.numeric(go_islands[,1]))

	###FET adjusted
	for(go_is in 1:nrow(go_islands))
	{
	aa = as.numeric(go_islands[go_is,1])
	b = sum(as.numeric(go_islands[,1]))
	cc = as.numeric(go_all[[go]][go_all[[go]][,2] == go_islands[go_is,2],1])
	d = sum(as.numeric(go_all[[go]][,1]))
	go_islands[go_is,5] = signif(fisher.test(cbind(c(aa,b),c(cc,d)), alternative = "g")$p.value,2)
	if((length(cc) != 0) & length((d) != 0)) go_islands[go_is,4] = cc / d
	}
	go_islands[,6] =  p.adjust(as.numeric(go_islands[,5]),"fdr")

	###get the description
	for(des in 1:nrow(go_islands))
		{
		x = description_go[,4] %in% gsub("^..","",go_islands[des,2])
		if(length(x[x == T]) ==1) go_islands[des,7] = description_go[x == T,2]
		}


threshold = 0.05

go_islands_temp = go_islands[as.numeric(go_islands[,6]) < threshold,]

if((go == 1) & (length(go_islands_temp) > 7)) {islands_go[c,7] = nrow(go_islands_temp); islands_go[c,5] =  paste(go_islands_temp[,2],sep = "",collapse = ","); islands_go[c,6] =  paste(go_islands_temp[,7],sep = "",collapse = ",")}
if((go == 1) & (length(go_islands_temp) == 7)) {islands_go[c,7] = 1; islands_go[c,5] =  go_islands_temp[2]; islands_go[c,6] =  go_islands_temp[7]}

}
}






###genome wide distribution scale on the number of genes allowed to sample###
###genome wide distribution scale on the number of genes allowed to sample###
###genome wide distribution scale on the number of genes allowed to sample###
setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

size_all_all = read.table("results/size_all_all",  header = T, stringsAsFactors = F)

if(a == 1) {fst_map_match = read.delim("results/ann_pet_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {fst_map_match = read.delim("results/ann_deb_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {fst_map_match = read.delim("results/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F);c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {fst_map_match = read.delim("results/pet_deb_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {fst_map_match = read.delim("results/pet_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {fst_map_match = read.delim("results/deb_arg_fst_map", header = T, sep = " ", stringsAsFactors = F); c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

#islands_go = cbind(size_all_all[size_all_all[,1] == a,][,c(2,3,4)],0,0,0,0,0,0)
colnames(islands_go)[4:9] = c("number_genes","GO_id","GO_description","GO_over_obs","GO_exp","pvalue")

col = ncol(fst_map_match)
	
unique_genes  = unique(fst_map_match[order(fst_map_match[,1]),c(1,(col-2))]); unique_genes = unique_genes[order(unique_genes[,2]),]
#for(c in 1:nrow(islands_go))
islands_go[,4] = ifelse(islands_go[,4] == 0,1,islands_go[,4])
for(c in 1:nrow(islands_go))
{
go_islands_temp = NULL
number_perm = floor(nrow(unique_genes)/islands_go[c,4])
	for(perm in 0:(number_perm-1))
	{
		island_genes = unique_genes[(perm*islands_go[c,4]+1):((perm+1)*islands_go[c,4]),1]

		island_genes_GO = NULL
		for(j in 1:length(island_genes))
			{
			x = go_annotation_expr[,1] %in% island_genes[j]
			if(length(x[x == T]) ==1) island_genes_GO = paste(island_genes_GO,go_annotation_expr[x == T,4], collapse = ",", sep = ",")
			}

		for(go in 1) #1 is bp, 2 is mf and three is cc
			{
			go_islands = strsplit(paste(island_genes_GO, collapse = ","),",")[[1]];go_islands = go_islands[regexpr(go_category[go],go_islands) > 0];go_islands = cbind(rle(sort(go_islands))$lengths, rle(sort(go_islands))$values,0,0,0,0,0)
			if(length(go_islands) == 5) go_islands = matrix(0,nrow= 2, ncol = 7)
			colnames(go_islands) = c("number","GO","proportion_observed","expected","pval","adjust_pvalBH","GO_name")

			go_islands[,3] = as.numeric(go_islands[,1]) / sum(as.numeric(go_islands[,1]))

		###FET adjusted
			for(go_is in 1:nrow(go_islands))
				{
				aa = as.numeric(go_islands[go_is,1])
				b = sum(as.numeric(go_islands[,1]))
				cc = as.numeric(go_all[[go]][go_all[[go]][,2] == go_islands[go_is,2],1])
				d = sum(as.numeric(go_all[[go]][,1]))
				go_islands[go_is,5] = signif(fisher.test(cbind(c(aa,b),c(cc,d)), alternative = "g")$p.value,2)
				if((length(cc) != 0) & length((d) != 0)) go_islands[go_is,4] = cc / d
				}
			go_islands[,6] =  p.adjust(as.numeric(go_islands[,5]),"fdr")

		###get the description
		threshold = 0.05

		go_islands_temp = c(go_islands_temp,length(go_islands[as.numeric(go_islands[,6]) < threshold,1]))

		}
	
	}
	islands_go[c,8] = mean(go_islands_temp) 
	islands_go[c,9] = length(go_islands_temp[go_islands_temp > islands_go[c,7]]) / perm
	print(islands_go[c,1:3])
}

if(a == 1) islands_go_1 = islands_go
if(a == 2) islands_go_2 = islands_go
if(a == 3) islands_go_3 = islands_go
if(a == 4) islands_go_4 = islands_go
if(a == 5) islands_go_5 = islands_go
if(a == 6) islands_go_6 = islands_go

}




