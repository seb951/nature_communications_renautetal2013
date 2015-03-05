#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
a = as.numeric(args[1])

library(hierfstat)

aa = Sys.time()
setwd("~/Documents/speciation_islands_dec2011") #set up working directory 


###kick out SNP which have more than 10% XX data missing. Allow more for the annuus comparisons given that there is more 454 libraries in annuus which contain much more missing data.strrrr
if(a ==1) {sample = "snp_table_ann_pet"; mis = 0.45; one = "ann_"; two = "pet_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)}
if(a ==2) {sample = "snp_table_ann_deb"; mis = 0.45; one = "ann_"; two = "deb_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)}
if(a ==3) {sample = "snp_table_ann_arg"; mis = 0.37; one = "ann_"; two = "arg_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)} #with 0.35, you get 180K snps....
if(a ==4) {sample = "snp_table_pet_deb"; mis = 0.2; one = "pet_"; two = "deb_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)}
if(a ==5) {sample = "snp_table_pet_arg"; mis = 0.15; one = "pet_"; two = "arg_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)}
if(a ==6) {sample = "snp_table_deb_arg"; mis = 0.07; one = "deb_"; two = "arg_"; name_ind = read.table(paste("results/",gsub("s","s",sample),sep = ""), nrow = 1, header = T)}



system(paste("wc -l results/",sample, " >wordcount1", sep = ""))
wordcount = read.table("wordcount1")

con <- file(paste("results/", gsub("s","s",sample), sep = ""))

open(con) ; file = readLines(con,1); n_ind = length(strsplit(file," ")[[1]]) -3;close(con); write.table(paste(file, "fst fit fis"),paste(gsub("snp_table_","results/",sample),"_comp4", sep = ""), row.names = F, col.names = F, quote = F)
 
con <- file(paste("results/", sample, sep = "")) 
open(con) ; file = readLines(con,1)

for(i in 1:(as.numeric(wordcount[1])-1) )
#for(i in 1:100000)
		{
		file = readLines(con,1)
		file2 = strsplit(file," ")[[1]]
		missing = length(c(1:n_ind)[grepl("XX",file2[4:(n_ind+3)]) == T]) / n_ind

##############################
### trim based on 2pq (Ht) ###
##############################
#counting the alleles.
	counter_ACGTX = rep(0,5)
	names(counter_ACGTX) = c("A","C","G","T","X")
	ht = 0

		x = paste(file2[4:n_ind], collapse = "")
		xx = strsplit(x, split = "")
		
		counter_ACGTX[1] = length(c(1:((n_ind)*2))[grepl("A", xx[[1]])])
		counter_ACGTX[2] = length(c(1:((n_ind)*2))[grepl("C", xx[[1]])])
		counter_ACGTX[3] = length(c(1:((n_ind)*2))[grepl("G", xx[[1]])])
		counter_ACGTX[4] = length(c(1:((n_ind)*2))[grepl("T", xx[[1]])])
		counter_ACGTX[5] = length(c(1:((n_ind)*2))[grepl("X", xx[[1]])])
		p = sort(counter_ACGTX[1:4])[4] 
		q = sort(counter_ACGTX[1:4])[3] 
		al = sum(counter_ACGTX[1:4])
		ht = 2 * (p / al) * (q / al)

###################################
### trim based on Ho (paralogs) ###
###################################	
		ho  = 0
		a1 = substring(file2[4:n_ind],1,1)
		a2 = substring(file2[4:n_ind],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho = (ho + 1) #count the heterozygotes
		}
	ho = ho / n_ind # observed heterozygosity


######## 

if((missing < mis) & (ht > 0.2) & (ho < 0.6)) {

########################
### hierfstat format ###
########################
#fstat_results = cbind(c(1:nrow(comp4)),0,0,0)

#colnames(fstat_results) = c("snp","fst","fit","fis")

	hf = t( rbind( rep("empty",   (length( file2)-3)), file2[4:length(file2)]))
	rownames(hf) = colnames(name_ind)[4:ncol(name_ind)]

	hf[regexpr(one,rownames(hf)) > 0,1] = "1";  hf[regexpr(two,rownames(hf)) >0,1] = "2"

	hf  = gsub("U|R|Y|M|K|S|W|B|D|H|V|N","X",hf);	hf = gsub("X","",hf); hf = gsub("0","",hf);	hf = gsub("A",1, hf);hf = gsub("C",2, hf);hf = gsub("G",3, hf);	hf = gsub("T",4, hf)

	colnames(hf) = c("pop","snp1")

	if(length(hf[hf[,2] == "",2])/length(hf[,2]) < 0.5) {xfst = varcomp(data.frame(as.integer(hf[,1]), as.integer(hf[,2]))); file = paste(file,signif(xfst$F[1,1],4), signif(xfst$F[1,2],4), signif(xfst$F[2,2],4), collapse = " ")} else file2 = c(file2,0,0,0) #Fst #Fit #Fis

}

if((missing < mis) & (ht > 0.2) & (ho < 0.6)) cat(file,file = paste(gsub("snp_table_","results/",sample),"_comp4", sep = ""),append = T, fill = 1)

if(i %% 10000 ==0) print(paste(i, Sys.time()))
		}

close(con)



bb = Sys.time()

