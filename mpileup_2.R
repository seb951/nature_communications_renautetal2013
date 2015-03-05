#!/usr/bin/Rscript --verbose

#
###
###calling SNP###
###



#index the reference file per gene chunks so that you can parralelize the alignments#

#system("vcfutils.pl splitchr -l 100000 reference/HA412_trinity_noAltSplice_400bpmin.fa.fai >reference/split")

###mpileup###

split = read.table("reference/split", stringsAsFactors = F)

cpu = 10 # how many cpus can you afford

for(j in 1:length(split))
#for(j in 1001:2000)
{
system("ps -u renaut | grep 'samtools' | wc -l >prog")
mpileup = paste("samtools mpileup -C50 -I -ugf reference/HA412_trinity_noAltSplice_400bpmin.fa  alignments/*sorted.bam -r ",split[j,1]," | bcftools view -cvg - > mpileup/variants",j,".raw.vcf", sep = "")
Sys.sleep(ifelse(read.table("prog") <= cpu,0, (read.table("prog") - cpu)^4    ))
mpileup_exe = paste("mpileup/cmd/mpileup_",j,sep = "")
write.table(mpileup,mpileup_exe, row.names = F, col.names = F, quote = F)
system(paste("chmod +x",mpileup_exe))
system(paste("nohup ./", mpileup_exe, ">log&",sep = ""))
}

#cat the variants_j.raw.vcf files. 
Sys.sleep(600) #tidy up everything
system("echo -n '' >mpileup/0-11k_variants.raw.vcf") #final file
system("cat mpileup/variants?.raw.vcf  mpileup/variants??.raw.vcf  mpileup/variants???.raw.vcf mpileup/variants????.raw.vcf mpileup/variants10???.raw.vcf >> mpileup/0-11k_variants.raw.vcf") # cat everything
system("grep '#' -v mpileup/0-11k_variants.raw.vcf >mpileup/0-11k_variants.clean.vcf") #tidy up
system("rm mpileup/variants*")
system("rm mpileup/all_variants.raw.vcf")


###

###Step 0.1 set up working directory ###
	setwd("~/Documents/speciation_islands_dec2011")

	system("wc -l mpileup/0-51k_variants.clean.vcf >mpileup/wordcount1")
	wordcount = read.table("mpileup/wordcount1")
	system("awk 'NR == 27' mpileup/variants19022.raw.vcf >mpileup/header")
	
	header = read.delim("mpileup/header", header = F, stringsAsFactors = F)
	header = gsub("alignments/","",header[10:length(header)], fixed = T)
	header = gsub(".fq.sorted.bam","",header, fixed = T)
	write.table(t(as.matrix(c("reference", "position","ref_allele",header))),"mpileup/snp_table",row.names = F, col.names = F, quote = F, sep = " ")

	mis = NULL
	con = file("mpileup/0-51k_variants.clean.vcf")
	open(con)
	for(i in 1:10000) 
	#for(i in 1:as.numeric(wordcount[1]))
	{
		x = readLines(con,1)
		xx = strsplit(x, split = "\t")[[1]]
		xxx = xx[(length(xx)-(length(header) - 1)):length(xx)]
		xxx = 	gsub("0/",paste(xx[4],"/" ,sep = ""),xxx, fixed = T)
		xxx = 	gsub("/0",paste("/",xx[4],sep = ""),xxx, fixed = T)
		xxx = 	gsub("1/",paste(xx[5],"/" ,sep = ""),xxx, fixed = T)
		xxx = 	gsub("/1",paste("/",xx[5],sep = ""),xxx, fixed = T)

			for(q in 1:length(xxx)) # this is to replace low quality calls (maximum Phred-scaled genotype likelihoods below 30). ie. essentially, you need at least 2 reads to call a SNP!
			{
			if(max(as.numeric(strsplit(substring(xxx[q],(gregexpr(":",xxx[q], fixed = T)[[1]][1]+1),(gregexpr(":",xxx[q], fixed = T)[[1]][2]-1)), split = ",")[[1]])) < 30) xxx[q] = "X/X"
			}

		xxx = substring(xxx,1,3)
		if(nchar(xx[5]) != 1) xxx[4:length(xxx)] = "XX"
		xxx = c(xx[c(1,2,4)],gsub("/","",xxx,fixed = T))
		if(length(xxx[xxx == "XX"]) < 80)	cat(t(as.matrix(xxx)),file = "mpileup/snp_table",append = T, fill = F, "\n") else mis = c(mis, i) #only cat loci with less than 75% missing data. 

		if(i %% 1000 == 0) print(paste(i,"of", wordcount[1], Sys.time()))
}
close(con)

system("cat mpileup/snp_table | sed 's/[ \t]*$//' >mpileup/snp_table_2")
system("rm mpileup/snp_table")



####create 6 snp_table for each pairwise###
setwd("/home/seb/Documents/speciation_islands_dec2011") #set up working directory 

snp_table = as.matrix(read.delim("mpileup/snp_table_2", header = T, sep = " ", stringsAsFactors = F))
#snp_table = as.matrix(read.delim("mpileup/snp_table_2", nrow = 10000, header = T, sep = " "))
snp_table[,2] = gsub(" ","",snp_table[,2])

snp_table = snp_table[,-c(16, 27,31, 95,99)]#kick out individuals 81 (PL109.white) and 124 (PI586932b) because they didnt sequence good. In addition, kick out btm15.2, and ARG1820.white which are not what they claim they are. kick out 19-2 which is just a repeat (file copied twice) of another sample. 
colnames(snp_table) = gsub("X14","14",colnames(snp_table)); colnames(snp_table) = gsub("X2O","2O",colnames(snp_table)); 

###
### KICK OUT INDIVIDUALS THAT HAVE MORE THAN 50% MISSING DATA (POOR SEQUENCING RESULTS)
###
#too_much_missing = rep(0,124)
#	for(j in 4:ncol(snp_table))
#		{too_much_missing[(j)-3]  = round(length(c(1:nrow(snp_table))[snp_table[,j] == "XX"  ]) / nrow(snp_table),3)}

################
### create file for each of the 3 comparisons ###
################
all = read.delim("reference/4_species", header = T, stringsAsFactors = F)
ind = colnames(snp_table)

for(i in 4:length(ind))
	{ind[i] = all[regexpr(colnames(snp_table)[i],all[,1]) > 0,2]} # there is a warning because PL109 is present twice. but it doesnt really matter.

colnames(snp_table)[4:ncol(snp_table)] = paste(ind[4:ncol(snp_table)],colnames(snp_table)[4:ncol(snp_table)],sep = "_")

ann_pet =  cbind(snp_table[,1:3],snp_table[,regexpr("ann_|pet_", colnames(snp_table))> 0])
ann_deb = cbind(snp_table[,1:3],snp_table[,regexpr("ann_|deb_", colnames(snp_table))> 0])
ann_arg = cbind(snp_table[,1:3],snp_table[,regexpr("ann_|arg_", colnames(snp_table))> 0]) 
pet_deb = cbind(snp_table[,1:3],snp_table[,regexpr("pet_|deb_", colnames(snp_table))> 0]) 
pet_arg = cbind(snp_table[,1:3],snp_table[,regexpr("pet_|arg_", colnames(snp_table))> 0]) 
deb_arg = cbind(snp_table[,1:3],snp_table[,regexpr("deb_|arg_", colnames(snp_table))> 0]) 

write.table(ann_pet, "results/snp_table_ann_pet", row.names = F, col.names = T, quote = F)
write.table(ann_deb, "results/snp_table_ann_deb", row.names = F, col.names = T, quote = F)
write.table(ann_arg, "results/snp_table_ann_arg", row.names = F, col.names = T, quote = F)
write.table(pet_deb, "results/snp_table_pet_deb", row.names = F, col.names = T, quote = F)
write.table(pet_arg, "results/snp_table_pet_arg", row.names = F, col.names = T, quote = F)
write.table(deb_arg, "results/snp_table_deb_arg", row.names = F, col.names = T, quote = F)





