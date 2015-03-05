#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to count stats on. You should do all that you aligned otherwise, bug!
to = as.numeric(args[2])

setwd("~/Documents/speciation_islands_dec2011/") ### setup the working directory in the right format


##########################################
### count the number of aligned sequences (samtools idxstats) ###
##########################################
idxstats = function(list_ind= "reference/4_species",total_idxstats = NULL, from = 1, to = 1) {

individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals))
	{
	individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
	if(length(grep("fq",individuals[p,1])) != 1) individuals[p,5] = paste(individuals[p,5],".fq",sep = "")
	individuals[p,5] = paste("/home/transfer/Documents/new/bam_files_jan2012/",individuals[p,5],".sorted.bam",sep = "")
	} # proper format to further process with idxstats

	for(i in from:to) {

	#step 1: generate the alignment statistics#
	
	sam_index = paste("samtools index ",individuals[i,5],sep = "") #index bam files
	
	command_idxstats = paste("samtools idxstats ", individuals[i,5]," >idxstats_results",  sep = "")
	system(sam_index)
	system(command_idxstats)
	
	#step 2: read the covage and the alignment stats and then update the result_matrix file
	idxstats_results = read.delim("idxstats_results", stringsAsFactors = F, header = F) #don't read last line.
	
	#total_seq[i,4] = sum(idxstats_results[,3]); total_seq[i,5] = total_seq[i,4] / (total_seq[i,3] * 2) 
	
	if(i == 1) {total_idxstats =  idxstats_results[1:(nrow(idxstats_results)-1),c(1:3)]; colnames(total_idxstats) = c("name","length",paste(strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])],"_number_reads",sep = "") )}
	if(i > 1) {total_idxstats = cbind(total_idxstats,idxstats_results[1:(nrow(idxstats_results)-1),3]); colnames(total_idxstats)[ncol(total_idxstats)] = paste(strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]
,"_number_reads",sep = "")}
	system("rm /home/transfer/Documents/new/bam_files_jan2012/*bai")
	}
	system("rm idxstats_results")
	colnames(total_idxstats) = gsub(".fq_number_reads","",colnames(total_idxstats))
	write.table(total_idxstats,"results/coverage_per_gene_individuals.txt",row.names = F,col.names = T)

}

############
### running ###
############
idxstats(from = 1, to = 108) ### count the number of aligned sequences



############
### extra stats you may want ###
############



if(file.exists("results/total_idxstats.txt")) {
total_idxstats = read.table("results/total_idxstats.txt",header = T, stringsAsFactors = F)

mean_med = NULL
for(i in 3:ncol(total_idxstats))
{
mean_med = rbind(mean_med,paste(signif(mean(total_idxstats[,i]),4)," (",median(total_idxstats[,i]),")",sep = ""))
}

write.table(mean_med,"results/mean_med", row.names = F, col.names = F, quote = F)
}




