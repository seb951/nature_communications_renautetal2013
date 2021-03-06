#!/usr/bin/Rscript --verbose

args = commandArgs(TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

for(c in from:to)
{
###########################
### genome wide MK test ###
###########################
### loading required files ###
setwd("~/Documents/speciation_islands_dec2011")

ann_pet_fst = read.table("results/ann_pet_comp4", header = T, stringsAsFactors = F)
ann_deb_fst = read.table("results/ann_deb_comp4",header = T, stringsAsFactors = F)
ann_arg_fst = read.table("results/ann_arg_comp4", header = T, stringsAsFactors = F)
pet_deb_fst = read.table("results/pet_deb_comp4", header = T, stringsAsFactors = F)
pet_arg_fst = read.table("results/pet_arg_comp4",   header = T, stringsAsFactors = F)
deb_arg_fst = read.table("results/deb_arg_comp4", header = T, stringsAsFactors = F)
all_comparisons_4 = list(ann_pet_fst,ann_deb_fst,ann_arg_fst,pet_deb_fst,pet_arg_fst,deb_arg_fst)

mk = NULL

ann_deb_fst = all_comparisons_4[[c]]
col = ncol(ann_deb_fst )
ann_deb_fst[,(col-2)] = jitter(ann_deb_fst[,(col-2)],0.1) #add/remove between -/+ 0-2e-05. to all fst values 
top3 = quantile(ann_deb_fst[,(col-2)],0.97,na.rm = T)
#top5 = quantile(ann_deb_fst[,(col-2)],0.95,na.rm = T)
###
#create a table of the consensus call for annuus and debilis.
###
counter_ACGTX = matrix(0,nrow = nrow(ann_deb_fst), ncol = 5)
counter_ACGTX1 = matrix(0,nrow = nrow(ann_deb_fst), ncol = 5)
colnames(counter_ACGTX) = c("A","C","G","T","X")
colnames(counter_ACGTX1) = c("A","C","G","T","X")
consensus_final =  NULL

ceil = ceiling(nrow(ann_deb_fst)/1000)

for(a in 1:ceil)
#for(a in 1:5)
{
ann_deb_fst_sub = ann_deb_fst[((a*1000 )-999):(a*1000),]
if(a == ceil) ann_deb_fst_sub = ann_deb_fst_sub[!is.na(ann_deb_fst_sub[,1]),]
consensus = cbind(ann_deb_fst_sub[,c(1,2)],0,0,0,0,0,0,0,0,0)
colnames(consensus) = c("contig","site","tex","deb","fixed", "syn","nonsyn","shared","private","noncoding","top3")
col = ncol(ann_deb_fst_sub)
if(c == 1) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("ann_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("pet_", colnames(ann_deb_fst_sub))> 0]}
if(c == 2) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("ann_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("deb_", colnames(ann_deb_fst_sub))> 0]}
if(c == 3) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("ann_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("arg_", colnames(ann_deb_fst_sub))> 0]}
if(c == 4) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("pet_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("deb_", colnames(ann_deb_fst_sub))> 0]}
if(c == 5) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("pet_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("arg", colnames(ann_deb_fst_sub))> 0]}
if(c == 6) {one = c(1:ncol(ann_deb_fst_sub))[regexpr("deb_", colnames(ann_deb_fst_sub))> 0]; two = c(1:ncol(ann_deb_fst_sub))[regexpr("arg_", colnames(ann_deb_fst_sub))> 0]}

	for(i in 1:nrow(ann_deb_fst_sub))
		{
		if(!is.na(ann_deb_fst_sub[i,(col - 2)]) && (ann_deb_fst_sub[i,(col - 2)] >= 0.9)) consensus[i,5] = 1 # filter for fixed snps... FST###
		if(!is.na(ann_deb_fst_sub[i,(col - 2)]) && (ann_deb_fst_sub[i,(col - 2)] >= top3)) consensus[i,11] = 1 # filter for top3 snps... FST###
		
		#majority base call for ONE. 	
		x = paste(ann_deb_fst_sub[i,one], collapse = "")
		xx = strsplit(x, split = "")
		counter_ACGTX[i,1] = length(c(1:nchar(x))[grepl("A", xx[[1]])])
		counter_ACGTX[i,2] = length(c(1:nchar(x))[grepl("C", xx[[1]])])
		counter_ACGTX[i,3] = length(c(1:nchar(x))[grepl("G", xx[[1]])])
		counter_ACGTX[i,4] = length(c(1:nchar(x))[grepl("T", xx[[1]])])
		counter_ACGTX[i,5] = length(c(1:nchar(x))[grepl("X", xx[[1]])])
		if(consensus[i,5] == 1) consensus[i,3] = paste(names(sort(counter_ACGTX[i,1:4])[4]),names(sort(counter_ACGTX[i,1:4])[4]),sep = "") 
		if(sort(counter_ACGTX[i,1:4])[3] == 0) consensus[i,3] = paste(names(sort(counter_ACGTX[i,1:4])[4]),names(sort(counter_ACGTX[i,1:4])[4]),sep = "")
		if(sort(counter_ACGTX[i,1:4])[3] != 0) consensus[i,3] = paste(names(sort(counter_ACGTX[i,1:4])[3:4]),collapse = "")
		
		#majority base call for TWO. 
		x1 = paste(ann_deb_fst[i,two], collapse = "")
		xx1 = strsplit(x1, split = "")
		counter_ACGTX1[i,1] = length(c(1:nchar(x1))[grepl("A", xx1[[1]])])
		counter_ACGTX1[i,2] = length(c(1:nchar(x1))[grepl("C", xx1[[1]])])
		counter_ACGTX1[i,3] = length(c(1:nchar(x1))[grepl("G", xx1[[1]])])
		counter_ACGTX1[i,4] = length(c(1:nchar(x1))[grepl("T", xx1[[1]])])
		counter_ACGTX1[i,5] = length(c(1:nchar(x1))[grepl("X", xx1[[1]])])
		if(consensus[i,5] == 1) consensus[i,4] = paste(names(sort(counter_ACGTX1[i,1:4])[4]),names(sort(counter_ACGTX1[i,1:4])[4]),sep = "")
		if(sort(counter_ACGTX1[i,1:4])[3] == 0) consensus[i,4] = paste(names(sort(counter_ACGTX1[i,1:4])[4]),names(sort(counter_ACGTX1[i,1:4])[4]),sep = "")
		if(sort(counter_ACGTX1[i,1:4])[3] != 0) consensus[i,4] = paste(names(sort(counter_ACGTX1[i,1:4])[3:4]),collapse = "")

		#shared or private snps?
		if(paste(sort(strsplit(consensus[i,3],"")[[1]]),collapse = "") == paste(sort(strsplit(consensus[i,4],"")[[1]]), collapse = "")) consensus[i,8] = 1 #shared
		if((paste(sort(strsplit(consensus[i,3],"")[[1]]),collapse = "") != paste(sort(strsplit(consensus[i,4],"")[[1]]), collapse = "")) & ((substring(consensus[i,3],1,1) == substring(consensus[i,3],2,2)) | (substring(consensus[i,4],1,1) == substring(consensus[i,4],2,2)))) consensus[i,9] = 1 #private
		if((consensus[i,4] == consensus[i,4]) &   (substring(consensus[i,3],1,1) == substring(consensus[i,3],2,2)) &    ( substring(consensus[i,4],1,1) == substring(consensus[i,4],2,2)) ) consensus[i,5:10] = 0  #not polymorphic		
		}
		if((a %% 1) == 0) print(paste(a, Sys.time()))
		consensus_final = rbind(consensus_final,consensus)

		}
	
#####################
### FIND BEST ORF ###
#####################
		
### RUN GETORF PACKAGE FROM THE EMBOSS PIPELINE ###
reference_transcriptome = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F,sep = "\t"))
reference = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F,sep = "\t"))
reference_transcriptome_unique_ID = reference_transcriptome
reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1] = paste(">",c(1:length(reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1])), sep = "") 
write.table(reference_transcriptome_unique_ID,"reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa", row.names = F, col.names = F, quote = F)

### run getorf ###
system("/usr/local/bin/getorf -sequence reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa -minsize 300 -find 2 -outseq reference/HA412_trinity_noAltSplice_orf.fasta")

### PARSE THE OUTPUT OF GETORF###
out = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_orf.fasta", header = F, sep = " "))
out = out[,1:5]
out[,2] = gsub("^.","", out[,2])
out[,4] = gsub(".$","", out[,4])
out = gsub("SENSE)","", out, fixed = T)
out = cbind(gsub("(REVERSE","", out, fixed = T), 0)
out = rbind(out,">")

x = c(1:nrow(out))[(regexpr(">",out[,1],fixed = T) > 0)]

for(i in 1: (length(x)-1)) {out[x[i],5]  = paste(out[(x[i]+1): (x[(i+1)]-1),1], collapse = "")}

out_seq = out[(regexpr(">",out[,1],fixed = T) > 0),]
out_seq[,3] =  abs(as.numeric(out_seq[,2]) - as.numeric(out_seq[,4]))
out_seq[,6] = apply(out_seq,2,substring,(regexpr(">", out_seq[,1],fixed = T)+1),(regexpr("_", out_seq[,1],fixed = T)-1))[,1]
unique_gene = unique(out_seq[,6])
unique_orf = NULL

for(i in 1:length(unique_gene))
	{
	x = out_seq[,6] %in% unique_gene[i]
	temp = out_seq[x == T, ]
	if(length(temp) == 6) longest_temp = temp[as.numeric(temp[3]) == max(as.numeric(temp[3]))] else longest_temp = temp[as.numeric(temp[,3]) == max(as.numeric(temp[,3])),]
	if(length(longest_temp) == 6) unique_orf = rbind(unique_orf, longest_temp) else unique_orf = rbind(unique_orf,longest_temp[1,])
	}

	unique_orf = unique_orf[1:(nrow(unique_orf)-1),] #These are the longest ORF.

### get the proper names back
	names = reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1] 
	for(i in 1:nrow(unique_orf)) {unique_orf[i,1] = names[as.numeric(unique_orf[i,6])]}
	unique_orf[,6] = gsub(">","",unique_orf[,1], fixed = T)

###	
#create a objet reference_matrix which contains the annotated consensus sequences. 
###
reference_vector = c(1:nrow(reference))[(regexpr(">", reference[,1],fixed = T) > 0)]
reference_vector = cbind(reference[(regexpr(">", reference[,1],fixed = T) > 0),1], reference_vector,0)
reference_vector = rbind(reference_vector,c("null",nrow(reference),0))

for(i in 1:(nrow(reference_vector)-1))
{reference_vector[i,3] = paste(reference[(as.numeric(reference_vector[i,2])+1): (as.numeric(reference_vector[(i+1),2])-1),1],collapse = "")}

reference_matrix = cbind(gsub(">","",reference_vector[1:nrow(reference_vector),1]), reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3])
colnames(reference_matrix) = c("name","consensus","ONE_haplo1","ONE_haplo2","TWO_haplo1","TWO_haplo2")

consensus_final_per_gene_template = cbind(reference_matrix[,1],0,0,0,nchar(reference_matrix[,2]),0)
colnames(consensus_final_per_gene_template) = c("gene","s","ns","nc","length_total","length_ORF")

write.table(consensus_final_per_gene_template,"results/consensus_final_per_gene_template.txt", row.names =F, col.names = T)


###
### annotate reference matrix with the consensus object.
###

for(i in 1:nrow(reference_matrix))

#for(i in 1:500)
	{
	x = consensus_final[,1] %in% reference_matrix[i,1]
	y = consensus_final[x == T, ]

	if(nrow(y) == 1) 
	{
	substring(reference_matrix[i,3], as.numeric(y[2]), as.numeric(y[2])) = substring(y[3],1,1);
	substring(reference_matrix[i,4], as.numeric(y[2]), as.numeric(y[2])) = substring(y[3],2,2);
	substring(reference_matrix[i,5], as.numeric(y[2]), as.numeric(y[2])) = substring(y[4],1,1);
	substring(reference_matrix[i,6], as.numeric(y[2]), as.numeric(y[2])) = substring(y[4],2,2)
	}

	if(nrow(y) > 1) for(j in 1:nrow(y))
	{
	substring(reference_matrix[i,3], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,3],1,1)
	substring(reference_matrix[i,4], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,3],2,2)
	substring(reference_matrix[i,5], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,4],1,1)
	substring(reference_matrix[i,6], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,4],2,2)
	}

	if((i %% 1000) == 0) print(paste(i, Sys.time()))
}

###
#Interrogate each SNP, find what codon it belongs to and see whether it is syn or nonsyn. add this information in the consensus table. 
###
aa = as.matrix(read.delim("reference/aa_code.txt", header = T))
aa_temp = c(0,0,0,0)
symbols =      c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B")
replacements = c("t","g","c","a","k","y","w","s","r","m","b","d","h","v")
codon = NULL
for(i in 1:nrow(consensus_final))
#for(i in 305:344)
{
	x_seq = reference_matrix[,1] %in% consensus_final[i,1]
	x_orf = unique_orf[,6] %in% consensus_final[i,1]
	y_seq = reference_matrix[x_seq == T, ]
	y_orf = unique_orf[x_orf == T, ]

		if(as.numeric(y_orf[2]) < as.numeric(y_orf[4]) & (length(y_orf) > 0))  {# the easy case
			codon_positions = seq(as.numeric(y_orf[2]), as.numeric(y_orf[4]), 3); #all start positions of the codons. 
			pos_temp = codon_positions[codon_positions <= as.numeric(consensus_final [i,2])]; 
			pos = pos_temp[length(pos_temp)];
			if(length(pos) > 0) codon = substring(y_seq[3:6], pos, (pos +2))
			if(length(pos) == 0) {codon = c(0,0,0,0); consensus_final [i,10] = 87} #you are before the start position
			if((length(pos) != 0) && (pos < (as.numeric(consensus_final[i,2])-2))) {codon = c(0,0,0,0); consensus_final [i,10] = 88}} #you are after the end position
		
		if(as.numeric(y_orf[2]) > as.numeric(y_orf[4]) & (length(y_orf) > 0))  {# the hard reverse complement case
			codon_positions = seq(as.numeric(y_orf[2]), as.numeric(y_orf[4]), -3); #all start positions of the codons. 
			pos_temp = codon_positions[codon_positions >= as.numeric(consensus_final [i,2])]; 
			pos = pos_temp[length(pos_temp)];
			if(length(pos) == 0) {codon = c(0,0,0,0); consensus_final [i,10] = 98} #after the end
			if(length(pos) > 0) {codon = substring(y_seq[3:6], (pos - 2), (pos));			
			codon[1] = paste(rev(strsplit(codon,"")[[1]]), collapse = "");
			codon[2] = paste(rev(strsplit(codon,"")[[2]]), collapse = "");
			codon[3] = paste(rev(strsplit(codon,"")[[3]]), collapse = "");
			codon[4] = paste(rev(strsplit(codon,"")[[4]]), collapse = "");
			for(s in 1:length(symbols)) codon = gsub(symbols[s], replacements[s], codon)} #complement sequence.
			codon = gsub("(\\w)", "\\U\\1", codon, perl=TRUE) #gsub for small caps to capital letters. 
			if((length(pos) != 0) && (pos > (as.numeric(consensus_final[i,2])+2))) {codon = c(0,0,0,0); consensus_final[i,10] = 99} #before the start
			}
		
		if(length(y_orf) == 0) 	 {codon = c(0,0,0,0); consensus_final [i,10] = 1}
		# find what are the amino acid encoded. Don't use codons with ambuiguities in the reference sequence...
		if((length(unique(codon)) > 1) & (unique(codon)[1] != 0) & (-grepl("M|R|W|S|Y|K|V|H|D|B|N", codon[1]) == 0) & (-grepl("M|R|W|S|Y|K|V|H|D|B|N", codon[2]) == 0) & (-grepl("M|R|W|S|Y|K|V|H|D|B|N", codon[3]) == 0) & (-grepl("M|R|W|S|Y|K|V|H|D|B|N", codon[4]) == 0)) { 
			aa_temp[1] = aa[aa[,2] == codon[1],1];
			aa_temp[2] = aa[aa[,2] == codon[2],1];
			aa_temp[3] = aa[aa[,2] == codon[3],1];
			aa_temp[4] = aa[aa[,2] == codon[4],1]} else aa_temp = c(0,0,0,0)
		
		if((length(unique(aa_temp)) == 1) & (aa_temp[1] != 0) ) consensus_final [i,6] = 1 #syn snp
		if((length(unique(aa_temp)) > 1) & (aa_temp[1] != 0) ) consensus_final [i,7] = 1 #nonsyn snp?
		if(aa_temp[1] == "END" | aa_temp[2] == "END" | aa_temp[3] == "END" | aa_temp[4] == "END")consensus_final [i,10] = 2 #ALTERNATIVE STOPS ARE CODED AS 2 IN THE NONCODING COLUMN
		
		if((i %% 1000) == 0) print(paste(i, "of", nrow(consensus_final),Sys.time()))
}

if(c == 1) write.table(consensus_final,"results/consensus_final_ann_pet.txt", row.names =F)
if(c == 2) write.table(consensus_final,"results/consensus_final_ann_deb.txt", row.names =F)
if(c == 3) write.table(consensus_final,"results/consensus_final_ann_arg.txt", row.names =F)
if(c == 4) write.table(consensus_final,"results/consensus_final_pet_deb.txt", row.names =F)
if(c == 5) write.table(consensus_final,"results/consensus_final_pet_arg.txt", row.names =F)
if(c == 6) write.table(consensus_final,"results/consensus_final_deb_arg.txt", row.names =F)


###
### LIST OF THE MUTATIONS BASED ON CODING, SYN, NONSYN
###
if(c == 1) consensus_final = read.delim("results/consensus_final_ann_pet.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 2) consensus_final = read.delim("results/consensus_final_ann_deb.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 3) consensus_final = read.delim("results/consensus_final_ann_arg.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 4) consensus_final = read.delim("results/consensus_final_pet_deb.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 5) consensus_final = read.delim("results/consensus_final_pet_arg.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 6) consensus_final = read.delim("results/consensus_final_deb_arg.txt", header = T,stringsAsFactors = F, sep = " ")

consensus_final_per_gene = read.delim("results/consensus_final_per_gene_template.txt", header = T,stringsAsFactors = F, sep = " ")
consensus_final_per_gene = cbind(consensus_final_per_gene,0)
colnames(consensus_final_per_gene)[7] = "SNP_per_ncsite"
#tex_deb_res = read.delim("tex_deb_res", header = T,stringsAsFactors = F, sep = " ")
#tex_deb_pet_fst = read.delim("../pileup_full/tex_deb_pet_fst", header = T,stringsAsFactors = F, sep = "\t")
for(i in 1:nrow(consensus_final_per_gene))
#for(i in 1:100)
{
	x_snp = consensus_final[consensus_final[,1] == consensus_final_per_gene[i,1],]
	consensus_final_per_gene[i,4] = sum(x_snp[,10]) #non coding
	consensus_final_per_gene[i,2] = sum(x_snp[,6]) #synonymous
	consensus_final_per_gene[i,3] = sum(x_snp[,7]) #non-syn
	consensus_final_per_gene[i,7] =round(consensus_final_per_gene[i,4] / (consensus_final_per_gene[i,5] - consensus_final_per_gene[i,6]),6) #non coding polymorphism per non coding site. 
	if(i %% 1000 == 0) print(paste(i,"of",nrow(consensus_final_per_gene), Sys.time()))
}

if(c == 1) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_ann_pet.txt", row.names =F)
if(c == 2) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_ann_deb.txt", row.names =F)
if(c == 3) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_ann_arg.txt", row.names =F)
if(c == 4) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_pet_deb.txt", row.names =F)
if(c == 5) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_pet_arg.txt", row.names =F)
if(c == 6) write.table(consensus_final_per_gene,"results/consensus_final_per_gene_deb_arg.txt", row.names =F)


####
#####MK test!!!!!!!
####
#
####BOOTSTRAP function###VERSION 2
library(boot)
boot.alpha =  function(alpha_gene,num) {
alpha_sample = alpha_gene[num,]
b.alpha = 1-(mean(alpha_sample[,2]) / mean(alpha_sample[,3])  * mean(alpha_sample[,6]) ) #calculate alpha
if(b.alpha != -Inf) return(b.alpha) #remove -Inf if there woudl be some.
}
#
mk = NULL
for(c in 5:5)
{
if(c == 1) consensus_final = read.delim("results/consensus_final_ann_pet.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 2) consensus_final = read.delim("results/consensus_final_ann_deb.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 3) consensus_final = read.delim("results/consensus_final_ann_arg.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 4) consensus_final = read.delim("results/consensus_final_pet_deb.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 5) consensus_final = read.delim("results/consensus_final_pet_arg.txt", header = T,stringsAsFactors = F, sep = " ")
if(c == 6) consensus_final = read.delim("results/consensus_final_deb_arg.txt", header = T,stringsAsFactors = F, sep = " ")

main = list(c("ann-pet(sympatric-ancient)"),c("ann-deb(parapatric-ancient)"),c("ann-arg(allopatric-recent)"),c("pet-deb (allopatric-recent)"),c("pet-arg(allopatric-ancient)"),c("deb-arg (allopatric-ancient)") )

aa = sum(as.numeric(consensus_final[consensus_final[,11] != 1,6])) #syn non-fixed top3
bb = sum(as.numeric(consensus_final[consensus_final[,11] != 1,7])) #non-syn non-fixed
cc = sum(as.numeric(consensus_final[consensus_final[,11] == 1,6])) #syn fixed
dd = sum(as.numeric(consensus_final[consensus_final[,11] == 1,7])) #=non-syn fixed
#
source(url("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r"))
x_mk = g.test(matrix(c(aa,bb,cc,dd), nrow = 2)) #real mk test
#
###per gene###
alpha_gene = data.frame(unique(consensus_final[,1]),cbind(rep(0,length(unique(consensus_final[,1]))),0,0,0,0), stringsAsFactors = F)
colnames(alpha_gene) = c("gene","Ds","Dn","Ps","Pn", "Pn / (Ps+1))")
#
for(i in 1: nrow(alpha_gene))
#for(i in 1:100)
{
x = consensus_final[,1] %in% alpha_gene[i,1]; temp = consensus_final[x == T, ]
alpha_gene[i,2] = sum(temp[temp[,5] == 1,6])
alpha_gene[i,3] = sum(temp[temp[,5] == 1,7])
alpha_gene[i,4] = sum(temp[temp[,5] != 1,6])
alpha_gene[i,5] = sum(temp[temp[,5] != 1,7])
alpha_gene[i,6] = as.numeric(alpha_gene[i,5]) / (1+as.numeric(alpha_gene[i,4]))
if(i %% 2000 == 0) print(paste(i,"of", nrow(alpha_gene),Sys.time()))
}

boot_Rfunction = boot(alpha_gene,boot.alpha, R=1000)
boot_ci = boot.ci(boot_Rfunction, type = "basic")
#
mk = cbind(mk,c(aa,bb,cc,dd,x_mk$statistic, x_mk$p.value, (1 - (bb * cc) / (aa*dd)),  boot_ci$t0,  boot_ci$basic[4],boot_ci$basic[5]))
colnames(mk)[ncol(mk)] =  main[[c]]
rownames(mk)  = c("syn","nonsyn","syn_top","nonsyn_top","gstat","pval", "alpha","alpha_boot","CImin", "CIplus")
print(paste(c,Sys.time()))
}
mk2 =  as.data.frame(cbind(t(mk),c(1:6)))

colnames(mk2)[11] = "comp"
#write.table(mk2, "results/mk3_200boot",row.names = F, col.names = T)

#
####PLOT
#library(ggplot2)
#library(gplots)
#mk2 = as.matrix(read.delim("results/mk3_200boot", sep = " "))
#mk3 = as.data.frame(mk2[c(1,2,6,5),])
#
#plotCI(c(1:4),mk3[,8],ui = mk3$CIplus,li=mk3$CImin,  ylim = c(0,0.6), xlab = "comparisons", ylab = "alpha",xaxt= "n", lwd = 5)
#axis(1,at =  c(1:4), labels = c(1,2,6,5))
#
##limits = aes(ymax = CIplus, ymin= CImin)
##p = ggplot(mk3,aes(x= 1:4,y = alpha_boot) )
##p + geom_bar(stat = "identity", aes(fill = factor(1:4))) + geom_errorbar(limits,width=0.25)
#
#dev.print(device=svg, "~/Desktop/boot200.svg", onefile=FALSE)
#
#dev.off()

}


