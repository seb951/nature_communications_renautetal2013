#!/usr/bin/Rscript --verbose

### loading required files ###
setwd("~/Documents/speciation_islands_dec2011")

ann_pet_fst = read.table("results/ann_pet_comp4", header = T, stringsAsFactors = F)
ann_deb_fst = read.table("results/ann_deb_comp4",header = T, stringsAsFactors = F)
ann_arg_fst = read.table("results/ann_arg_comp4", header = T, stringsAsFactors = F)
pet_deb_fst = read.table("results/pet_deb_comp4", header = T, stringsAsFactors = F)
pet_arg_fst = read.table("results/pet_arg_comp4",   header = T, stringsAsFactors = F)
deb_arg_fst = read.table("results/deb_arg_comp4", header = T, stringsAsFactors = F)
all_comparisons_4 = list(ann_pet_fst,ann_deb_fst,ann_arg_fst,pet_deb_fst,pet_arg_fst,deb_arg_fst)

reference =  as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa",header = F, sep = "\t")) #reference sequences

for(c in 4:6)
{
setwd("~/Documents/speciation_islands_dec2011")

ann_deb_fst = all_comparisons_4[[c]]
ann_deb_fst = ann_deb_fst[1:10000,]

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

top3 = quantile(ann_deb_fst[,(col-2)],0.97,na.rm = T)#watch this. THERE IS NO JITTER FUNCTION HERE###

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
		if(consensus[i,5] == 1) consensus[i,3] = names(sort(counter_ACGTX[i,1:4])[4])
		if(sort(counter_ACGTX[i,1:4])[3] == 0) consensus[i,3] =names(sort(counter_ACGTX[i,1:4])[4])
		if(sort(counter_ACGTX[i,1:4])[3] != 0) consensus[i,3] =names(sort(counter_ACGTX[i,1:4])[4])
	
		#majority base call for TWO. 
		x1 = paste(ann_deb_fst[i,two], collapse = "")
		xx1 = strsplit(x1, split = "")
		counter_ACGTX1[i,1] = length(c(1:nchar(x1))[grepl("A", xx1[[1]])])
		counter_ACGTX1[i,2] = length(c(1:nchar(x1))[grepl("C", xx1[[1]])])
		counter_ACGTX1[i,3] = length(c(1:nchar(x1))[grepl("G", xx1[[1]])])
		counter_ACGTX1[i,4] = length(c(1:nchar(x1))[grepl("T", xx1[[1]])])
		counter_ACGTX1[i,5] = length(c(1:nchar(x1))[grepl("X", xx1[[1]])])
		if(consensus[i,5] == 1) consensus[i,4] = names(sort(counter_ACGTX1[i,1:4])[4])
		if(sort(counter_ACGTX1[i,1:4])[3] == 0) consensus[i,4] = names(sort(counter_ACGTX1[i,1:4])[4])
		if(sort(counter_ACGTX1[i,1:4])[3] != 0) consensus[i,4] = names(sort(counter_ACGTX1[i,1:4])[4])

		#shared or private snps?
		if(paste(sort(strsplit(consensus[i,3],"")[[1]]),collapse = "") == paste(sort(strsplit(consensus[i,4],"")[[1]]), collapse = "")) consensus[i,8] = 1 #shared
		if((paste(sort(strsplit(consensus[i,3],"")[[1]]),collapse = "") != paste(sort(strsplit(consensus[i,4],"")[[1]]), collapse = "")) & ((substring(consensus[i,3],1,1) == substring(consensus[i,3],2,2)) | (substring(consensus[i,4],1,1) == substring(consensus[i,4],2,2)))) consensus[i,9] = 1 #private
		if((consensus[i,4] == consensus[i,4]) &   (substring(consensus[i,3],1,1) == substring(consensus[i,3],2,2)) &    ( substring(consensus[i,4],1,1) == substring(consensus[i,4],2,2)) ) consensus[i,5:10] = 0  #not polymorphic		
		}
		if((a %% 100) == 0) print(paste(a, Sys.time()))
		consensus_final = rbind(consensus_final,consensus)
}
		
	
	
#####################
### FIND BEST ORF ###
#####################
		
### RUN GETORF PACKAGE FROM THE EMBOSS PIPELINE ###
reference_transcriptome = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F,sep = "\t"))
reference_transcriptome_unique_ID = reference_transcriptome
reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1] = paste(">",c(1:length(reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1])), sep = "") 
write.table(reference_transcriptome_unique_ID,"reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa", row.names = F, col.names = F, quote = F)

### run getorf ###
#system("/usr/local/bin/getorf -sequence reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa -minsize 300 -find 2 -outseq reference/HA412_trinity_noAltSplice_orf.fasta")

### PARSE THE OUTPUT OF GETORF###
out = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_orf.fasta", header = F, sep = " "))
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
	{
	x = consensus_final[,1] %in% reference_matrix[i,1]
	y = consensus_final[x == T, ]

	if(nrow(y) == 1) 
	{
	substring(reference_matrix[i,3], as.numeric(y[2]), as.numeric(y[2])) = substring(y[3],1,1);
	substring(reference_matrix[i,5], as.numeric(y[2]), as.numeric(y[2])) = substring(y[4],1,1);
	}

	if(nrow(y) > 1) for(j in 1:nrow(y))
	{
	substring(reference_matrix[i,3], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,3],1,1)
	substring(reference_matrix[i,5], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,4],1,1)
		}

	if((i %% 10000) == 0) print(paste(i, Sys.time()))
}

###
#####create a new matrix which contains the ORF consensus...
symbols =      c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B")
replacements = c("t","g","c","a","k","y","w","s","r","m","b","d","h","v")

orf_positions = unique_orf[,c(6,2,4)] #These are the position of the ORF.
orf_consensus = cbind(reference_matrix[,1],0,0)#This object will contain the SNP annotated ORFs.
colnames(orf_consensus) = c("reference","ann_cons1","deb_cons")

for(i in 1:nrow(reference_matrix))
{

x = unique_orf[,6] %in% reference_matrix[i,1]
x_positions = unique_orf[x == T, c(1,2,4)]

if((length(x_positions) == 3) & (as.numeric(x_positions[2]) < as.numeric(x_positions[3]))) orf_consensus[i,2:3] =  substring(reference_matrix[i,c(3,5)],as.numeric(x_positions[2]),as.numeric(x_positions[3]))
if((length(x_positions) == 3) & (as.numeric(x_positions[2]) > as.numeric(x_positions[3]))) # the orf for the reverse complement case.
	{
		orf_consensus[i,2:3] = substring(reference_matrix[i,c(3,5)],as.numeric(x_positions[3]),as.numeric(x_positions[2])) #substring the ORF
		orf_consensus[i,3] = paste(rev(strsplit(orf_consensus[i,3],"")[[1]]),collapse = "")  # ORF in the reverse order
		orf_consensus[i,2] = paste(rev(strsplit(orf_consensus[i,2],"")[[1]]),collapse = "")  # ORF in the reverse order

	for(s in 1:length(symbols)) orf_consensus[i,2:3] = gsub(symbols[s], replacements[s], orf_consensus[i,2:3]) #Complement sequence. 
	}
}
orf_consensus[,2:3] = toupper(orf_consensus[,2:3]) #gsub for small caps to capital letters. 


####This is to check for stop codons and cut the sequence until the STOP codon appears### ####STOP CODONS: TAA TRA TAG TGA#### 
####This is to check for stop codons and cut the sequence until the STOP codon appears### ####STOP CODONS: TAA TRA TAG TGA####
####This is to check for stop codons and cut the sequence until the STOP codon appears### ####STOP CODONS: TAA TRA TAG TGA####

premature_stop = NULL

for(i in 1: nrow(orf_consensus))
{
	for(t in 1: (nchar(orf_consensus[i,3])/3))
		{
			for(k in 2:3)
			{
		if((substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAA") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TRA")  | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAG") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TGA")) {premature_stop = c(premature_stop, orf_consensus[i,1]); orf_consensus[i,2:3] = substring(orf_consensus[i,2:3],1,((t-1)*3)) }
  			}
  		}	
  if((i %% 10000) == 0) print(c(i,"of", nrow(orf_consensus)))
}

#####PAML FORMAT######
#####PAML FORMAT######
#####PAML FORMAT######
#####PAML FORMAT######
#####PAML FORMAT######

temp = c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3]))

orf_length = nchar(orf_consensus[,3]) #number of character in each sequence
orf_paml = NULL # this will be the object which contains all the sequences
name_paml = NULL # this will contain all the names 

for(i in 1: nrow(orf_consensus))
	{
	if(orf_length[i] > 300) temp = c(paste(2,orf_length[i]),c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3])))
	if(orf_length[i] > 300) orf_paml = c(orf_paml, temp)
	if(orf_length[i] > 300) name_paml = c(name_paml, orf_consensus[i,1])
	}


setwd("/home/seb/Documents/PAML/ann_deb_pet_arg_RESULTS/")

if(c == 1) {write.table(c(orf_paml,"//end", name_paml),"ann_pet/ann_pet", row.names = F, col.names = F, quote = F); sp = "ann_pet"; setwd(sp)}
if(c == 2) {write.table(c(orf_paml,"//end", name_paml),"ann_deb/ann_deb", row.names = F, col.names = F, quote = F); sp = "ann_deb"; setwd(sp)}
if(c == 3) {write.table(c(orf_paml,"//end", name_paml),"ann_arg/ann_arg", row.names = F, col.names = F, quote = F); sp = "ann_arg"; setwd(sp)}
if(c == 4) {write.table(c(orf_paml,"//end", name_paml),"pet_deb/pet_deb", row.names = F, col.names = F, quote = F); sp = "pet_deb"; setwd(sp)}
if(c == 5) {write.table(c(orf_paml,"//end", name_paml),"pet_arg/pet_arg", row.names = F, col.names = F, quote = F); sp = "pet_arg"; setwd(sp)}
if(c == 6) {write.table(c(orf_paml,"//end", name_paml),"deb_arg/deb_arg", row.names = F, col.names = F, quote = F); sp = "deb_arg"; setwd(sp)}

system(paste("awk '{sub(/XXX1/,",length(orf_paml[orf_paml == "ann_cons1"]),")}; 1'"," ../","codeml.ctl  >x_nd", sep = "")) # 	ndata = number of sequences
system(paste("awk '{sub(/XXX2/,","\"","/home/seb/Documents/PAML/ann_deb_pet_arg_RESULTS/",sp,"/",sp,"\"",")}; 1' x_nd >x_ndseq", sep = "")) # 	where is sequence data
system(paste("awk '{sub(/XXX3/,","\"","/home/seb/Documents/PAML/ann_deb_pet_arg_RESULTS/ann_deb.tree","\"",")}; 1' ","x_ndseq >x_ndseqtree", sep = "")) # 	where is tree data
system(paste("awk '{sub(/XXX4/,","\"","/home/seb/Documents/PAML/ann_deb_pet_arg_RESULTS/",sp, "/","out","\"",")}; 1' ","x_ndseqtree >codeml.ctl", sep = "")) # 	where is out file

system("rm x*")

######RUN PAML########
#This allow to get estimates of dn and ds. 
system(paste("../codeml codeml.ctl",sep = ""))


#you need the codeml.ctl file + tree.file + sequence file.

#####PAML OUTPUT######
###Parse the output of PAML###
#setwd("~/Desktop/paml44/ann_deb_RESULTS")
ds = as.matrix(read.delim("2NG.dS", header = F))
dn = as.matrix(read.delim("2NG.dN", header = F))
####fixed_contig_names = as.matrix(read.delim("fixed_contig_names.txt", header = F)) # This is the list of fixed contigs. 
orf_paml = as.matrix(read.delim(sp, header = F)) # names


dnds = matrix(0, nrow = (length(ds)/3), ncol = 4)
colnames(dnds) = c("name","dn","ds","dnds")

for(i in 1: (length(ds)/3))
	{
	dnds[i,2] = tail(strsplit(dn[(i*3),1],split = " ")[[1]],1)
	dnds[i,3] =  tail(strsplit(ds[(i*3),1],split = " ")[[1]],1)
	}

	dnds[,4] = as.numeric(dnds[,2]) / as.numeric(dnds[,3])
	dnds[,1] = orf_paml[(c(1:nrow(orf_paml))[orf_paml[,1] == "//end"]+1):nrow(orf_paml),1]
	dnds[,4] = gsub("Inf","NaN", dnds[,4])

write.table(dnds,paste("dnds_",sp, sep = ""), row.names = F, col.names = T, quote = F)

}
#######################
########SANDBOX########
#######################
