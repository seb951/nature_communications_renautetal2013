#!/usr/bin/Rscript --slave

setwd("~/Documents/speciation_islands_dec2011/")

###############################
#### NEIGHBOUR JOINING TREE ###
###############################

#all_comparisons_h  = as.matrix(read.delim("mpileup/snp_table_50k_head", header = T, sep = " ")) #head20k SNP
#all_comparisons_t  = as.matrix(read.delim("mpileup/snp_table_50k_tail", header = T, sep = " ")) #tail 40k SNP
#all_comparisons = rbind(all_comparisons_h,all_comparisons_t)

all_comparisons  = as.matrix(read.delim("mpileup/snp_table_2", header = T, sep = " ")) #tail 40k SNP

colnames(all_comparisons) = gsub("X2","2",colnames(all_comparisons))
colnames(all_comparisons) = gsub("X1","1",colnames(all_comparisons))
colnames(all_comparisons) = gsub("PET.","PET-",colnames(all_comparisons), fixed = T)

four_species = read.delim("reference/4_species", header = T, stringsAsFactors = F)
ind = colnames(all_comparisons)

for(i in 4:length(ind))
	{ind[i] = four_species[regexpr(colnames(all_comparisons)[i],four_species[,1]) > 0,2]} # there is a warning because PL109 is present twice. but it doesnt really matter.

 colnames(all_comparisons)[4:ncol(all_comparisons)] = paste(ind[4:ncol(all_comparisons)],colnames(all_comparisons)[4:ncol(all_comparisons)],sep = "_")
 
################
### TRIMMING ###
################

###kick out SNP which have more than 10% XX data missing
too_much_missing = c(1:nrow(all_comparisons))

	for(i in 1:nrow(all_comparisons))
		{
		n_ind = ncol(all_comparisons)-3
		too_much_missing[i] = length(c(1:n_ind)[grepl("XX",all_comparisons[i,4:ncol(all_comparisons)]) == T]) / n_ind
		}

all_comparisons_2 = all_comparisons[too_much_missing < 0.1,]
all_comparisons = NULL
##############################
### trim based on 2pq (Ht) ###
##############################
all_comparisons_3 = 1

#counting the alleles.
	counter_ACGTX = matrix(0,nrow = nrow(all_comparisons_2), ncol = 5)
	colnames(counter_ACGTX) = c("A","C","G","T","X")
	ht = c(1:nrow(all_comparisons_2))
	n_ind = ncol(all_comparisons_2)
	
	for(i in 1:nrow(all_comparisons_2))
		{
		x = paste(all_comparisons_2[i,4:n_ind], collapse = "")
		xx = strsplit(x, split = "")
		
		counter_ACGTX[i,1] = length(c(1:((n_ind-3)*2))[grepl("A", xx[[1]])])
		counter_ACGTX[i,2] = length(c(1:((n_ind-3)*2))[grepl("C", xx[[1]])])
		counter_ACGTX[i,3] = length(c(1:((n_ind-3)*2))[grepl("G", xx[[1]])])
		counter_ACGTX[i,4] = length(c(1:((n_ind-3)*2))[grepl("T", xx[[1]])])
		counter_ACGTX[i,5] = length(c(1:((n_ind-3)*2))[grepl("X", xx[[1]])])
		p = sort(counter_ACGTX[i,1:4])[4] 
		q = sort(counter_ACGTX[i,1:4])[3] 
		al = sum(counter_ACGTX[i,1:4])
		ht[i] = 2 * (p / al) * (q / al)
		}
	all_comparisons_3 = all_comparisons_2[ht > 0.05,]


###################################
### trim based on Ho (paralogs) ###
###################################

all_comparisons_4 = list(1)

ho  = rep(0,nrow(all_comparisons_3))

for(i in 1:nrow(all_comparisons_3))
	{
		a1 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],1,1)
		a2 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho[i] = (ho[i] + 1) #count the heterozygotes
		}
	}
	ho = ho / (ncol(all_comparisons_3)-3) # observed heterozygosity
	all_comparisons_4 = all_comparisons_3[ho < 0.6,]

###################################
### create a single consensus of both alleles. ###
###################################	
	
col = ncol(all_comparisons_4)
all_comparisons_4_nj_preformat = all_comparisons_4

for(i in 1:nrow(all_comparisons_4))
	{
	temp1 = all_comparisons_4[i,4:col]
	temp2 = temp1
	for(j in 1:length(temp2))
		{
		if(substring(temp1[j],1,1) == substring(temp1[j],2,2)) temp2[j] = substring(temp1[j],1,1)
		if(substring(temp1[j],1,1) == "X" & substring(temp1[j],2,2) != "X" ) temp2[j] = substring(temp1[j],2,2)
		if(substring(temp1[j],2,2) == "X" & substring(temp1[j],1,1) != "X" ) temp2[j] = substring(temp1[j],1,1)
	
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "C") | (substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "A"))) temp2[j] = "M"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "A"))) temp2[j] = "T"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "A"))) temp2[j] = "W"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "C"))) temp2[j] = "S"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "C"))) temp2[j] = "Y"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "G"))) temp2[j] = "K"
		
		if(temp2[j] == "X") temp2[j] = "N"
		}
	if((i %% 500) == 0) print(paste(i, "of", nrow(all_comparisons_4),Sys.time()))
	all_comparisons_4_nj_preformat[i,4:col] = temp2
	}

all_comparisons_4_nj_preformat = all_comparisons_4_nj_preformat[,-c(95,99)] # these two hicks are out
write.table(all_comparisons_4_nj_preformat,"all_comparisons_4_nj_preformat", row.names = F, col.names = T, quote = T)	

############
####NJ TREE###
############
	
library(ape)
all_comparisons_4_nj_preformat = as.matrix(read.delim("all_comparisons_4_nj_preformat", header  = T, sep = " "))

all_comparisons_4_nj_format = t(tolower(all_comparisons_4_nj_preformat))[4:ncol(all_comparisons_4_nj_preformat),]
	
colnames(all_comparisons_4_nj_format) = apply(cbind(all_comparisons_4_nj_preformat[,1],as.numeric(all_comparisons_4_nj_preformat[,2])),1,paste,collapse = "_")
	
all_comparisons_4_nj_format_1 = as.DNAbin(all_comparisons_4_nj_format) # transform DNA matrix to DNAbin object
all_comparisons_4_nj_format_2 = dist.dna(all_comparisons_4_nj_format_1, pairwise.deletion = T) # distance matrix

#all_comparisons_4_nj_format_2[(all_comparisons_4_nj_format_2) == Inf] = 10
#all_comparisons_4_nj_format_2[is.na(all_comparisons_4_nj_format_2)] = 10

all_comparisons_4_nj_format_3 = nj(all_comparisons_4_nj_format_2) # neighbour joining tree

###################
###PLOTTING TREE###
###################
edge.col = cbind(all_comparisons_4_nj_format_3$tip.label,0)

#plot.phylo
#tips color
#all_comparisons_pop = cbind(all_comparisons_pop,c(rep("darkred",length(ann)),rep("darkblue",length(pet)),rep("darkgreen",length(deb)), rep("black",length(arg)) ))


for(i in 1:length(all_comparisons_4_nj_format_3$tip.label))
	{
	if(regexpr("ann_",edge.col[i,1]) > 0) edge.col[i,2] = "darkred"
	if(regexpr("pet_",edge.col[i,1]) > 0) edge.col[i,2] = "darkblue"
	if(regexpr("deb_",edge.col[i,1]) > 0) edge.col[i,2] = "darkgreen"
	if(regexpr("arg_",edge.col[i,1]) > 0) edge.col[i,2] = "black"
	}
	
all_comparisons_4_nj_format_3$tip.label = gsub("ann_|pet_|deb_|arg_|inv_|wee_|.white|tex_|pra_|.454Reads","",all_comparisons_4_nj_format_3$tip.label); all_comparisons_4_nj_format_3$tip.label[13] = "arg_arg1820.white"

###boot the phylogeny

 f = function(x) nj(dist.dna(x, pairwise.deletion = T))
 tr <- f(all_comparisons_4_nj_format_1)
b = boot.phylo(tr, all_comparisons_4_nj_format_1, f, quiet = F, B = 1000); b = round(b /10)
bb = ifelse(b >80,b,"")
t = c("phylogram", "cladogram", "fan", "unrooted", "radial" )
	
#plot(all_comparisons_4_nj_format_3, type = "fan",font = 2, cex = 1, tip.col = edge.col[,2], lwd = 4, edge.width = 4  ) #plotting tree. (plot.phylo) 
plot(all_comparisons_4_nj_format_3, type = t[1],font = 2, cex = 0.6, tip.col = edge.col[,2], edge.width = 4 , label.offset = 2, node.pos = 1 ,use.edge.length = F) #plotting tree
nodelabels(bb,frame = "n")
dev.print(device=svg, "~/Desktop/njtree.svg", onefile=FALSE)
dev.off()
	
#


plot(all_comparisons_4_nj_format_3, type = "fan",font = 2, cex = 1, tip.col = edge.col[,2], lwd = 4, edge.width = 4  ) #plotting tree. 
	
##############
####NEXUS FILE FOR SPLITSTREE###
##############	
library(ape)
all_comparisons_4_nj_preformat = as.matrix(read.delim("all_comparisons_4_nj_preformat", header  = T, sep = " "))
all_comparisons_4_nj_preformat = all_comparisons_4_nj_preformat[c(1:10000),]
sequences = c(1:(ncol(all_comparisons_4_nj_preformat) -3))
colnames(all_comparisons_4_nj_preformat) = gsub(".454Reads","",colnames(all_comparisons_4_nj_preformat),fixed = T)
colnames(all_comparisons_4_nj_preformat) = gsub(".white","",colnames(all_comparisons_4_nj_preformat),fixed = T)
colnames(all_comparisons_4_nj_preformat) = gsub("arg_|deb_|ann_|pet_","",colnames(all_comparisons_4_nj_preformat))
for(j in 4:ncol(all_comparisons_4_nj_preformat))
{
temp = colnames(all_comparisons_4_nj_preformat)[j]
temp = paste(temp,paste(all_comparisons_4_nj_preformat[,j], collapse = ""), sep = "     ")
sequences[j-3] = temp
}

nrow(all_comparisons_4_nj_preformat)
ncol(all_comparisons_4_nj_preformat)

all_helianthus_speciation_island = c("#NEXUS","BEGIN taxa;","DIMENSIONS ntax=106;","TAXLABELS",colnames(all_comparisons_4_nj_preformat)[4:ncol(all_comparisons_4_nj_preformat)],";","END;","BEGIN characters;","DIMENSIONS nchar=10000;","FORMAT","datatype=DNA","missing=N","gap=-","symbols=\"A C G T\"","labels","interleave",";","MATRIX",sequences, ";","END;")

write.table(all_helianthus_speciation_island,"results/all_helianthus_speciation_island_10k.nex", sep = "", quote = F, row.names = F,col.names = F)

#############
##SAND AND SEX AND ROCK AND ROLL# 



