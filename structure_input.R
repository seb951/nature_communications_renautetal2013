#!/usr/bin/Rscript --verbose

setwd("~/Documents/speciation_islands_dec2011/")

snp_table  = read.delim("mpileup/snp_table_2", header = T, sep = " ", stringsAsFactors = F, nrows = 10000)

snp_table[,2] = gsub(" ","",snp_table[,2])

#snp_table = snp_table[,-c(27,95,99)]#kick out individuals 81 (PL109.white) and 124 (PI586932b) because they didnt sequence good. In addition, kick out btm15.2, btm30.4 and ARG1820 which are not what they claim they are. 
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


ind[regexpr("ann",ind) > 0] = 1   
ind[regexpr("pet",ind) > 0] = 2
ind[regexpr("deb",ind) > 0] = 3
ind[regexpr("arg",ind) > 0] = 4
structure_matrix = matrix(sort(rep(c(1:108),2)), nrow = 108*2, ncol = 1)
structure_matrix = cbind(structure_matrix,0)
structure_matrix[seq(1,108*2,by = 2),2] = ind[4:111]
structure_matrix[seq(2,108*2,by = 2),2] = ind[4:111]




rand_1000 = sample(c(1:nrow(snp_table)),1000)  # random set of 200 nuclear loci
for(i in rand_1000)
{
	hf = cbind(rep(c(1:108),2),c(substring(as.character(snp_table[i,4:111]),1,1),substring(as.character(snp_table[i,4:111]),2,2))) #population info
	hf = hf[order(as.numeric(hf[,1])),]

	structure_matrix = cbind(structure_matrix,hf[,2])
}


for(i in 1:216)
	{
		for(j in 3:ncol(structure_matrix))
		{
	if(structure_matrix[i,j] == "0") structure_matrix[i,j] = "-9"
	if(structure_matrix[i,j] == "X") structure_matrix[i,j] = "-9"

	if(structure_matrix[i,j] == "") structure_matrix[i,j] = "-9"

	if(structure_matrix[i,j] == "A") structure_matrix[i,j] = 1
	if(structure_matrix[i,j] == "C") structure_matrix[i,j] = 2
	if(structure_matrix[i,j] == "T") structure_matrix[i,j] = 3
	if(structure_matrix[i,j] == "G") structure_matrix[i,j] = 4
		}
	}


#write.table(structure_matrix,"/usr/share/structure/speciation_islands/speciation_islands_1000_108.input", row.names = F, col.names = F, quote = F)
write.table(paste(c(1:length(rand_1000)),collapse = "\t"),"/usr/share/structure/speciation_islands/speciation_islands_1000_108.input", row.names = F, col.names = F, quote = F)

cat(apply(structure_matrix, 1,paste, collapse = "\t"),file = "/usr/share/structure/speciation_islands/speciation_islands_1000_108.input",append = T, fill = 1)

#run this shit in /usr/share/structure. make sure you update the mainparam file

###############
###Run structure
###############

setwd("/usr/share/structure/speciation_islands")


cpu = 5 #how many cpu (max)

#all
for(i in 1:20)
#for(i in 26:100)
{

system("ps -u seb | grep 'structure' | wc -l >prog")
print("pogo")
Sys.sleep(ifelse(read.table("prog") < cpu,0, (read.table("prog") - cpu + 4)^4    ))

for(j in 1:6)
{
system("ps -u seb | grep 'structure' | wc -l >prog")

Sys.sleep(ifelse(read.table("prog") < cpu,0, (read.table("prog") - cpu + 3)^4    ))

command_1 = paste("nohup ./../structure -i speciation_islands_1000_108.input -e extraparams -m mainparams_", j," -o output/outfile_k",j, "_",i," >output/log",j,"_",i, " &", sep = "")

system(command_1)
}
print(paste("repetition_",i," the time is:", Sys.time(), sep = ""))
}
################
###structureHarverster 
################
system("./../harvester/structureHarvester.py --dir output/ --out harvest_6 --evanno --clump")





#read structure input#

#setwd("/usr/share/structure/")

#k_out = read.delim("outfile_5000_k4_f", header = T, stringsAsFactors = F); k = 4

#clusters = k_out[(grep("90% probability intervals", k_out[,1])+1): (grep("Estimated Allele Frequencies in each cluster",k_out[,1])-1 ),1]
#clusters_m = matrix(0,nrow = length(clusters), ncol = k +1)

#for(i in 1:length(clusters))
#	{
#	clusters_m[i,2:(k+1)] = as.numeric(strsplit(clusters[i]," ")[[1]][strsplit(clusters[i]," ")[[1]] != ""][(6:(k+5))])
#	}

#clusters_m[,1] = populations

#populations_2 = rep(0,29)

#for(i in 1:length(populations))
#	{
#	x = tex_deb_pet[,1] %in% colnames(tex_deb_pet_fst)[(i+3)]
#	if(length(x[x == T]) == 1) populations_2[i] = tex_deb_pet[x == T,1] 
#	}

#clusters_m =  clusters_m[order(populations),][-c(25),]

#par(mfrow = c(3,1))
#barplot(t(clusters_m[,2:(k+1)]), col = c("red","green","blue","yellow","black"), space = 0, border = NA)
#text(c(0:27),0,populations_2[order(populations)][-c(25)], srt = 90, pos = 4)
#dev.print(device=svg, "structure_k3-5.svg", onefile=FALSE)
#dev.off()

################
###clummpppp
################
#summary matrix
#delta_k = matrix(0,ncol = 6, nrow = 5)

###delta K
#k2_res_prob = k2_res[grep("Estimated Ln Prob of Data",k2_res)]
#k2_res_prob = as.numeric(gsub("Estimated Ln Prob of Data   = ","",k2_res_prob))

#k3_res_prob = k3_res[grep("Estimated Ln Prob of Data",k3_res)]
#k3_res_prob = as.numeric(gsub("Estimated Ln Prob of Data   = ","",k3_res_prob))

#k4_res_prob = k4_res[grep("Estimated Ln Prob of Data",k4_res)]
#k4_res_prob = as.numeric(gsub("Estimated Ln Prob of Data   = ","",k4_res_prob))

#k5_res_prob = k5_res[grep("Estimated Ln Prob of Data",k5_res)]
#k5_res_prob = as.numeric(gsub("Estimated Ln Prob of Data   = ","",k5_res_prob))

#delta_k[1,1] = mean(k2_res_prob); delta_k[1,2] = sd(k2_res_prob)
#delta_k[2,1] = mean(k3_res_prob); delta_k[2,2] = sd(k3_res_prob)
#delta_k[3,1] = mean(k4_res_prob); delta_k[3,2] = sd(k4_res_prob)
#delta_k[4,1] = mean(k5_res_prob); delta_k[4,2] = sd(k5_res_prob)


#write.table(delta_k, "delta_k.csv")








