

###calculate recombination rates###

setwd("~/Documents/speciation_islands_dec2011") #set up working directory 

system("sort -k 4 reference/combinedPhysicalGeneticMap.out.placedBacs.out | awk '{ print $4}' | uniq -c > reference/tags_total")

system(" awk 'NR > 1' reference/gMap.280x801.30112011.tab | sort -nk 2 -nk 3 | awk '{print $3 }' | uniq >reference/uniq_position")

gMap_uniq = read.table("reference/uniq_position", header = F, stringsAsFactors = F)
gMap = read.table("reference/gMap.280x801.30112011.tab", header = T, stringsAsFactors = F)

tags_total = read.table("reference/tags_total", header = F, stringsAsFactors = F)


gMap_uniq = cbind(1,gMap_uniq,0,0,0)
colnames(gMap_uniq) = c("LG","cm","number_contigs","number_tags","recomb_rate")

for(i in 1:(nrow(gMap_uniq)-1)) #regain the linkage groups
{
gMap_uniq[i+1,1] = gMap_uniq[i,1]
if(gMap_uniq[i,2] <= gMap_uniq[i+1,2]) gMap_uniq[i+1,1] = gMap_uniq[i+1,1] else gMap_uniq[i+1,1] = gMap_uniq[i,1] + 1 
}

#number of contigs per unique position

for(chr in 1:17)
{
chr_temp = gMap[gMap[,2] == chr,]
gMap_uniq_temp = gMap_uniq[gMap_uniq[,1] == chr,]

for(i in 1:nrow(gMap_uniq_temp))
{
	x = chr_temp[,3] %in% gMap_uniq_temp[i,2]

	contig_temp = chr_temp[x == T,1]
	gMap_uniq_temp[i,3] = length(contig_temp)

	for(cc in 1:length(contig_temp))
	{
	tags_total_temp = tags_total[tags_total[,2]  == contig_temp[cc],1]	
	if(length(tags_total_temp) > 0) gMap_uniq_temp[i,4] = gMap_uniq_temp[i,4] + tags_total_temp
	}}
	gMap_uniq[as.numeric(rownames(gMap_uniq_temp)[1]): as.numeric(tail(rownames(gMap_uniq_temp),1)),] = gMap_uniq_temp
	
	print(paste(chr,"of 17 chromosomes",Sys.time()))
	}
	
	
###calculate recombination rate on a sliding window approach (1cM)
###based on a distance of 6kb per tags###
for(chr in 1:17)
	{
	gMap_uniq_temp = gMap_uniq[gMap_uniq[,1] == chr,]
	
		for(i in 1:nrow(gMap_uniq_temp))
			{
			temp = gMap_uniq_temp[(gMap_uniq_temp[,2] >= (gMap_uniq_temp[i,2] - 0.5)) & (gMap_uniq_temp[,2] <= (gMap_uniq_temp[i,2] + 0.5)), ]
			gMap_uniq_temp[i,5] = (max(temp[,2]) - min(temp[,2])) / (sum(temp[,4])*0.006) #centimorgan per megabase
			}
			gMap_uniq[as.numeric(rownames(gMap_uniq_temp)[1]): as.numeric(tail(rownames(gMap_uniq_temp),1)),] = gMap_uniq_temp
	}
gMap_uniq[,5] = ifelse(gMap_uniq[,5] == Inf, NA,gMap_uniq[,5]) #remove the infinite recombination rates
gMap_uniq[,5] = ifelse(gMap_uniq[,5] == 0,NA,gMap_uniq[,5]) #remove the zero recombination rates
	
###	
###plot the recombination rates back to the transcriptome genetic map
###the map I got from the island.R script. 

unique_map_transcript = cbind(map[duplicated(map[,4]) == F,c(1,2,3,4,5)],0)
colnames(unique_map_transcript)[6] = "recomb_rate"

for(chr in 1:17)
{
gMap_temp = gMap_uniq[gMap_uniq[,1] == chr,]
unique_map_transcript_temp = unique_map_transcript[unique_map_transcript[,2] == chr, ]

unique_map_transcript_vector = c(1:nrow(unique_map_transcript))[unique_map_transcript[,2] == chr]


for(i in 1:nrow(unique_map_transcript_temp))
	{
	x = gMap_temp[,2] == unique_map_transcript_temp[i,3] 
	unique_map_transcript_temp[i,6] = gMap_temp[x == T,5]
	}

unique_map_transcript[unique_map_transcript_vector,] = unique_map_transcript_temp
}
mean(unique_map_transcript[,6],na.rm = T) #median recombination rate
median(unique_map_transcript[,6],na.rm = T) #median recombination rate


write.table(unique_map_transcript,"reference/recombination_rates",row.names = F, col.names = T, quote = T)


###
###plot map (both genetic and physical)
###
library(qtl)
unique_map_transcript_all = read.delim("reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
unique_map_transcript_all[is.na(unique_map_transcript_all[,6]),6] =  1
dup = duplicated(unique_map_transcript_all[,4]); unique_map_transcript_all = unique_map_transcript_all[dup == F, ]

write.table(cbind(c(1,"",""),t(cbind(unique_map_transcript_all[,1],unique_map_transcript_all[,2],unique_map_transcript_all[,3]))), "results/annuus_genetic_map.csv",row.names = F,col.names = F, quote = F, sep = ",")
write.table(cbind(c(1,"",""),t(cbind(unique_map_transcript_all[,1],unique_map_transcript_all[,2],unique_map_transcript_all[,7]))), "results/annuus_physical_map.csv",row.names = F,col.names = F, quote = F, sep = ",")
genetic_map = read.cross("csv", file ="results/annuus_genetic_map.csv" )
physical_map = read.cross("csv", file ="results/annuus_physical_map.csv" )

plot.map(genetic_map,physical_map) 

