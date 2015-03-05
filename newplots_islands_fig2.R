
setwd("~/Documents/speciation_islands_dec2011/")

###
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",]
map = cbind(map,0,0,0,0)
colnames(map) = c("name","LG","map_centiMorgan","unique_position","nbnuc","dn","ds")

#unique position for the markers#
max_position = cbind(c(1:17), 0,0) # the size of each chromosome and the cumulative number of centimorgans.

for(i in 1:17)
	{
	max_position[i,2] = max(as.numeric(map[as.numeric(map[,2]) == i,3]))
	max_position[i,3] = sum(max_position[1:i,2])-max_position[i,2]
	}


###genome scans
dev.new(width=18, height=6)
par(mfrow = c(4,1))

for(a in c(1,2,6,5))
{
################
### PLOTTING ###
################
#if(a == 7) c1 =  read.table("results_6species/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)
if(a == 1) {c1 = read.table("results/ann_pet_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_PET"}
if(a == 2) {c1 =  read.table("results/ann_deb_fst_map_cluster2", header = T, stringsAsFactors = F);main = "ANN_DEB"}
if(a == 3) {c1 =  read.table("results/ann_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "ANN_ARG"}
if(a == 4) {c1 =  read.table("results/pet_deb_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_DEB"}
if(a == 5) {c1 =  read.table("results/pet_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "PET_ARG"}
if(a == 6) {c1 = read.table("results/deb_arg_fst_map_cluster2",  header = T, stringsAsFactors = F);main = "DEB_ARG"}

color = c("#f8756b","#7aae00","#7aae00","#00bdc2","#c77aff","#00bdc2")
###fst outlier regions
par(mar = c(3,4,1.5,2))
for(i in 1:17) #17 linkage groups
{
	cc1= c1[c1[,2] == i, ] # plot one linkage group at a time.

	if(i != 17) cc2 = c1[c1[,2] == (i+1), ] else {cc2 = cc1[cc1[,1] == max(cc1[,1]),]; cc2 = rbind(cc2,cc2)}
	if(i == 1) plot(0,1,  type = "l", xlim = c(-20,max(c1[,1])+21), ylim = c(0,1.09), xaxt = "n", yaxt = "n", lwd = 1,col = "black", main = main, xlab = "", ylab = "",xaxs = "i")
	if(i == 1) rect(0,0,min(cc2[cc2[,3]>0,1]),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))
	if(i != 1) rect(min(cc1[cc1[,3]>0,1]),0,min(cc2[cc2[,3]>0,1]),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))

	points(cc1[cc1[,3]>0,1],cc1[cc1[,3] > 0 ,3], type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(0,1.09), lwd  = 1, col =  color[a]) #just the points OR next line...
	points(lowess(cc1[cc1[,3] > 0 ,3]~cc1[cc1[,3]>0,1], f = 1 / 100 ), type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(0,1.09), lwd  = 4, col =  color[a]) # some amount of smoothing, f = 1 / 100

}
points(c1[c1[,16]<0.001,1],rep(1.09,length(c1[c1[,16]<0.001,1])), col = color[a], pch = 19) # significant peaks.
axis(1, at = max_position[,3], label = c(1:17), line = 0.5) #x axis
axis(2, at = c(0,0.2,0.4,0.6,0.8,1),  label = c(0,0.2,0.4,0.6,0.8,1),line= 0.5, las = 1) #y axis
}

dev.print(device=svg, "~/Desktop/fst.svg", onefile=FALSE)
dev.off()



