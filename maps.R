library(maps)
library(mapproj)

setwd("~/Documents/speciation_islands_dec2011/")

latlong = as.matrix(read.delim("reference/latlong_final.csv", header = F, sep = "\t"))

#change alberta to go further south

latlong[3,6:7] = c("-113","49.25")



######################
####	ALL POP #########
######################

map('usa', proj='bonne', param=35)
map("state", lty = 2, add = T, col = "#00000080",boundary = T,proj='bonne', param=35)

title("Distribution of populations", line = 2)


for(i in 1:length(latlong[,2]))
	{
		if(regexpr("ann",latlong[i,2]) > 0) {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkred", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1), mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])),cex =0.7, col = "darkred",font = 2, labels = latlong[i,3])}
		if(regexpr("pet",latlong[i,2]) > 0)  {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkblue", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1),mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])),cex =0.7, col = "darkblue",font = 2, labels = latlong[i,3])}
		if(regexpr("arg",latlong[i,2]) > 0)  points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "black", pch = 21, lwd = 4) #;text(adj = c(-0.1,-0.1),as.numeric(latlong[i,4]), as.numeric(latlong[i,5]), cex =0.5, col = "black",font = 2, labels = latlong[i,1])}
		if(regexpr("deb",latlong[i,2]) > 0)  points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkgreen", pch = 21, lwd = 4)  #;text(adj = c(-0.1,-0.1),as.numeric(latlong[i,4]), as.numeric(latlong[i,5]), cex =0.5, col = "darkgreen",font = 2, labels = latlong[i,1])}
	}
	
	
legend(-123,22, legend = c(expression(~italic(H.annuus)~"(18 individuals)"), expression(~italic(H.petiolaris)~"(24 individuals)"), expression(~italic(H.argophyllus)~"(22 individuals)"),expression(~italic(H.debilis)~"(10 individuals)")),pch = 21, col = c("darkred", "darkblue","black","darkgreen"), pt.lwd = 4)

	
dev.print(device=svg, "~/Desktop/USA.svg", onefile=FALSE)
dev.off()


#####################
#TEXAS ONLY #########
#####################
map('state', 'texas',proj='bonne', param=35)
map('county', 'texas',col = "#00000080", lty = 2, add = T, boundary = F,proj='bonne', param=35)

title("Distribution of populations TEXAS ", line = 2)


for(i in 1:length(latlong[,2]))
	{
		if(regexpr("ann",latlong[i,2]) > 0) {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkred", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1), mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])),cex = 1, col = "darkred",font = 2, labels = latlong[i,3])}
		if(regexpr("pet",latlong[i,2]) > 0)  {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkblue", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1),mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])),cex = 1, col = "darkblue",font = 2, labels = latlong[i,3])}
		if(regexpr("arg",latlong[i,2]) > 0)  {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "black", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1),mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), cex =1, col = "black",font = 2, labels = latlong[i,3])}
		if(regexpr("deb",latlong[i,2]) > 0)  {points(mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), col = "darkgreen", pch = 21, lwd = 4);text(adj = c(-0.1,-0.1),mapproject(as.numeric(latlong[i,6]), as.numeric(latlong[i,7])), cex = 1, col = "darkgreen",font = 2, labels = latlong[i,3])}
	}


dev.print(device=svg, "~/Desktop/texas.svg", onefile=FALSE)
dev.off()

