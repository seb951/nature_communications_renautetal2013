
res = seq(1,10,by = 0.1)
res = cbind(res,0)

for(i in seq(1,10,by = 0.1))
{
res[res[,1] == i ,2] = chisq.test(c(59*i,53*i,44*i,54*i))$p.value
}


library(pwr)
pwr.chisq.test()

ch = chisq.test(c(59,53,44,54))
effect_size = sqrt((ch$statistic*ch$statistic)   / (210*(4-1)). ) ###kramer's V

pwr.chisq.test(w = effect_size,N = 210, df = 3, sig.level=0.05, power = NULL)




###
###how many more islands would I need to find to see a significant difference in size and number?
###
res = NULL

for(sim in 1:100)
{

for(n in seq(0,5,by = 0.1))
{
one = (size_all_all[size_all_all[,1] == 1,])
two = (size_all_all[size_all_all[,1] == 2,])
fiv = (size_all_all[size_all_all[,1] == 5,])
six = (size_all_all[size_all_all[,1] == 6,])

one_s = rbind(one,one[sample(c(1:nrow(one)),ceiling(nrow(one) * n),replace = T),])  #add samples according to n 
two_s = rbind(two,two[sample(c(1:nrow(two)),ceiling(nrow(two) * n),replace = T),])  #add samples according to n 
fiv_s = rbind(fiv,fiv[sample(c(1:nrow(fiv)),ceiling(nrow(fiv) * n),replace = T),])  #add samples according to n 
six_s = rbind(six,six[sample(c(1:nrow(six)),ceiling(nrow(six) * n),replace = T),])  #add samples according to n 

size_all_all_s = rbind(one_s,two_s,fiv_s,six_s)
k_pval = kruskal.test(x = size_all_all_s[,6],g = as.factor(size_all_all_s[,1]))$p.value #island size
c_pval = chisq.test(c(nrow(one_s),nrow(two_s),nrow(fiv_s),nrow(six_s)))$p.value #island number

res = rbind(res, c(n,c_pval,k_pval))
}
if(sim %% 10 == 0) print(paste(sim,"done, time is:",Sys.time()))
}

res_final = cbind( seq(0,5,by = 0.1) * 100,0,0); n = seq(0,5,by = 0.1)
for(i in 1:length(n)) {res_final[i,2] = mean(res[res[,1] == n[i],2]); res_final[i,3] = mean(res[res[,1] == n[i],3])}


plot(res_final[,1],res_final[,2],lwd = 10,col = "darkred")
points(res_final[,1],res_final[,3],lwd = 10,col = "darkblue", ylim = c(0,0.6))
abline(h = 0.05, lty = 2)
dev.print(device=svg, "~/Desktop/power2.svg", onefile=FALSE)
dev.off()



