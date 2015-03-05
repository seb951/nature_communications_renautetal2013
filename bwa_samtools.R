#!/usr/bin/Rscript --verbose

###Re do alignments of all fastq files for ALL argophyllus, petiolaris, debilis and annuus. 

###
###Step 1 set up working directory and individuals of interest##
###
setwd("/SciBorg/array0/renaut/speciation_islands_individuals")

individuals = read.delim("individuals/4_species_fuck_ups", header = T, stringsAsFactors = F)
individuals = cbind(individuals,0)

for(i in 1:nrow(individuals)) {individuals[i,5] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}

###
### Step 2. Index reference
###
#system("bwa index reference/HA412_trinity_noAltSplice_400bpmin.fa")
#system("samtools faidx reference/HA412_trinity_noAltSplice_400bpmin.fa")

###
###Step 3. BWA ALIGNMENTS create fai files.
###

for(i in 1: nrow(individuals)
#for(i in c(1:4))
{
bwa_aln1 = bwa_aln2 = "ls"

	if((individuals[i,3] == "ill") & (individuals[i,4] == "sanger")) # new illumina quality format
{bwa_aln1 = paste("bwa aln -t 14 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 14 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

	if((individuals[i,3] == "ill") & (individuals[i,4] == "Ill1.3")) # old illumina quality format
{bwa_aln1 = paste("bwa aln -t 14 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 14 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

system(bwa_aln1)
system(bwa_aln2) 

}

###
###Step 3. BWA ALIGNMENTS create sam files.
###
for(i in 1: nrow(individuals)
#for(i in c(1:4))
{
bwa_sampe = bwasw = "ls"

if(individuals[i,3] == "ill")
bwa_sampe = paste("bwa sampe reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai"," alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai ",sub(".fq","_1.fq",individuals[i,1], fixed = T)," ", sub(".fq","_2.fq",individuals[i,1], fixed = T)," >alignments/",individuals[i,5], ".sam", sep = "")

if(individuals[i,3] == "454") bwasw = paste("bwa bwasw -t 12 reference/HA412_trinity_noAltSplice_400bpmin.fa ",individuals[i,1]," -f alignments/", individuals[i,5],".sam", sep = "")

#if(i %% 4 != 0) system(paste("nohup ",bwa_sampe, " &", sep = "")) else {Sys.sleep(300);system(bwa_sampe)} #run in background 3 jobs. 

system(bwa_sampe)
system(bwasw) #run in foreground. 

}

###
###Step 4. SAMTOOLS VARIANT AND CONSENSUS CALLING
###
for(i in 1: nrow(individuals)
{
sam_view = paste("samtools view -bS -o alignments/",individuals[i,5], ".bam", " alignments/",individuals[i,5], ".sam", sep = "")
#if(i %% 4 != 0) system(paste("nohup ",sam_view, " >view_",i, "_log &", sep = "")) else {Sys.sleep(300);
system(sam_view)
}
#
for(i in 1: nrow(individuals)
{
sam_sort = paste("samtools sort alignments/",individuals[i,5], ".bam ","alignments/",individuals[i,5], ".sorted",sep = "")
#if(i %% 4 != 0) system(paste("nohup ",sam_sort, " >view_",i, "_log &", sep = "")) else {Sys.sleep(300);
system(sam_sort)
}
#index bam files
for(i in 1: nrow(individuals)
{
sam_index = paste("samtools index alignments/",individuals[i,5],".sorted.bam",sep = "")
#if(i %% 10 != 0) system(paste("nohup ",sam_index, " >index_",i, "_log &", sep = "")) else {Sys.sleep(30);
system(sam_index)
}


#index the reference file per gene chunks so that you can parralelize the alignments#
#system("vcfutils.pl splitchr -l 100000 reference/HA412_trinity_noAltSplice_400bpmin.fa.fai >reference/split")
###mpileup###
split = read.table("reference/split", stringsAsFactors = F)
cpu = 12 # how many cpus can you afford

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
system("echo -n '' >mpileup/all_variants.raw.vcf") #final file
system("cat mpileup/variants* >> mpileup/all_variants.raw.vcf") # cat everything
system("grep '#' -v mpileup/all_variants.raw.vcf >mpileup/all_variants.clean.vcf") #tidy up
system("rm mpileup/variants*")
system("rm mpileup/all_variants.raw.vcf")








