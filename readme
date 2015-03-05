###READ ME ###
###note from MARCH 2015: I decided to add all my scripts used for this paper.
###Renaut, S., et al. "Genomic islands of divergence are not affected by geography of speciation in sunflowers." Nature communications 4 (2013): 1827.
###however, use at your own risks. These scripts may not be properly documented. They also need various bioinformatics softwares to run, some of which may be outdated today. (Most bioinfo software date from late 2011, 2012).
###Sebastien Renaut 2015###

###Step 1.
bwa_samtools.R
Aligns using bwa and generated a bam files using samtools


###Step 2. 
mpileup_2.R
Generates a large SNP table using mpileup and bcftools. 
This is then split into 6 different tables for each comparison to generate six files "snp_table_ann_pet"


###Step 3. 
snp_parser_4species_v3.R
This filters the snp_table_* files and calculates Fst to generate six "*_comp4" files.


###Step 4.
PAML_consensus_calling_from_SNPdata.R
Takes in the "*_comp4" files to generate "dnds_*" files.


###Step 5. 
island.R
Takes in the "*_comp4" and "dnds_*" files. , then adds the map information to each locus. Generates the "*_fst_map" files
Then calculates the clustering of divergent sites and generates six "*_cluster2" files. 


###Step 6
Moran.R
takes in the "*_fst_map" files to calculate measures of spatial autocor per chromosome to generate "*_moran" files.

###Step 7
plotting.R
takes in the "*_cluster2" files and plot histogram, number of islands, etc...
 
