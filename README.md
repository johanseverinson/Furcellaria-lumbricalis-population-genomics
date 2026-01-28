# Furcellaria-lumbricalis-population-genomics
Scrits and data to accompany the manuscript "Drifting Furcellaria lumbricalis Exhibits Genetic Differentiation and Increased Clonality Compared to Attached Conspecifics" by Severinson, J., Moksnes, P-O. and Jahnke, M.

In this study, we compare attached and free-drifting forms of the perennial red algae Furcellaria lumbricalis, to investigate their reproductive mode, stability and population structure. Samples were collected from four sites and sequencing was performed using 2bRAD sequencing.

Fasta files and sample information can be found at NCBI BioProject: PRJNA1387154. The naming convention for fasta files are: site-FL-type-number, where sites are Askö (A), Kalvöfjorden (K), Nordön (N) and Ryskärsfjorden (R), while types are attached to hard substrate (A) and drifting on soft bottom sediments (M and T). For example, an attached sample from Askö is called A-FL-A-5. You will see this naming convention also used throughout scripts (as opposed to the naming in the manuscript).

For contact: johan.severinson@gu.se

Files:

IBS_matrices.zip  
IBS matrices produced by ANGSD.

Protocol_for_Illumina_2b_RAD_Ries_2023.pdf  
Lab protocol for library preparations for 2b-RAD.

admixture_all_sites.zip  
Files produced by NGSadmix for evaluation of K and for admixture of all sites together.

admixture_each_site.zip  
Files produced by NGSadmix for admixture of sites separately.

bayescan_files.zip  
Files produced by BayeScan for identifying outlier loci. Needed are also vcf files, can be converted from beagle files in calculating_metrics_files.zip.

bioinfo_pipeline_Furcellaria2026.sh  
Bioinformatics pipeline used to produce bam files and run bioinformatics on cluster.

calculating_metrics_files.zip  
Files produced by ANGSD, PCAngsd and realSFS to calculate R, B, D*, Ho, He, Fis and pi.

meta_data.zip  
Meta data for all samples or samples retained after clone correction.

plot_R.r  
Function to plot BayeScan results.

plot_admixture_v5_function.R  
Function to plot admixture.

popgen_analyses_R.R  
R script for performing clone correction, PCoA, admixture, test isolation by distance, find outlier loci from BayeScan and the selective sweep, calculate and test diversity metrics and produce maps.

selective_sweep_files.zip  
Files produced by PCAngsd for selective sweep.
