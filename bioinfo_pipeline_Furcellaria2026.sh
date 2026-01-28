# 2bRAD pipeline log for Furcellaria lumbricalis
# Edited by Johan Severinson 2026
# Based on version by Stefanie Ries 2023

###############################################
### 	     	Rackham		                      ###
###############################################
# user guides
# https://uppmax.uu.se/support-sv/user-guides/slurm-user-guide/
# https://uppmax.uu.se/support-sv/user-guides/rackham-user-guide/

uquota #file system usage on Dardel
projinfo #CPU hrs
jobinfo -u xsevej
finishedjobinfo -j
squeue -u xsevej
projsummary[project id] #overview

#log in to Rackham
ssh -X xsevej@rackham.uppmax.uu.se

# log into Dardel (new cluster for 2025)
Add current IP-adress at: https://loginportal.pdc.kth.se/ by clicking Add address under Dardel.
Activate network identity manager by renewing credentials for jseve@NADA.KTH.SE, using the PDC kerberos password
Load the Dardel profile/session in PuTTY and hit open.

# Get interactive node with salloc. Once you get the node, connect with: ssh <node>
salloc --nodes=1 -t 8:00:00 -A naiss2025-22-1128 -p shared
ssh XXX # Connect to assignet allocation.

# Load modules on Dardel
module load bioinfo-tools
ml PDC
ml nano

# From windows-created .txt files, use this command to remove hidden new-line characters:
sed -i -e "s/\r//g" file.txt

# Transfer files to and from Dardel. Input this command into cmd on Windows. Reverse local and remote file depending on it sending to or from.
"C:\Program Files (x86)\PuTTY\pscp.exe" -load dardel jseve@dardel.pdc.kth.se:/cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/Fastqc_Results/*html C:\test

###########
## There's a problem using wildcards (*) for scp on a Mac. In short OsX thinks * means something else.
## You can easily fix this by putting quotation marks around the filepath.
## scp 'usr@rackham.uppmax.uu.se:/dir/path/*.gz' .
## More details here: https://unix.stackexchange.com/questions/27419/how-to-use-wildcards-when-copying-with-scp (bearbeitet) 
###########

cd /cfs/klemming/projects/supr/naiss2023-23-410 # my storage

# info given by UPPMAX verify the integrity of the downloaded files with checksums
# for each checksum file, go to the folder where it is located and verify the checksums within
#mv checksums.md5 one folder up
cp checksums.md5 ..
md5sum -c checksums.md5

# -n nodes
# -t h:mm:ss
# -A project
# -p partitioning

# for a job submission

#!/bin/bash -l
#SBATCH -A snic2023-22-797
#SBATCH -p node -n 64
#SBATCH -t 2:00:00
#SBATCH -J jobname
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

# to submit
# sbatch <filename>


#####################################
###     DOWNLOAD 2bRAD SCRIPTS    ###
#####################################

# downloading and installing all 2bRAD scripts in $HOME/bin (or change to whatever directory you want)
cd
mkdir bin
cd ~/bin
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
#git clone https://github.com/z0on/2bRAD_GATK.git
# move scripts to ~/bin from sub-directories
#mv 2bRAD_GATK/* .
mv 2bRAD_denovo/* .
# remove now-empty directory
rm -rf 2bRAD_denovo
rm -rf 2bRAD_GATK

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl
chmod +x *.py
chmod +x *.R

# adding ~/bin to your $PATH
cd
nano .bashrc
# paste this where appropriate (note: .bashrc configuration might be specific to your cluster, consult your sysadmin if in doubt)
export PATH=$HOME/bin:$PATH

# Ctl-o, Ctl-x  (to save and exit in nano)
# log out and re-login to make sure .bashrc changes took effect

# does it work?
# try running a script from $HOME:
cd
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong

#####################################
### DOWNLOAD READS FROM BASESPACE ###
#####################################

#organize fastqs
mkdir fastqs
cd fastqs

# copy them like this
rsync --progress /proj/naiss2023-23-410/onlyr1/230712_A00605_0587_AH77HGDRX3/Sample_WF-3727*/* /proj/naiss2023-23-410/fastqs
# For Dardel 2025
rsync --progress /cfs/klemming/projects/supr/naiss2023-23-410/YC-4275/250430_A00605_0755_AHYFGJDRX5/Sample_YC-4275*/* /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025
# Count number of files
ls -1q YC-4275* | wc -l # 424
# Remove R2
rm /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*_R2_*

ls ./*L001_R1_001.fastq.gz | wc -l # 212

# Below only relevant if you have reference genome, otherwise go to concatenate lane duplicates
############################
####### INDEX GENOME #######
############################

# ran this interactively
module load bioinfo-tools
module load bowtie2
module load samtools
module load picard/2.23.4

export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"
bowtie2-build $REFERENCE_GENOME $REFERENCE_GENOME
samtools faidx $REFERENCE_GENOME
java -jar $PICARD_ROOT/picard.jar CreateSequenceDictionary R=$REFERENCE_GENOME  O=$REFERENCE_GENOME_DICT


######################################
#### CONCATENTATE LANE DUPLICATES ####
######################################

>concatenate.sh
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J concat
#SBATCH -o concat.%J.out
#SBATCH -e concat.%J.err
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
#forward
for file in /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*L001_R1_001.fastq.gz
do echo `cat ${file} ${file/_L001_/_L002_} > ${file/_L001_R1_001.fastq.gz/}_R1.fq.gz`
done


#CHECK CONCATENATENATION RESULTS MAKE SENSE
ls *R1.fq.gz | wc -l # 106

# Move concatenated files to new folder
mkdir fastqsconcat_2025
mv /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*R1.fq.gz /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/

################################
### parallel gzip and gunzip ###
################################
>unpigz.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J unpigz
#SBATCH -o unpigz.%J.out
#SBATCH -e unpigz.%J.err
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

unpigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*gz'>>unpigz.sh

>pigz.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J pigz
#SBATCH -o pigz.%J.out
#SBATCH -e pigz.%J.err
#SBATCH -t 03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

pigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*fq
# wait
# pigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*tr0
# wait
#pigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*trim
# wait
#pigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*sam
# wait
#pigz -p 8 -v -f /cfs/klemming/projects/supr/naiss2023-23-410/fastqs_2025/*bam
'>>pigz.sh

######################################
############## RUN FASTQ #############
######################################

#RUN FASTQC
module load bioinfo-tools
module load FastQC
mkdir Fastqc_Results/
> runFQC
echo '#!/bin/bash -l' > runFQC
for file in /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/*.fq
do echo "fastqc -o /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/Fastqc_Results -f fastq $file " >> runFQC
done

sh runFQC

ls Fastqc_Results/*zip | wc -l
ls Fastqc_Results/*html | wc -l

# copy fastQC results to local computer

############################
###     PREP READS       ###
############################

# Trimming low-quality bases and remove pcr duplicates. SBATCH command below creates trimming.sh which contains commands for each fastq file. Run this file using: bash trimming.sh
# make sure to use newest trimming scripts, for me that was 2bRAD_trim_launch_dedup_PDW.pl, in which I needed to go in and change to correct filepath to trim2bRAD_dedup_new.pl

# This method is sequential, far from optimal. When submitting job with sbatch I used -p main since shared crashed due to "out of memory".
# takes roughly 1min and 30 sek per command
# ==========================
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p main
#SBATCH -J trimming
#SBATCH -o trimming.%J.out
#SBATCH -e trimming.%J.err
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

ml bioinfo-tools
ml perl' > trimming.sh

2bRAD_trim_launch_dedup_PDW.pl fq >> trimming.sh

# do we have expected number of *.tr0 files created?
ls *R1.fq.tr0 | wc -l # 106

# ===============================
> trimse.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p main
#SBATCH -J trims
#SBATCH -o trims.%J.out
#SBATCH -e trims.%J.err
#SBATCH -t 05:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
ml bioinfo-tools
ml PDC
ml cutadapt' > trimse.sh

# How many reads left?
grep @A00 WF-3727-R-FL-T-7_S61_R1.fq | wc -l
# 3869510
grep @A00 WF-3727-R-FL-T-7.tr0 | wc -l
# 1879806
# So about 49% of reads left
grep @A00 YC-4275-A-FL-A-11-rep2_S104_R1.fq | wc -l
# 9913192
grep @A00 YC-4275-A-FL-A-11-rep2_S104_R1.fq.tr0 | wc -l
# 5580088
# 56% of reads left

# for de novo analysis: removing reads with qualities at ends less than Q15
module load bioinfo-tools
moudle load cutadapt

>trimse
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 36 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse.sh;
done

# for reference-based analysis: trimming poor quality bases off ends:
# removing reads with qualities at ends less than Q15
# mishas script says you can use -m 25 (minimum read length) if you are aligning to genome

for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse.sh;
done

# execute all commands in trimse file (serial or parallel using Launcher, if your system allows)

####### HOR-11 concat ############
#cat HOR-11_aa.trim HOR-11_ab.trim HOR-11_ac.trim HOR-11_ad.trim HOR-11_ae.trim HOR-11_af.trim HOR-11_ag.trim HOR-11_ah.trim > UJ-3099-HOR-11_S63_R1.fq.trim
#gzip HOR-11_a*
##################################

# do we have expected number of *.trim files created?
ls -l *.trim | wc -l # Now together with 2023 samples, 106 (2025) + 63 (2023) = 169 samples in total.

#==============
# Mapping reads to reference (reads-derived fake one, or real) and formatting bam files

# for reference-based:
export REFERENCE_GENOME="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.fasta"
export REFERENCE_GENOME_DICT="/proj/snic2022-23-656/2bRAD_analysis/genome/Zostera_marina.mainGenome.dict"

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
> maps.sh
echo '#!/bin/bash -l
#SBATCH -A snic2022-22-1045
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J maps
#SBATCH -o maps.%J.out
#SBATCH -e maps.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=stefanie.ries@gu.se
module load bioinfo-tools
module load bowtie2' > maps.sh

2bRAD_bowtie2_launch.pl '\.trim$' $REFERENCE_GENOME >> maps.sh
# execute all commands written to mapsls 

ls *.sam | wc -l  # number should match number of trim files
grep 'unpaired' maps | wc -l # just make sure this is also the same as .sam files are written part by part

#============================
####################### extracting unique tags

# 'uniquing'/stacking individual trimmed fastq files:
# step one: create a long list of commands to run in parallel, out of list of files, using 'perl -pe' substitutor:
###########Uniquing script
cd /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025
ls /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/*.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii.sh 
# You can add a multithreading flag (-p 16) to the sge script with:
#sed -i 's/-L 16/-L 16 -p 16/g' unii.sge

#Then open the file and make sure the path is correct, and add a header: # this might run over the weekend, if faster need, seperate to several files
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J unii
#SBATCH -o unii.%J.out
#SBATCH -e unii.%J.err
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl

# submit
sbatch unii.sh

# Done! do you have .uni for all your samples?... 
ls *.uni | wc -l  # 169

# how does it look? Lets display entries from 49,960 to 50,000 :
head -5000 YC-4275-A-FL-A-11-rep2_S104_R1.fq.trim.uni | tail -40


##################### merging uniqued files: 
# again, a three-stage process even though we only have one command to run
###### For tag sharing analysis: #DP: read depth, Ind: individuals, this is very strong filter and would be too much for AFS, #but good for most other purposes, maxInd: set this for 2x the no of ind per pop, #informative tags

#Then open the file and make sure the path is correct, and add a header: 
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p main
#SBATCH -J unimer
#SBATCH -o unimer.%J.out
#SBATCH -e unimer.%J.err
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl

###### collecting nearly all tags (found in at least 2 individuals) for actual genotyping, so here we do another one to have 2 different filters, all tags
echo "mergeUniq.pl uni minDP=5 minInd=2 >all.uniq" > mer2.sh                                                                            

###### merging uniqued files, creating mergedUniqTags.fasta for clustering 
# set minInd to 25% of all your individuals 
# (assuming we want to assemble fake reference genome out of these tags, we only
# need to capture major alleles at each locus)
# set minDP to even higher value (total depth), something like 2 x minInd
echo "mergeUniq.pl uni minDP=84 minInd=42 >all_ref.uniq" > mer.sh                                                           

# how does the merged file look? Lets look at entries from 1800 to 2000:
head -2000 all_ref.uniq | tail -200 | less -S
head -2000 all.uniq | tail -200 | less -S
# Q to exit

# how many different kinds of reads we have 
cat all_ref.uniq | wc -l  # 514527 
cat all.uniq | wc -l  # 10433102

######### discarding reads that have more than 7 observations without reverse-complement
# (do it yourself: calculate random chance that this will happen!), awk programming language for unix
#! these end up with 2 header lines, remove the first manually
head -1 all.uniq >all.tab
awk '!($3>7 && $4==0)' all.uniq >>all.tab  #of this file I removed 2nd header 
    
head -1 all_ref.uniq >all_ref2.tab
awk '!($3>7 && $4==0)' all_ref.uniq >>all_ref2.tab   #of this file I removed 2nd header line manually

######### creating fasta file out of merged and filtered tags:
Recommendation: if you want to alter tables search for awk solution
awk '{print ">"$1"\n"$2}' all.tab | tail -n +3 > all.fasta

awk '{print ">"$1"\n"$2}' all_ref2.tab | tail -n +3 > all_ref2.fasta

# Move on to genotyping

#===============================
#############  actual DE-NOVO GENOTYPING using "fake reference genome" or "reads-derived reference": 

############ Which alleles (tags) belong to the same loci?
# clustering tags allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference, #make sure to have all fasta files here!
# make cd-hit-est.sh script
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH -J ref
#SBATCH -o ref.%J.out
#SBATCH -e ref.%J.err
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl 
module load cd-hit
cd-hit-est -i /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/all.fasta -o /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 10000 -T 0 #was M=0 and T=0    #done

#and same for
cd-hit-est -i /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/all_ref2.fasta -o /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/cdh_all_reftags2.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 10000 -T 0  #done

# do we have new huge cdh_alltags files created?
# list the newest 10 files in the directory:
ls -t | head

############## making fake reference genome (of 34 chromosomes, number of chromosomes for Furcellaria lumbricalis according to Austin 1960) out of cd-hit cluster representatives, ref.sh
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J refgen
#SBATCH -o refgen.%J.out
#SBATCH -e refgen.%J.err
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl 
module load cd-hit

concatFasta.pl fasta=cdh_alltags.fas num=34 
concatFasta.pl fasta=cdh_all_reftags2.fas num=34  #ref2.sh

# from here on, all the procedures are the same as when mapping to a reference genome
## just go on with cdh_alltags.fas, as that is more relevant one

########### formatting / indexing the fake genome
#export GENOME_FASTA=cdh_alltags_cc.fasta
#export GENOME_DICT=cdh_alltags_cc.dict
##########1step
> bowtie.sh
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J bowtie
#SBATCH -o bowtie.%J.out
#SBATCH -e bowtie.%J.err
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl 
module load cd-hit
bowtie2-build cdh_alltags_cc.fasta cdh_alltags_cc.fasta

##########2step
> samtools.sh
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J index
#SBATCH -o index.%J.out
#SBATCH -e index.%J.err
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load samtools
samtools faidx cdh_alltags_cc.fasta

##########3step
> picard.sh
#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J picard
#SBATCH -o picard.%J.out
#SBATCH -e picard.%J.err
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load perl 
module load java
java -jar /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/picard.jar CreateSequenceDictionary R=cdh_alltags_cc.fasta  O=cdh_alltags_cc.dict

################# mapping to the (fake) reference genome and creating BAM files
# for denovo: map reads to fake genome with bowtie2: 
>maps
for F in `ls /cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/*.trim`; do 
done

#add stuff from previous Bowtie script to maps to make it into .sh
sbatch maps #THIS WILL TAKE LONG!! A few hrs   

# checking
ls *.trim.sam > sams
cat sams | wc -l 
# do you have sams for all your samples?... If not, rerun the chunk above

# remove old sam files after checking you have all the files
ll *rg.sam | wc -l
ll *bt2.sam | wc -l
rm *bt2.sam

# Check alignment rates of files (ran as sbatch job). Doesn't seem to work... Manual check majority >80-90% but a few lower. Lowest 49% but is massive outlier.
>alignmentRatesTest
for F in `ls *trim`; do
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRatesTest;
done

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
# adding read groups, validating bams
#-------------------------------

# Add read groups to sam files before converting to bamnano. I can't find this addReadGroup.pl so I did the same procedure but with AddOrReplaceReadGroups from picard (see commands below)
echo '#!/bin/bash -l '> addrg.sh # 1 min for 10 files
for i in *.sam; do
	newfile="$(basename $i .fq.trim.sam)"
	echo "addReadGroup.pl -i $i -o ${newfile}.rg.sam -r ${newfile} -s ${newfile} -l ${newfile} -p Illumina -c NGI" >>addrg.sh;
done
 # Using picard instead. Ran as a job using sbatch, took ~1 hour.
 > addrg.sh
for i in *.sam; do
	newfile="$(basename $i .fq.trim.sam)"
	echo "java -jar picard.jar AddOrReplaceReadGroups -I $i -O ${newfile}.rg.sam --RGLB ${newfile} --RGPL Illumina --RGPU NGI --RGSM ${newfile}" >>addrg.sh;
done

ls *.sam | wc -l
ls *rg.sam | wc -l
# validate sam files
# this takes a long time so skip this step unless you run into an error

#>validateSams.sh # 3 min fr 10 samples
#for i in *rg.sam; do
#	echo "java -jar $PICARD_ROOT/picard.jar ValidateSamFile -I $i -MODE SUMMARY 1>>sam_validationsummary.out 2>>sam_validationsummary.err" >>validateSams.sh;
#done
#>sam_validationsummary.out
#>sam_validationsummary.err
# sh validateSams.sh

# convert, sort and index from sams2b to bam file
export REFERENCE_GENOME="/cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/cdh_alltags_cc.fasta"
export REFERENCE_GENOME_DICT="/cfs/klemming/projects/supr/naiss2023-23-410/fastqsconcat_2025/cdh_alltags_cc.dict"

echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J sambam
#SBATCH -o sambam.%J.out
#SBATCH -e sambam.%J.err
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools' > s2b.sh

for file in *rg.sam; do
echo "samtools sort --reference $REFERENCE_GENOME -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b.sh;
done

# Do we have bam files for all samples?
ls *.bam | wc -l # 169

#################################

# Validate bam files
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J valbam
#SBATCH -o valbam.%J.out
#SBATCH -e valbam.%J.err
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load picard' >validateBams.sh # 1 min 30 sec for 10 samples
>bam_validationsummary.out
>bam_validationsummary.err
for i in *.bam; do
	echo "java -jar picard.jar ValidateSamFile -I $i -R $REFERENCE_GENOME -MODE SUMMARY 1>>bam_validationsummary.out 2>>bam_validationsummary.err" >>validateBams.sh;
done

cat bam_validationsummary.out # No errors for any files

# If everything looks good, remove all sam files
# rm *bt2.sam
# rm *rg.sam


###########################################################################################################################################
# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them. #
###########################################################################################################################################

# Next we check the average depth to see if we proceed with GATK (>10x) or ANGSD (<10x)
# if your coverage is >10x, go to GATK section below

module load bioinfo-tools
module load samtools

> average_depth
for file in *.bam; do
echo "working on $file"
D=`samtools depth $file | awk '{sum+=$3} END { print "Average = ",sum/NR}'`
echo "$file $D" >> average_depth;
done

# listing all bam filenames
ls *bam >bams

# Final amount of bams: 
ls *bam | wc -l
# 169
ls *bam.bai | wc -l
# 169
#----------- assessing base qualities and coverage depth

module load bioinfo-tools
module load samtools
module load ANGSD

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option)
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples (169 samples x 10 = 1690)
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1690 -minInd 10"

# T O   D O :
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd

# scp dd to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ, -minIndDepth and -minInd filters in subsequent ANGSD runs

> basemapping_coveragedepth.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2024-22-1048
#SBATCH -p shared
#SBATCH -J basemap
#SBATCH -o basemap.%J.out
#SBATCH -e basemap.%J.err
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se
module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1690 -minInd 10"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd' > basemap.sh

sbatch basemap.sh

# summarizing results (using modified script by Matteo Fumagalli)
module load RStudio
module load R_packages/4.1.1

Rscript ~/bin/plotQC.R prefix=dd

# copy dd files it in my laptop folder

# proportion of sites covered at >5x:
cat quality.txt

# Coverage for each chromosome of each sample. Results in sbatch output file. Samples in order of ls *bam

# Our depth is not >10x -> using ANGSD

# --- ANGSD pipeline ---

# --- Population structure ---
# Generating genotype likelihoods from highly confident (non-sequencing-error) SNPs
# set minInd to 75-80% of your total number of bams
# if you expect very highly differentiated populations with nearly fixed alternative alleles, remove '-hwe_pval 1e-5' form FILTERS
# -doGeno 8 : genotype likelihood format setting for ngsLD; if you want to run PCA, use -doGeno 32 (but I recommend using ibsMat for all ordination work)

# bams_noclones is the list of all BAM files except for clones and the 7 identified bad samples.
# This leaves us with 152 samples.

# Notes on bam files: minQ can be 30 (99.9% accuracy), as base Q scores are good for all samples. Can try to set 20% (~12 samples lost) and 30% (~25 samples lost) as cutoffs for quality (minQ), at 40% we loose ~80 samples!
# minInd we can set to 80% of number of samples (calculate based on if repeats are included or not). if included 80% of 169 is 135.
# Do 4 sets of filters combining versions of minQ and minMapQ at 20 or 30.
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 20 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 30 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1"

> bamfilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -J filtertest
#SBATCH -o filtertest.%J.out
#SBATCH -e filtertest.%J.err
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 135 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 $FILTERS1 $TODO -P 1 -out result.filter1
angsd -b bams -GL 1 $FILTERS2 $TODO -P 1 -out result.filter2
angsd -b bams -GL 1 $FILTERS3 $TODO -P 1 -out result.filter3
angsd -b bams -GL 1 $FILTERS4 $TODO -P 1 -out result.filter4' > bamfilters.sh

sbatch bamfilters.sh

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 230 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -HWE_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"

# This took a very long time. Over 2 hours for 1 filter...

# Try another round with 2 filters (removing -HWE_pval as Furcellaria is asexual) for bams excluding the bad (bams_nobad) and bads and technical reps (bams_noclones).
> bammorefilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J filtertest2
#SBATCH -o filtertest2.%J.out
#SBATCH -e filtertest2.%J.err
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 107 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -HWE_pval 1e-5"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 107 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -HWE_pval 1e-5"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1"

angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out result.filter9 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS2 $TODO -P 1 -out result.filter10 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS3 $TODO -P 1 -out result.filter11 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS4 $TODO -P 1 -out result.filter12 -nThreads 20' > bammorefilters.sh

# One angsd took ~4 hours to run. Did however not include -nThreads to set number of cores in angsd command. Will do that now and see.
# Including #SBATCH -p 10 and -nThreads 10 decreased time to ~30 minutes per ANGSD run.

# Meeting with Marlene 2025-09-05. Suggested using filters of result.filter10 as base but look into excluding -HWE_pval and -minMaf, and relaxing -hetbias_pval and -minInd. Also to maybe add -C 50, -baq 1 and a -setMaxDepth.
# Also to include -doDepth 1 to get a file of snp depths.
# After researching these filters, I will not use -C but will try -baq 1.
# So will try another round of different filters.

> bammorefilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J filtertest2
#SBATCH -o filtertest2.%J.out
#SBATCH -e filtertest2.%J.err
#SBATCH -t 56:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -hetbias_pval 1e-3 -skipTriallelic 1"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -skipTriallelic 1"
FILTERS5="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 107 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1"
FILTERS6="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -baq 1"
FILTERS7="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -setMaxDepth 100"
FILTERS8="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -setMaxDepth 300"
FILTERS9="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.05 -sb_pval 1e-5 -hetbias_pval 1e-2 -skipTriallelic 1 -setMaxDepth 500"
FILTERS10="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 107 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -baq 1"
FILTERS10="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -baq 1"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1"

angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out result.filter14 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS2 $TODO -P 1 -out result.filter15 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS3 $TODO -P 1 -out result.filter16 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS4 $TODO -P 1 -out result.filter17 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS5 $TODO -P 1 -out result.filter18 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS6 $TODO -P 1 -out result.filter19 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS7 $TODO -P 1 -out result.filter20 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS8 $TODO -P 1 -out result.filter21 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS9 $TODO -P 1 -out result.filter22 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS10 $TODO -P 1 -out result.filter23 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS11 $TODO -P 1 -out result.filter24 -nThreads 20' > bammorefilters.sh

# How many SNPs?
zcat result.filter1.mafs.gz | wc -l # 141 SNPs, very low
zcat result.filter5.mafs.gz | wc -l # 1366 SNPs, still low. Hopefullt the less harsh filtering will provide good results.
zcat result.filter6.mafs.gz | wc -l # 2478 SNPs
zcat result.filter7.mafs.gz | wc -l # 1881 SNPs
zcat result.filter8.mafs.gz | wc -l # 3164 SNPs
zcat result.filter9.mafs.gz | wc -l # 1009 SNPs
zcat result.filter10.mafs.gz | wc -l # 1923 SNPs
zcat result.filter11.mafs.gz | wc -l # 3266 SNPs, Huh, we get more sites when applying the HWE filter?
zcat result.filter12.mafs.gz | wc -l # 5691 SNPs

# Even another round of filters to try to estimate -baq 1 and -C 50. Plus try my suggested filters.
> bammorefilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J filtertest4
#SBATCH -o filtertest4.%J.out
#SBATCH -e filtertest4.%J.err
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700 -baq 1 -ref cdh_alltags_cc.fasta"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700 -baq 1 -ref cdh_alltags_cc.fasta -C 50"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"
FILTERS5="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1"

angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out result.filter25 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS2 $TODO -P 1 -out result.filter26 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS3 $TODO -P 1 -out result.filter27 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS4 $TODO -P 1 -out result.filter28 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS5 $TODO -P 1 -out result.filter29 -nThreads 20' > bammorefilters.sh


# Another round to examine quality scores using -doQsDist 1 and -dumpCounts 2.
> bammorefilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J filtertest4
#SBATCH -o filtertest6.%J.out
#SBATCH -e filtertest6.%J.err
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 124 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 107 -snp_pval 1e-5 -minMaf 0.01 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5700"

TODO="-doMajorMinor 1 -doMaf 1 -doHWE 1 -doQsDist 1 -doCounts 1 -dumpCounts 2 -doDepth 1 -dosnpstat 1"

angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out result.filter30 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS2 $TODO -P 1 -out result.filter31 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS3 $TODO -P 1 -out result.filter32 -nThreads 20' > bammorefilters.sh

# Another round to examine metrics with minInd at 70%, 75% and 80%. After these, look at total depth distribution to set -setMaxDepth.
> bammorefilters.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J filtertest4
#SBATCH -o filtertest7.%J.out
#SBATCH -e filtertest7.%J.err
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 107 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 115 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 122 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1"

TODO="-doMajorMinor 1 -doMaf 1 -doHWE 1 -doQsDist 1 -doCounts 1 -dumpCounts 2 -doDepth 1 -dosnpstat 1"

angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out result.filter33 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS2 $TODO -P 1 -out result.filter34 -nThreads 20
angsd -b bams_nobad -GL 1 $FILTERS3 $TODO -P 1 -out result.filter35 -nThreads 20' > bammorefilters.sh

# After meeting 23/9-2025 we decided on the following final filters:
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 122 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 25000"

# ALWAYS record your filter settings and explore different combinations to confirm that results are robust.
# Our filters will be:
# -minMapQ 30 only highly unique mappings
# -minQ 30 only highly confident base calls
# -minInd 127 the site must be genotyped in at least 80% of individuals. This still gave us plenty of SNPs while increasing depth.
# -snp_pval 1e-5 high confidence that the SNP is not just sequencing error 
# -setMaxDepth 25000 removes exteme outliers in terms of SNP total depth, this value equals mean depth + 4*SD. Removes extreme outliers in terms of depth ~1% of SNPs.
# We thus excluded the following filters: minMaf (barley affected results, and since Furcellaria is largely clonal, we want to capture presence of rare somatic mutations),
# hetbias_pval (had a very large impact on results but drastically reduced SNP and sample depths, plus we are not expecting HWE due to high clonality),
# baq (reduced SNP and sample depths and is not an important filter for 2bRAD data, plus SNP quality is covered by having minQ and minMapQ at 30).

# Starting angsd with -P the number of parallel processes.
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 122 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 25000"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1"
angsd -b bams_noclones -GL 1 $FILTERS1 $TODO -P 1 -out myresult

# how many SNPs?
NSITES=`zcat myresult.mafs.gz | wc -l`
echo $NSITES
# 2324 SNPs

# Removed 2 samples, which were both acting strange and which had been next to each other on the same strip during library preps. Suspect they had been cross-contaminated. R_FL_T_14 and R_FL_A_8.
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 122 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 25000"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1"
angsd -b bams_nobad -GL 1 $FILTERS1 $TODO -P 1 -out myresult2

# how many SNPs?
NSITES=`zcat myresult2.mafs.gz | wc -l`
echo $NSITES
# 18989 SNPs

# Output is myresult2
# Running this took ~1 hour

# Using the ibsMat from ANGSD output we visually decided on a cutoff for clonality based on hierarchical clustering and the distance of technical replicates. 
# We decided to move forward with 2 cutoffs, 0.141 and 0.149. Re-run ANGSD with the bam-list with just one individual per MLL after cutoffs.
# Bam lists are bams-cutoff0141 and bams-cutoff0149. We adjust minInd to 80% of MLLs.
# Run this once without -setMaxDepth in order to be able to calculate mean depth + 4*SD, and then rerun with this filter on. 
# Initial runs showed that mean depth + 4*SD = 9456.549 for cutoff = 0.141 and = 7497.514 for cutoff = 0.149. Run ANGSD again with -setMaxDepth filter on.
> angsdNoClones.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J angsdNoClones
#SBATCH -o angsdNoClones.%J.out
#SBATCH -e angsdNoClones.%J.err
#SBATCH -t 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load samtools
module load ANGSD

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 45 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 9457"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 34 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 7498"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1 -doSaf 1"

angsd -b bams-cutoff0141 -GL 1 $FILTERS1 $TODO -P 1 -out myresult4_noclones_cutoff0141 -nThreads 20 -anc cdh_alltags_cc.fasta
angsd -b bams-cutoff0149 -GL 1 $FILTERS2 $TODO -P 1 -out myresult4_noclones_cutoff0149 -nThreads 20 -anc cdh_alltags_cc.fasta' > angsdNoClones.sh

# How many SNPs?
zcat myresult4_noclones_cutoff0141.mafs.gz | wc -l # 35173 SNPs
zcat myresult4_noclones_cutoff0149.mafs.gz | wc -l # 31118 SNPs


# ------------ Population structure using ngsAdmix -----------

# Run admixture using ngsAdmix
> NGSadmix.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J admix2
#SBATCH -o admix2.%J.out
#SBATCH -e admix2.%J.err
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load ANGSD

NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 2 -o myresult4_noclones_cutoff0141_ngsadmix_K2 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 3 -o myresult4_noclones_cutoff0141_ngsadmix_K3 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 4 -o myresult4_noclones_cutoff0141_ngsadmix_K4 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 5 -o myresult4_noclones_cutoff0141_ngsadmix_K5 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 6 -o myresult4_noclones_cutoff0141_ngsadmix_K6 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 7 -o myresult4_noclones_cutoff0141_ngsadmix_K7 -P 4
NGSadmix -likes myresult4_noclones_cutoff0141.beagle.gz -K 8 -o myresult4_noclones_cutoff0141_ngsadmix_K8 -P 4

NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 2 -o myresult4_noclones_cutoff0149_ngsadmix_K2 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 3 -o myresult4_noclones_cutoff0149_ngsadmix_K3 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 4 -o myresult4_noclones_cutoff0149_ngsadmix_K4 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 5 -o myresult4_noclones_cutoff0149_ngsadmix_K5 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 6 -o myresult4_noclones_cutoff0149_ngsadmix_K6 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 7 -o myresult4_noclones_cutoff0149_ngsadmix_K7 -P 4
NGSadmix -likes myresult4_noclones_cutoff0149.beagle.gz -K 8 -o myresult4_noclones_cutoff0149_ngsadmix_K8 -P 4' > NGSadmix.sh

# Will also do admixture within each site to see if we detect clearer separation of forms. setMaxDepth for each site is based on
# filter value for all populations together adjusted for number of individuals per site. minInd is set to 80% of N in each site.
# To remove annoying ^M from windows text editor: sed -i -e "s/\r//g" A.bamlist.cutoff0141.txt
# For cutoff=0.141 setMaxDepth is 1103 for A, 2679 for K, 3310 for N and 2364 for R.
# For cutoff=0.149 setMaxDepth is 682 for A, 2897 for K, 2727 for N and 1193 for R.

FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 2 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 349"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 13 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 3354"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 13 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 3156"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 2 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 592"

TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -dosnpstat 1 -doHWE 1"

angsd -b A.bamlist.cutoff0149.txt -GL 1 $FILTERS1 $TODO -P 1 -out A_cutoff0149 -nThreads 20
NGSadmix -likes A_cutoff0149.beagle.gz -K 2 -o A_cutoff0149_ngsadmix_K2 -P 4

angsd -b K.bamlist.cutoff0149.txt -GL 1 $FILTERS2 $TODO -P 1 -out K_cutoff0149 -nThreads 20
NGSadmix -likes K_cutoff0149.beagle.gz -K 2 -o K_cutoff0149_ngsadmix_K2 -P 4

angsd -b N.bamlist.cutoff0149.txt -GL 1 $FILTERS3 $TODO -P 1 -out N_cutoff0149 -nThreads 20
NGSadmix -likes N_cutoff0149.beagle.gz -K 2 -o N_cutoff0149_ngsadmix_K2 -P 4

angsd -b R.bamlist.cutoff0149.txt -GL 1 $FILTERS4 $TODO -P 1 -out R_cutoff0149 -nThreads 20
NGSadmix -likes R_cutoff0149.beagle.gz -K 2 -o R_cutoff0149_ngsadmix_K2 -P 4


# --------- Summary statistics ---------
# Fis using PCAngsd --inbreed-samples
pcangsd --beagle myresult4_noclones_cutoff0141.beagle.gz --admix -o pcangsd_cutoff0141 --inbreed-samples --selection --threads 12
pcangsd --beagle myresult4_noclones_cutoff0149.beagle.gz --admix -o pcangsd_cutoff0149 --inbreed-samples --selection --threads 12

# ----- Nucelotide diversity (both variable and invariable sites) for each population using ANGSD. See https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests for info -----
> calculateThetas.sh
echo "#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J calcThetas
#SBATCH -o calcThetas.%J.out
#SBATCH -e calcThetas.%J.err
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load angsd

# Set reference genome
REF="cdh_alltags_cc.fasta"

# List of population files
POPS=("AA.bamlist.cutoff0141" "AM.bamlist.cutoff0141" "KA.bamlist.cutoff0141" "KM.bamlist.cutoff0141" "NA.bamlist.cutoff0141" "NM.bamlist.cutoff0141" "RA.bamlist.cutoff0141" "RM.bamlist.cutoff0141")

# Number of threads
THREADS=4

# Loop through each population
for POP in "${POPS[@]}"; do
  echo "Processing $POP..."

  # Step 1: SAF calculation
  angsd -b ${POP}.txt \
    -ref $REF \
    -anc $REF \
    -out $POP \
    -doSaf 1 \
    -GL 1 \
    -doCounts 1 \
    -minMapQ 30 \
    -minQ 30 \
    -uniqueOnly 1 \
    -remove_bads 1

  # Step 2: Estimate SFS
  realSFS ${POP}.saf.idx -P $THREADS > ${POP}.sfs

  # Step 3: Calculate theta (π)
  realSFS saf2theta ${POP}.saf.idx -sfs ${POP}.sfs -outname ${POP}

  # Step 4: Summarize theta stats
  thetaStat do_stat ${POP}.thetas.idx > ${POP}_theta_summary.txt

  echo "$POP done."
done" > calculateThetas.sh

# For both cutoffs: Job failed (out of memory). AA, AM, KA and KM were completed. Job failed on step 2 for NA. Will resubmit job with NM, RA and RM. And run step 2-4 as seperate job for only NA.

# ------------- Analysis of outlier loci --------------

# --- Bayescan for outliers in Fst ---
# Converting vcf or bcf (using PGDspider) to Bayescan format: 

# make tab-delimited file called bspops LISTING assignments of individuals (as they are named in the vcf file) to populations, for example:
ind1	pop0
ind2	pop0
ind3	pop1
ind4	pop1

# create a file called vcf2bayescan.spid containing this text:
echo "############
# VCF Parser questions (replace VCF with BCF in the line below if you used ANGSD to make bcf instead of vcf)
PARSER_FORMAT=VCF
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=true
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=true
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=true
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=./bspops
# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# GESTE / BayeScan Writer questions
WRITER_FORMAT=GESTE_BAYE_SCAN
# Specify which data type should be included in the GESTE / BayeScan file  (GESTE / BayeScan can only analyze one data type per file):
GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP
############" >vcf2bayescan.spid

# We want to run analyses for all samples (without clone-correction but without technical replicates), all MLLs, all samples within each site and all MLLs within each site. 
# We have the bcf-files for both analyses with MLLs from previous admixture runs above but need to re-run ANGSD with all samples and all samples within each site.
# ANGSD for all samples without technical replicates
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 106 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 25000"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -doHWE 1 -dosnpstat 1 -doDepth 1 -dumpCounts 1"
angsd -b bams_noclones -GL 1 $FILTERS1 $TODO -P 1 -out myresult_5

# ANGSD for all samples without technical replicates within each site. minInd is set to 80% of N, setMaxDepth is based on value for all samples corrected for number of samples in site.
FILTERS1="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 23 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 5451"
FILTERS2="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 19 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 4511"
FILTERS3="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 29 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 6766"
FILTERS4="-uniqueOnly 1 -remove_bads 1 -minMapQ 30 -minQ 30 -minInd 35 -snp_pval 1e-5 -sb_pval 1e-5 -skipTriallelic 1 -setMaxDepth 8270"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doBcf 1 -doPost 1 -doGlf 2 -dosnpstat 1 -doHWE 1"
angsd -b A.bamlist -GL 1 $FILTERS1 $TODO -P 1 -out myresult5_A -nThreads 20
angsd -b K.bamlist -GL 1 $FILTERS2 $TODO -P 1 -out myresult5_K -nThreads 20
angsd -b N.bamlist -GL 1 $FILTERS3 $TODO -P 1 -out myresult5_N -nThreads 20
angsd -b R.bamlist -GL 1 $FILTERS4 $TODO -P 1 -out myresult5_R -nThreads 20

# Converting vcf to bayescan format. I did this locally by getting PGDSpider for windows, converting my bcf file to vcf (using bcftools, see below) and then doing the vcf to bayescan using PGDSpider.
module load bcftools

# MLL bcf files
bcftools convert -O v -o myresult4_noclones_cutoff0141.vcf myresult4_noclones_cutoff0141.bcf
bcftools convert -O v -o myresult4_noclones_cutoff0149.vcf myresult4_noclones_cutoff0149.bcf
bcftools convert -O v -o A_cutoff0141.vcf A_cutoff0141.bcf
bcftools convert -O v -o A_cutoff0149.vcf A_cutoff0149.bcf
bcftools convert -O v -o K_cutoff0141.vcf K_cutoff0141.bcf
bcftools convert -O v -o K_cutoff0149.vcf K_cutoff0149.bcf
bcftools convert -O v -o N_cutoff0141.vcf N_cutoff0141.bcf
bcftools convert -O v -o N_cutoff0149.vcf N_cutoff0149.bcf
bcftools convert -O v -o R_cutoff0141.vcf R_cutoff0141.bcf
bcftools convert -O v -o R_cutoff0149.vcf R_cutoff0149.bcf

# All sample bcf files
bcftools convert -O v -o myresult5.vcf myresult_5.bcf
bcftools convert -O v -o myresult5_A.vcf myresult5_A.bcf
bcftools convert -O v -o myresult5_K.vcf myresult5_K.bcf
bcftools convert -O v -o myresult5_N.vcf myresult5_N.bcf
bcftools convert -O v -o myresult5_R.vcf myresult5_R.bcf

# Transfer these vcf files to local computer and use PGDSpider to convert to bayescan format.

# launching bayescan (this might take 12-24 hours), make sure you have installed bayescan (see top of script) since this was not available on cluster by default.
> bayescanFl.sh
echo '#!/bin/bash -l
#SBATCH -A naiss2025-22-1128
#SBATCH -p shared
#SBATCH -n 10
#SBATCH -J bayescanFl
#SBATCH -o bayescanFl.%J.out
#SBATCH -e bayescanFl.%J.err
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=johan.severinson@gu.se

module load bioinfo-tools
module load PDC
module load bayescan

bayescan best.bayescan.all -o best.bayescan -threads=20

bayescan best.bayescan.cutoff0141 -o best.bayescan.cutoff0141 -threads=20
bayescan best.bayescan.cutoff0149 -o best.bayescan.cutoff -threads=20

bayescan best_A -threads=20
bayescan best_K -threads=20
bayescan best_N -threads=20
bayescan best_R -threads=20

bayescan best_A_cutoff0141 -threads=20
bayescan best_A_cutoff0149 -threads=20
bayescan best_K_cutoff0141 -threads=20
bayescan best_K_cutoff0149 -threads=20
bayescan best_N_cutoff0141 -threads=20
bayescan best_N_cutoff0149 -threads=20
bayescan best_R_cutoff0141 -threads=20
bayescan best_R_cutoff0149 -threads=20' > bayescanFl.sh

# use bayescan_plots.R to examine results


# ----- Outlier analysis using PCAngsd -----
# PCAngsd for running selection scan to find SNPs driving differentiation between populations along axes. Running in addition to Bayescan since that uses Fst (and we don't assume HWE since clonal organism).
# Will do same set of analyses as we did for Bayescan: all samples (before clone correction but without technical replicates), all MLLs, all samples per site, and all MLLs per site.
# We have all the needed beagle files from previous ANGSD runs. Runs fast so can run using salloc instead of sbatch.

# All samples but no technical replicates
pcangsd --beagle myresult_5.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_all_sel --threads 2

# All MLLs
pcangsd --beagle myresult4_noclones_cutoff0141.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_cutoff0141_sel --threads 2
pcangsd --beagle myresult4_noclones_cutoff0149.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_cutoff0149_sel --threads 2

# All samples per site
pcangsd --beagle myresult5_A.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_A_sel --threads 2
pcangsd --beagle myresult5_K.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_K_sel --threads 2
pcangsd --beagle myresult5_N.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_N_sel --threads 2
pcangsd --beagle myresult5_R.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_R_sel --threads 2

# All MLLs per site
pcangsd --beagle A_cutoff0141.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_A_cutoff0141_sel --threads 2
pcangsd --beagle A_cutoff0149.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_A_cutoff0149_sel --threads 2
pcangsd --beagle K_cutoff0141.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_K_cutoff0141_sel --threads 2
pcangsd --beagle K_cutoff0149.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_K_cutoff0149_sel --threads 2
pcangsd --beagle N_cutoff0141.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_N_cutoff0141_sel --threads 2
pcangsd --beagle N_cutoff0149.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_N_cutoff0149_sel --threads 2
pcangsd --beagle R_cutoff0141.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_R_cutoff0141_sel --threads 2
pcangsd --beagle R_cutoff0149.beagle.gz --selection --pcadapt --snp-weights --sites-save --maf-save --admix --out pcangsd_R_cutoff0149_sel --threads 2


# ----- Matching outlier loci to annotated reference genome -----
# Try to see if sites under diversifying selection match annotated regions on Chondrus crispus reference genome
# Use the following code for each run of bayescan or PCAngsd outlier analysis to get sequence thats +-100bp of outlier loci. Just update input lines.
module load bedtools

# Set input files
GENOME="cdh_alltags_cc.fasta"
FLANK=100
SNP_LIST="diversifying_bayescan_cutoff0141.txt"
BED_FILE="snps_flank_bayescan_cutoff0141.bed"
FASTA_OUT="snp_flanks_bayescan_cutoff0141.fasta"

# Step 1: Create BED file with flanking regions
awk -v flank=$FLANK 'BEGIN{OFS="\t"} {
  start = $2 - flank
  end = $2 + flank
  if (start < 0) start = 0
  print $1, start, end, $1 "_" $2
}' $SNP_LIST > $BED_FILE

# Step 2: Extract sequences using BEDtools
bedtools getfasta -fi $GENOME -bed $BED_FILE -fo $FASTA_OUT -name

# Next I transfer the $FASTA_OUT to my computer, go to NCBI and use BLASTn to see if any of the sequences match to the annotated Chondrus crispus reference genome.
