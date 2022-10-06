# ------------  INSTALLATIONS:

#----- anaconda install

cd $STOCKYARD
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
idev
bash Anaconda3-2021.11-Linux-x86_64.sh
# make sure to set $STOCKYARD/anaconda3 as target directory
exit
# re-login
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#create a new conda environment to house my r packages etc
conda create --name ReadProcessing

#activate the environment
conda activate ReadProcessing

#check that the environment is active,
#it should have a * next to it if it is active
conda env list

####################################################
########  INSTALL PACKAGES ########################
####################################################

conda install angsd
conda install bowtie2
conda install cutadapt
conda install fastqc
conda install samtools
conda install -c bioconda itsx
conda install -c bioconda multiqc
conda install -c conda-forge r-base=4.1.2

# activate conda
conda activate ReadProcessing

# download all read files (come as .fq.gzip)
# each sample has two files (forward and reverse reads)(R1 & R2)

# decompress all of the fq files
# -d tells gzip to decompress the files instead of compressing them
gzip -d  *.gz &

# downloading vimv (a program which lets us rename files using a text editor) to home
cd
git clone https://github.com/thameera/vimv.git

# copying the vimv file to bin and making it executable
cp vimv/vimv bin
chmod +x $HOME/bin/vimv

# setting nano as the default text editor
nano .bashrc

# scroll down to environment variables and add this in
export EDITOR=nano

# then exit and log back in

# move to where read files are, launch vimv, and edit file names
# use the find and replace function (ctrl-\) to change .fastq with .fq, get rid of unnecessary parts of the name and then replace the sample name with something that makes sense

vimv

# QUALITY CONTROL

# assessing quality measures with fastqc
fastqc -o fastqcout --extract *.fq

# summarising results with multiqc
multiqc .

# CUTADAPT - quality trimming

for file in *_R1.fq; do
echo "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -q 15,15 -m 25 -o ${file/.fq/}.trim -p ${file/_R1.fq/}_R2.trim $file ${file/_R1.fq/}_R2.fq > ${file/_R1.fq/}_trimlog.txt" >> filter;
done

# execute all commands in file filter (recommended: use LAUNCHER module if you are on a cluster)


# BOWTIE2 - aligning to genome

# download genome

# setting up shortcut to genome for ease of coding
export HONS_GENOME_FASTA=Amil.symb.cat.genome.fasta

# indexing genome for bowtie2 and samtools
bowtie2-build $HONS_GENOME_FASTA $HONS_GENOME_FASTA
samtools faidx $HONS_GENOME_FASTA

# mapping reads to genome

>mapping
for file in *R1.trim; do 
echo "bowtie2 --no-unal --local -p 8 -x $HONS_GENOME_FASTA -1 $file -2 ${file/R1.trim/}R2.trim -S ${file/_R1.trim/}.sam" >> mapping;
done

# execute all commands in file mapping (recommended: use LAUNCHER module if you are on a cluster)


# samtools sort to make bam files
>samtobam
for file in *.sam; do 
echo "samtools sort -O bam -o ${file/.sam}.bam $file" >> samtobam;
done

# execute all commands in file samtobam (recommended: use LAUNCHER module if you are on a cluster)

# samtools index to make compressed bam files (.bam.bai) for easy access

>bamindex
for file in *.bam; do
echo "samtools index -c $file" >> bamindex;
done

# execute all commands in file bamindex (recommended: use LAUNCHER module if you are on a cluster)

# identifying proportions of symbiont genera
# making new files with just chr11 and symbiont reads

for F in *.bam; do
echo $F;
samtools view $F chr11 symbiont1 symbiont2 symbiont3 symbiont4 -o ${F/.bam/}.c11s.bam;
done

# converting bam files into sam files because zooxtype needs sam files

for F in *.bam; do
echo $F;
samtools view -h $F -o ${F/.bam/}.sam;
done

# running zooxtype

chmod +x zooxType_symbiont1234.pl
zooxType_symbiont1234.pl >> zoox

### ANGSD

# making new .sam files with just symbiont3 (Cladocopium) 

>cladeC
for file in *.sam; do
echo "grep -e 'symbiont3' -e '@' $file > C$file" >> cladeC;
done

# execute all commands in file cladeC (recommended: use LAUNCHER module if you are on a cluster)

# using samtools sort to make .bam files from my .sam files which only contain clade C reads

>Cbam
for file in C*.sam; do 
echo "samtools sort -O bam -o ${file/.sam}.bam $file" >> Cbam;
done

# execute all commands in file Cbam (recommended: use LAUNCHER module if you are on a cluster)


# using samtools index to make indexed .bam.csi files

>Cindex
for file in C*.bam; do
echo "samtools index -c $file" >> Cindex;
done

# execute all commands in file samtobam (recommended: use LAUNCHER module if you are on a cluster)


##### ANGSD Quality Analysis

# making a file with a list of bam files

ls C*.bam > Cbams

# setting FILTERS and TODO to analyse quality of reads using angsd

FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

angsd -b Cbams -GL 1 $FILTERS $TODO -P 12 -out dd

# making the output file smaller so that the next script can handle it

mv dd.counts.gz dd.countsAll.gz
zcat dd.countsAll.gz | head -300000 | gzip >dd.counts.gz

# summarizing results (using cannibalized Matteo Fumagalli's script)

Rscript ~/bin/plotQC.R prefix=dd bams=Cbams

##### ANGSD to make IBS matrix

# running angsd with 80% of samples (minInd = 39)

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 39 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OKall

# running angsd again with 50% of samples (minInd = 25)

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 25 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OK50

# scp OK50.ibsMat, OKall.ibsMat and bams to laptop, use IBS.R to analyze

# making copy bam files to test for clones
# cp samplename.bam > clonesamplename.bam - do this for five samples, and then run angsd again to make a new ibs matrix, then use clone_check.bash to check for clones

### extracting ITS2 sequences

# scp published_div.fasta.gz from SymPortal

gzip -d published_div.fasta.gz

# converting bam files to fasta files for ITSx

for file in *.bam; do
samtools fasta $file > ${file/.bam}.fa
done

#running ITSx
# -t a: only search for alveolate sequences (makes it a lot quicker)

for file in C*.fa; do
echo "ITSx -i $file -o ${file/.fa} -t a" >>ITS
done

# execute all commands in file ITS (recommended: use LAUNCHER module if you are on a cluster)

# scp all ITS files to laptop and use MEGA to do phylogenetic analysis
