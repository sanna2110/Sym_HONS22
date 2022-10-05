

##### Everything in this file is in directory PCA (cd $SCRATCH/PCA)

# making new .sam files with just symbiont3 (Clade C) 
# made a new file called cladeC, wrote a grep command into this file and then ran it as a job
# from this, I have all new .sam files which start in C and contain only clade C reads

>cladeC
for file in *.sam; do
echo "grep -e 'symbiont3' -e '@' $file > C$file" >> cladeC;
done

ls6_launcher_creator.py -j cladeC -n cladeC -t 1:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 32
sbatch cladeC.slurm

# using samtools sort to make .bam files from my .sam files which only contain clade C reads
# made a new file called Cbam, wrote samtools sort command to this file and ran it as a job
# from this, I have .bam files containing only clade C reads

>Cbam
for file in C*.sam; do 
echo "samtools sort -O bam -o ${file/.sam}.bam $file" >> Cbam;
done

ls6_launcher_creator.py -j Cbam -n Cbam -t 1:00:00 -w 24 -a IBN21018 -e se990@uowmail.edu.au
sbatch Cbam.slurm

# using samtools index to make indexed .bam.csi files
# made a new file called Cindex, wrote samtools index command to this file and ran it as a job

>Cindex
for file in C*.bam; do
echo "samtools index -c $file" >> Cindex;
done

ls6_launcher_creator.py -j Cindex -n Cindex -t 1:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch Cindex.slurm

##### ANGSD Quality Analysis

# within directory PCA, there are two directories which contain failed angsd outputs (angsd and angsd1)

# making a file with a list of bam files (list all files beginning with C and ending in .bam, and redirect the ouptut of this into a file
# called Cbams

ls C*.bam > Cbams

# code used trying -r chr1:1-300000 (output of this is stored in directory angsd1)

# FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000"
# TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
 
# echo "angsd -b Cbams -r chr1:1-300000 -GL 1 $FILTERS $TODO -P 12 -out dd" > ang1

# ls6_launcher_creator.py -j ang1 -n ang1 -q normal -t 24:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
# sbatch ang1.slurm

# trying again without -r option (output of this is what I have used for the rest of the analysis, and is stored in PCA)

# setting FILTERS and TODO to analyse quality of reads using angsd

FILTERS="-uniqueOnly 1 -minMapQ 20 -maxDepth 10000"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# making file called ang2 with the instructions for angsd, and running the job

echo "angsd -b Cbams -GL 1 $FILTERS $TODO -P 12 -out dd" > ang2

ls6_launcher_creator.py -j ang2 -n ang2 -q normal -t 24:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch ang2.slurm

# as the next step didn't work with our original dd.counts file (because it was too large), we have made a new smaller version
# firstly, renaming dd.counts.gz as dd.countsALL.gz
# then, reading the file, piping that into head which will look at the first 300000 lines, and then piping this into gzip which will compress the 
# information and put it into a new file called dd.counts.gz
# we want this to be called the same as our original file, because the Rscript looks for this file name

mv dd.counts.gz dd.countsAll.gz
zcat dd.countsAll.gz | head -300000 | gzip >dd.counts.gz

# running plotQC (output of this will be dd.)

Rscript ~/bin/plotQC.R prefix=dd bams=Cbams
 
# copying dd.pdf to my desktop to look at the quality graphs

scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/HONS/dd.pdf .


##### ANGSD to make matrix

# running angsd with 80% of samples (minInd = 39) (output will be OKall and any files related to this are aa.)

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 39 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

echo "angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OKall" >aa

ls6_launcher_creator.py -j aa -n aa -t 2:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 1
sbatch aa.slurm 

# cat aa.e file (from using minInd=80% of samples)
#        -> Total number of sites analyzed: 332028959
#        -> Number of sites retained after filtering: 274
#        [ALL done] cpu-time used =  10333.45 sec
#        [ALL done] walltime used =  1887.00 sec


# running angsd again with minInd as 50% of our samples (output will be OK50 and any files related to this are bb.)

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 25 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

echo "angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OK50" >bb

ls6_launcher_creator.py -j bb -n bb -t 2:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 1
sbatch bb.slurm 

# cat bb.e file
#        -> Total number of sites analyzed: 332028959
#        -> Number of sites retained after filtering: 2862
#        [ALL done] cpu-time used =  9717.18 sec
#        [ALL done] walltime used =  1486.00 sec


# making a matrix with only autumn spawners

# ls C*A*.bam > 
# FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 23 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
# TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

# echo "angsd -b bams.qc -GL 1 $FILTERS $TODO -P 12 -out OK50" >bb

# ls6_launcher_creator.py -j bb -n bb -t 2:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 1
# sbatch bb.slurm 

# making some copy .bam files to test for clonality

# low coverage sites
cp CMA13.bam cloneMA13.bam
cp CSS8.bam cloneSS8.bam

cp CMS5.bam cloneMS5.bam
cp CSA5.bam cloneSA5.bam

cp CMA2.bam cloneMA2.bam
cp CSA3.bam cloneSA3.bam

ls *.bam > clonebams

FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 39 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

echo "angsd -b clonebams -GL 1 $FILTERS $TODO -P 12 -out OKclone" >cloneang

ls6_launcher_creator.py -j cloneang -n cloneang -t 2:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 1
sbatch cloneang.slurm

# mindInd 25% using six clones
FILTERS="-uniqueOnly 1 -minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 12 -snp_pval 1e-5 -minMaf 0.05 -maxHetFreq 0.5"
TODO="-doMajorMinor 1 -skipTriallelic 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1 -doGlf 2"

echo "angsd -b clonebams -GL 1 $FILTERS $TODO -P 12 -out OK25clone" >cloneang25

ls6_launcher_creator.py -j cloneang25 -n cloneang25 -t 2:00:00 -a IBN21018 -e se990@uowmail.edu.au -w 1
sbatch cloneang25.slurm

# copying clone files
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/OKclone.ibsMat .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/clonebams .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/OK25clone.ibsMat .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/clone25bams .



# copying ibs matrices, file with bam names and quailty.txt to my desktop

scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/OKall.ibsMat .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/OK50.ibsMat .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/bams.qc .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/quality.txt .

scp sa705@ls6.tacc.utexas.edu:/scratch/08912/sa705/HONS/OK.ibsMat .


# zooxType.pl

# making new files with just chr11 and symbiont reads

for F in *.bam; do
echo $F;
samtools view $F chr11 symbiont1 symbiont2 symbiont3 symbiont4 -o ${F/.bam/}.c11s.bam;
done

scp zooxType_symbiont1234.pl sanna21@ls6.tacc.utexas.edu:/home1/08908/sanna21/bin

chmod +x zooxType_symbiont1234.pl

# converting bam files into sam files because I think zooxtype needs sam

for F in *.bam; do
echo $F;
samtools view -h $F -o ${F/.bam/}.sam;
done


zooxType_symbiont1234.pl >> zoox

scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/HONS/zooxtype/zoox .


### FST

export HONS_GENOME_FASTA=$STOCKYARD/honsgen/Amil.symb.cat.genome.fasta

## LOCATION
# calculating saf files for each location population - needed some adjustments so the final file was locsaf4

FILTERS="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -minMaf 0.05 -SNP_pval 1e-5"
TODO="-doMaf 1 -doMajorMinor 1 -doSaf 1 -doCounts 1"

>locsaf4
for pop in *.bams; do
echo "angsd -b $pop -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop/.bams} $FILTERS $TODO -P 8" >>locsaf4
done

ls6_launcher_creator.py -j locsaf4 -n locsaf4 -q normal -t 24:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch locsaf4.slurm

# ashmore
#        -> Total number of sites analyzed: 112126809
#        -> Number of sites retained after filtering: 558

# barrow
#        -> Total number of sites analyzed: 129071281
#        -> Number of sites retained after filtering: 1640

# ningaloo
#        -> Total number of sites analyzed: 135342665
#        -> Number of sites retained after filtering: 989


# calculating FST

realSFS slope.saf.idx lagoon.saf.idx > slope.lagoon.m

# pop1=ashmore
# pop2=barrow
# pop3=ningaloo

realSFS ashmore.saf.idx barrow.saf.idx -P 24 >ashmore.barrow.ml
realSFS ashmore.saf.idx ningaloo.saf.idx -P 24 >ashmore.ningaloo.ml
realSFS barrow.saf.idx ningaloo.saf.idx -P 24 >barrow.ningaloo.ml

realSFS fst index ashmore.saf.idx barrow.saf.idx ningaloo.saf.idx -sfs ashmore.barrow.ml -sfs ashmore.ningaloo.ml -sfs barrow.ningaloo.ml -fstout Fst

realSFS fst stats Fst.fst.idx > fstoutput

# FST values for original version
#         -> FST.Unweight[nObs:324]:0.002504 Fst.Weight:0.002548
# 0.002504        0.002548
#         -> FST.Unweight[nObs:324]:-0.002495 Fst.Weight:-0.003213
# -0.002495       -0.003213
#         -> FST.Unweight[nObs:324]:-0.002606 Fst.Weight:-0.002919
# -0.002606       -0.002919
# pbs.pop1        0.001129
# pbs.pop2        0.001422
# pbs.pop3        -0.004337

# redoing it with new MinInd values set at 70% of each sample group
>locsaf5
for pop in ashmore; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 7 $TODO -P 8" >>locsaf5
done

for pop in barrow; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 17 $TODO -P 8" >>locsaf5
done

for pop in ningaloo; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 10 $TODO -P 8" >>locsaf5
done

ls6_launcher_creator.py -j locsaf5 -n locsaf5 -t 1:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch locsaf5.slurm


# ashmore
#        -> Total number of sites analyzed: 112126809
#        -> Number of sites retained after filtering: 776

# barrow
#        -> Total number of sites analyzed: 129071281
#        -> Number of sites retained after filtering: 370

# ningaloo
#        -> Total number of sites analyzed: 135342665
#        -> Number of sites retained after filtering: 589

realSFS ashmore70.saf.idx barrow70.saf.idx -P 24 >ashmore70.barrow70.ml
realSFS ashmore70.saf.idx ningaloo70.saf.idx -P 24 >ashmore70.ningaloo70.ml
realSFS barrow70.saf.idx ningaloo70.saf.idx -P 24 >barrow70.ningaloo70.ml

realSFS fst index ashmore70.saf.idx barrow70.saf.idx ningaloo70.saf.idx -sfs ashmore70.barrow70.ml -sfs ashmore70.ningaloo70.ml -sfs barrow70.ningaloo70.ml -fstout location70

realSFS fst stats location70.fst.idx > fstoutput

# FST for 70% MinInd location
#         -> FST.Unweight[nObs:226]:0.004741 Fst.Weight:0.004874
#         -> FST.Unweight[nObs:226]:-0.001728 Fst.Weight:-0.002761
#         -> FST.Unweight[nObs:226]:-0.002454 Fst.Weight:-0.002802
# pbs.pop1        0.002464
# pbs.pop2        0.002422
# pbs.pop3        -0.005221


## SAMOENSIS SPRING VS SAMOENSIS AUTUMN

>mixsaf
for pop in samspring; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 7 $TODO -P 8" >>mixsaf
done

for pop in samautumn; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 10 $TODO -P 8" >>mixsaf
done

for pop in samoensis; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 17 $TODO -P 8" >>mixsaf
done

for pop in millepora; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 17 $TODO -P 8" >>mixsaf
done

for pop in north; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 7 $TODO -P 8" >>mixsaf
done

for pop in south; do
echo "angsd -b $pop.bams -GL 1 -ref $HONS_GENOME_FASTA -anc $HONS_GENOME_FASTA -out ${pop}70 $FILTERS -minInd 27 $TODO -P 8" >>mixsaf
done


ls6_launcher_creator.py -j mixsaf -n mixsaf -t 1:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch mixsaf.slurm

# SAMOENSIS SPRING VS AUTUMN

realSFS samspring70.saf.idx samautumn70.saf.idx -P 24 >samspring70.samautumn70.ml

realSFS fst index samspring70.saf.idx samautumn70.saf.idx -sfs samspring70.samautumn70.ml -fstout samspawn70

realSFS fst stats samspawn70.fst.idx

#         -> FST.Unweight[nObs:189]:-0.001502 Fst.Weight:-0.002757


# SAM VS MILL

realSFS samoensis70.saf.idx millepora70.saf.idx -P 24 >samoensis70.millepora70.ml

realSFS fst index samoensis70.saf.idx millepora70.saf.idx -sfs samoensis70.millepora70.ml -fstout species70

realSFS fst stats species70.fst.idx

#         -> FST.Unweight[nObs:303]:0.000740 Fst.Weight:0.000827


# NORTH VS SOUTH

realSFS north70.saf.idx south70.saf.idx -P 24 >north70.south70.ml

realSFS fst index north70.saf.idx south70.saf.idx -sfs north70.south70.ml -fstout nthsth70

realSFS fst stats nthsth70.fst.idx

#         -> FST.Unweight[nObs:280]:0.008847 Fst.Weight:0.007180



scp published_div.fasta.gz sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA

gzip -d published_div.fasta.gz

for file in *.bam; do
samtools fasta $file > ${file/.bam}.fa
done

>ITS
for file in C*.fa; do
echo "ITSx -i $file -o ${file/.fa} -t a" >>ITS
done

ls6_launcher_creator.py -j ITS -n ITS -q normal -t 24:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch ITS.slurm

echo "ITSx -i testMA10.fa -o test" >testITS
ls6_launcher_creator.py -j testITS -n testITS -q normal -t 24:00:00 -w 12 -a IBN21018 -e se990@uowmail.edu.au
sbatch testITS.slurm


scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/test.graph .

scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/ITS2out .
ITS2out


>ITS2out
for file in C*ITS2*; do
echo "${file/.ITS2.fasta}" >> ITS2out
cat $file >> ITS2out
done

cat published_div.fasta | grep '>'  | cut -f 1 -d ' ' |  sed 's/>//g'   > list_of_geneIDs.txt
grep 'C' list_of_geneIDs.txt > subsetIDs.txt
seqtk subseq   published_div.fasta  subsetIDs.txt   > clad_subset.fasta

scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/*ITS*.fasta .
scp sanna21@ls6.tacc.utexas.edu:/scratch/08908/sanna21/PCA/clad_subset.fasta .
