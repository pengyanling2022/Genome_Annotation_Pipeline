## A Genome Annotation Pipeline Combining MAKER and Liftoff

### [Maker Pipeline (Holt et al. 2011)](https://github.com/Yandell-Lab/maker)
+ Prepare RNA-seq files in the same folder
```
ls * > list
sed -i "s/_trim_1.fq.gz//g" list
sed -i "s/_trim_2.fq.gz//g" list  
sort -n list | uniq > list2  ##recover rna list
awk '{print $0 "_hap1.bam"}' list2 > list_bam1
awk '{print $0 "_hap2.bam"}' list2 > list_bam2  ##create bam list
```

+ pasa preprocessing
```
for i in $(cat ./rna/list2);do hisat2 --dta -p 64 -x ./genome_hap1.fa  -1 ./rna/${i}_trim_1.fq.gz -2 ./rna/${i}_trim_2.fq.gz | samtools sort -@ 64 > ${i}_hap1.bam; done
samtools merge -@ 64  hap1.bam $(cat ./rna/list_bam1)
stringtie -p 64 -o hap1.gtf hap1.bam
gffread hap1.gtf -g genome_hap1.fa -w transcripts_hap1.fa
```
+ running pasa
```
~/PASApipeline-pasa-v2.5.2/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../genome_hap1.fa --ALIGNERS minimap2 -t ../transcripts_hap1.fa  --CPU 64
~/PASApipeline-pasa-v2.5.2/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../genome_hap1.fa --ALIGNERS minimap2 -t ../transcripts_hap1.fa  --CPU 1
```
+ runing maker
```
source activate maker2
export LIBDIR=~/miniconda3/envs/maker2/share/RepeatMasker/Libraries
mpiexec -n 64 maker -base species_rnd1 -fix_nucleotides maker_opts.round1.ctl maker_bopts.ctl maker_exe.ctl

genome=/public/home/wwr0530/anno/mushu/genome_hap1.fa
sp=Chardonnay
name=CH1
hap=hap1

source activate maker2
cd ~/${sp}/maker/${hap}/species_rnd1.maker.output
gff3_merge -s -d species_rnd1_master_datastore_index.log > species_rnd1.all.maker.gff
fasta_merge -d species_rnd1_master_datastore_index.log
gff3_merge -n -s -d species_rnd1_master_datastore_index.log > species_rnd1.all.maker.noseq.gff

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' species_rnd1.all.maker.noseq.gff | \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2- 1000, $3+1000 }' | \
bedtools getfasta -fi ${genome} -bed - -fo species_rnd1.all.maker.transcripts1000.fasta
awk '{ if ($2 == "est2genome") print $0 }' species_rnd1.all.maker.noseq.gff > species_rnd1.all.maker.est2genome.gff
awk '{ if ($2 == "protein2genome") print $0}' species_rnd1.all.maker.noseq.gff > speices_rnd1.all.maker.protein2genome.gff
awk '{ if ($2 ~ "repeat") print $0}' speices_rnd1.all.maker.noseq.gff > speices_rnd1.all.maker.repeats.gff

cd ~/${sp}/maker/${hap}/species_rnd1.maker.output
mkdir snap
cd snap
maker2zff -l 50 -x 0.5 ../species_rnd1.all.maker.gff
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl ${name} params > ${name}_rnd1.hmm

cd ~/${sp}/maker/${hap}/species_rnd1.maker.output
mkdir augustus
cd augustus
awk '{ if ($2 == "maker") print $0}' ../species_rnd1.all.maker.gff > maker_rnd1.gff
gff2gbSmallDNA.pl maker_rnd1.gff ${genome} 2000 ${name}_rnd1.gb
new_species.pl --species=${name}_rnd1
etraining --species=${name}_rnd1 ${name}_rnd1.gb
randomSplit.pl ${name}_rnd1.gb 200
mv ${name}_rnd1.gb.test ${name}_rnd1.gb.evaluation
augutus --species=${name}_rnd1 ${name}_rnd1.gb.evaluation >& first_evaluate.out
grep -A 22 Evaluation first_evaluate.out
randomSplit.pl ${name}_rnd1.gb 2000
optimize_augustus.pl --species=${name}_rnd1 --kfold=24 --cpus=64 --onlytrain=${name}_rnd1.gb.train
${name}_rnd1.gb.test
etraining --species=${name}_rnd1 ${name}_rnd1.gb
augustus --species=${name}_rnd1 ${name}_rnd1.gb.evalutation >& second_evaluate.out
grep -A 22 Evaluation second_evaluate.out
```
### [Liftoff Pipeline (Shumate et al. 2020)](https://github.com/agshumate/Liftoff)
+running liftoff
```
source activate liftoff
liftoff -g PN_T2T.gff3 CH_hap1.fa PN_T2T.fa -o CH_hap1.gff3
```

### Combing MAKER and Liftoff results
+ extract genes based on liftoff
```
ID1=L_CH_hap1
awk '$3 == "gene" {print $1"\t"$4"\t"$5"\t"$9}' $ID1.gff3 | sed 's/;.*//;s/ID=//' > $ID1.bed
```
+ extract genes based on maker
```
ID2=M_CH_hap1
awk '$3 == "gene" {print $1"\t"$4"\t"$5"\t"$9}' $ID2.gff| sed 's/;.*//;s/ID=//' > $ID2.bed
```
+ combine the genes and remove the redudant information
```
ID3=hap1
bedtools intersect -a $ID1.bed -b $ID2.bed -v >  $ID3.bed
awk '{print $0"\t"($3-$2)}' $ID3.bed > $ID3.bed2
cut -f 4 $ID3.bed2 >addlist1.txt
grep -f addlist1.txt $ID1.gff3 | grep -v "#" > ${ID3}_add.gff3
cat $ID2.gff ${ID3}_add.gff3 > $ID3.final.gff3
awk '$3 == "gene" ||$3 == "mRNA"|| $3 == "exon" || $3 == "five_prime_UTR" || $3 == "three_prime_UTR" || $3 == "CDS"' $ID3.final.gff3| sed 's/Liftoff/\./g'|sed 's/Maker/\./g' > $ID3.final1.gff3
source activate py3
gff3_sort -g $ID3.final1.gff3 -og $ID3.sorted.gff3
python3 filter.py $ID3.sorted.gff3 $ID3.filter.gff3
python3 rename.py -g $ID3.filter.gff3 -c key1.txt -p $ID3
gff3_sort -g $ID3.rename.gff3 -og $ID3.rename.sorted.gff3

rm addlist1.txt rmlist1.txt ${ID3}_rm.gff3 $ID3.final.gff3  ${ID3}_add.gff3 $ID3.bed2 $ID3.bed $ID3.final1.gff3 $ID3.sorted.gff3 $ID3.filter.gff3 $ID3.filter.sorted.gff3 $ID3.rename.gff3
```

