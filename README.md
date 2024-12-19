## A Genome Annotation Pipeline Combining MAKER and Liftoff

### Maker Pipeline
```
## Step 1. put all the RNA-seq data in one folder
ls * > list
sed -i "s/1.fq.gz//g" list
sed -i "s/2.fq.gz//g" list  
sort -n list | uniq > list2  ## recover the list of RNA-seq
awk '{print $0 "_hap1.bam"}' list2 > list_bam1 
awk '{print $0 "_hap2.bam"}' list2 > list_bam2  ## create bam file list
```
