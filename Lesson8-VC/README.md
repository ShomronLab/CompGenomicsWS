# Lesson 8
[Class Slides](Slides8.pdf)

## Tools to be used today
[Samtools](http://www.htslib.org/)</br>
[BWA](https://bio-bwa.sourceforge.net/)</br>
[GATK](https://gatk.broadinstitute.org/hc/en-us)</br>
[IGV](https://software.broadinstitute.org/software/igv/)</br>
[AWS-CLI](https://aws.amazon.com/cli/)</br>
[GSutil](https://cloud.google.com/storage/docs/gsutil)
```bash
conda create -c conda-forge -c bioconda -n variant_calling samtools bwa gatk4 awscli gsutil igv -y
conda activate variant_calling
```

## Assignment
In this assignment you will perform sequencing alignment and somatic variant calling for the data we downloaded on Lesson5. </br>
Let's start by creating the index files for the `bwa mem` algorthim that we will use for mapping and alignment. You can produce these yourself using the bellow command on the *hg38.fa* genome we previously downloaded
```bash
bwa index hg38.fa
```
If you don't want to wait for the above command to complete, you can download the index files from a publicly open AWS s3 bucket:
```bash
# This is the hg38.fa (no need to re-download if you already have it, just rename)
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa .
# These are the index files
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.amb .
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.ann .
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.bwt .
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.pac .
aws s3 --no-sign-request --region eu-west-1 cp s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa.sa .
```
## Alignment
Let's align! The below command will align the normal DNA sequencing sample given all files are present in the current working directory
```bash
bwa mem -t 8 -R "@RG\tID:normal\tSM:normal" genome.fa AD0699_18N_R1.fastq.gz AD0699_18N_R2.fastq.gz
```
We are executing the *mem* algorithm form the *bwa* suite using 8 threads (you can change this depending on your machine spec). We indicate the group and sample name as *normal* (this information is needed for the variant caller). We pass as arguments the fasta file of our genome follows by R1 and R2 fastqs. This command will output an **unsorted** SAM file to stdout. </br>
**Q1:** How can we save the SAM to file?</br>
**Q2:** What command will you use to sort by position the SAM file and output a BAM format?</br>
**Q3:** Write a command that will perform the alignment, sorting and output a BAM in a single execution.</br>
Next, we will index the final BAM to allow for random access
```bash
samtools index -@ 8 AD0699_18N.bam
```
**Q4:** What is the most common mapping quality (MAPQ) for the alignments in the BAM file? How many alignments have that mapping quality? What does that mapping quality reflect in terms of the estimated probability that the mapping location is wrong?</br>
**Q5:** The `samtools mpileup` command reports the depth of coverage and the alleles observed at each position in the genome. Using this command, figure out which column of the output represents the depth of sequencing coverage report which two positions in the genome have the highest depth of sequencing coverage in this BAM file.</br></br>
Load the BAM file into IGV using the File>Load (**make sure you choose refernce genome as hg38**)</br>
**Q6:** Using IGV, what is the total aligned sequencing depth on chromosome 8, position 109,557,831? How many G alleles and A alleles were observed?</br></br>
Align, sort and index the somatic sample in a similar manner.</br></br>
## Variant Calling
Before we can execute the variant caller we need to produce a few auxiliary files. </br>
We need the fasta file simply indexed in the Samtools indexing format
```bash
samtools faidx genome.fa
```
**Q7:** What information is kept in the samtools' index file?</br>
We need a *dict* file (internal format used almost exclusively by the GATK suite)
```bash
gatk CreateSequenceDictionary -R genome.fa -O genome.dict
```
Finally, we need to download a dataset containing high frequency SNPs and an aggregated normalization reference from multiple samples in order to better the performance of the variant caller.
```bash
gsutil -m cp \
  "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz" \
  "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" \
  "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz" \
  "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" \
  .
```
We will use [`GATK's Mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) to compare the normal to somatic genomes and to identify high likely small somatic variants.
```bash
gatk Mutect2 \
     -R genome.fa \
     -I AD0699_18N.bam \
     -I AD0729_18T.bam \
     -tumor tumor \
     -normal normal \
     -pon 1000g_pon.hg38.vcf.gz \
     --germline-resource af-only-gnomad.hg38.vcf.gz \
     -O somatic.vcf.gz \
     -bamout AD0699_Ad0729.bam
```
Once complete, we will end with default filtering of the somatic candidates
```bash
gatk FilterMutectCalls -R genome.fa -V somatic.vcf.gz -O somatic_filtered.vcf.gz
```
**Q8:** How many filtered candidates are present in the final VCF?</br>
The file [hg38_cosmic70_cervix.txt](hg38_cosmic70_cervix.txt) is a list of somatic mutations found in Cervical Adenocarcinoma patients.</br>
**Q9:** Are all mutations found in the database are also present in our current case? Give one assumption on why this is so?</br>
**Q10:** Find all mutations that **PASS** all filters in Mutect's final VCF and not present in the database. Attach the list and add an IGV screenshot of one of these candidates