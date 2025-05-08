# Lesson 4
[Class Slides](Lecture4_2025.pdf)

## Tools to be used today
[Samtools](http://www.htslib.org/)</br>
[AWS-CLI](https://aws.amazon.com/cli/)</br>
[BWA](https://bio-bwa.sourceforge.net/)</br>
You can create a new conda environment for this assignment:
```bash
conda create -c conda-forge -c bioconda -n read_alignment samtools awscli bwa -y
conda activate read_alignment
```
Or install these tools on your previous one.

## Assignment
In this assignment you will perform sequencing alignment for the data we downloaded on Lesson3. </br>
Let's start by creating the index files for the `bwa mem` algorthim that we will use for mapping and alignment. You can produce these yourself using the below command on the *hg38.fa* genome we previously downloaded:
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
Let's align! The below command will align the normal DNA sequencing sample given all files are present in the current working directory:
```bash
bwa mem -t 8 -R "@RG\tID:normal\tSM:normal" genome.fa AD0699_18N_R1.fastq.gz AD0699_18N_R2.fastq.gz
```
We are executing the *mem* algorithm from the *bwa* suite using 8 threads (you can change this depending on your machine specs). We indicate the group and sample name as *normal* (this information is needed for the future variant caller). We pass as arguments the fasta file of our genome followed by R1 and R2 fastqs. This command will output an **unsorted** SAM file to stdout. </br>
**Q1:** How can we save the SAM to file?</br>
**Q2:** What command will you use to sort by position the SAM file, and output a BAM format?</br>
**Q3:** Write a command that will perform the alignment, sorting and output of a BAM in a single execution.</br>
Next, we will index the final BAM to allow for random access:
```bash
samtools index -@ 8 AD0699_18N.bam
```
**Q4:** What is the most common mapping quality (MAPQ) for the alignments in the BAM file? How many alignments have that mapping quality? What does that mapping quality reflect in terms of the estimated probability that the mapping location is wrong?</br>
**Q5:** The `samtools mpileup` command reports the depth of coverage and the alleles observed at each position in the genome. Using this command, figure out which column of the output represents the depth of sequencing coverage. Report which two positions in the genome have the highest depth of sequencing coverage in this BAM file.</br></br>
Align, sort and index the somatic sample in a similar manner.</br></br>
