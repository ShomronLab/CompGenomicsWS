# Lesson 6
[Class Slides](Lecture6_2025.pdf)

## Tools to be used today
[Samtools](http://www.htslib.org/)</br>
[IGV](https://software.broadinstitute.org/software/igv/)</br>
[GATK](https://gatk.broadinstitute.org/hc/en-us)</br>
[GSutil](https://cloud.google.com/storage/docs/gsutil)
```bash
conda create -c conda-forge -c bioconda -n variant_calling samtools igv gatk4 gsutil -y
conda activate variant_calling
```

## Assignment
In this assignment you will perform somatic variant calling for the data we processed in the previous assignment. </br>
Load the BAM file into IGV using the File>Load (**make sure you choose refernce genome as hg38**)</br>
**Q1:** Using IGV, what is the total aligned sequencing depth on chromosome 8, position 109,557,831? How many G alleles and A alleles were observed?</br></br>

## Variant Calling
Before we can execute the variant caller we need to produce a few auxiliary files. </br>
We need the fasta file simply indexed in the Samtools indexing format
```bash
samtools faidx genome.fa
```
**Q2:** What information is kept in the samtools' index file?</br>
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
**Q3:** How many filtered candidates are present in the final VCF?</br>
The file [hg38_cosmic70_cervix.txt](hg38_cosmic70_cervix.txt) is a list of somatic mutations found in Cervical Adenocarcinoma patients.</br>
**Q4:** Are all mutations found in the database are also present in our current case? Give one assumption on why this is so?</br>
**Q5:** Find all mutations that **PASS** all filters in Mutect's final VCF and not present in the database. Attach the list and add an IGV screenshot of one of these candidates