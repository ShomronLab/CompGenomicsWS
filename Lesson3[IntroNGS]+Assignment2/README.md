# Lesson 3
[Class Slides](Lecture3_2025.pdf)

## Tools to be used today:
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
conda install -c bioconda fastqc
```
Their website has a tutorial video, documentation and example reports.
 
## Getting the data:
Through the next few lessons we will make use of publicly available data deposited at the Sequence Read Archive ( [SRA](https://www.ncbi.nlm.nih.gov/sra/?term=ERS5326207) ). We will download two paired-end sequence samples, one taken from normal tissue and the second from Cervical Adenocarcinoma tumor biopsy, both from the same patient, a 55 year old female. </br>
Let's grab all files simultaneously from the SRA ftp server:
```bash
curl -o AD0699_18N_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/007/ERR4833597/ERR4833597_1.fastq.gz &
curl -o AD0699_18N_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/007/ERR4833597/ERR4833597_2.fastq.gz &
curl -o AD0729_18T_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/001/ERR4833621/ERR4833621_1.fastq.gz &
curl -o AD0729_18T_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/001/ERR4833621/ERR4833621_2.fastq.gz &
```
As you can see, we have R1 & R2 for the normal (18N) and the tumor (18T) samples. </br></br>
**Q1:** The file [Metadata_18.csv](Metadata_18.csv) contains additional meta information about the files we are downloading. What does 'neoplasm' and 'normal tissue adjacent to neoplasm' mean? </br>
**Q2:** The 23rd column (LIBRARY_STRATEGY) is marked as WXS. What is Whole Exome Sequencing (WXS/WES)? What would be the fastq file sizes if we did Whole Genome Sequencing (WGS) instead?</br></br>
Once the downloads have completed start by simply looking at each file and making sure the data format is OK. Make sure the read titles make sense and match between R1 and R2.</br>
**Q3:** What is the difference between the header of the first read in R1 and the header of the first read in R2?</br>
**Q4:** What is the read length? How many reads are there per file?</br></br>
Run FastQC on all the files.</br>
Did you get the number of reads per file correct?</br>
**Q5:** Check the "Per base sequence quality". Do the overall base qualities look ok? Do you notice any strange behavior?</br>
**Q6:** Does the "Per base sequence content" behave as expected?</br>

In the next lesson we will learn all about alignment and mapping. Let us pull the human reference genome in a fasta format from the UCSC Genome Browser ( [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html#human) ):
```bash
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
**Q7:** What is hg38?</br>
**Q8:** How many sequences are present in the human reference genome you just downloaded?</br>
**Q9:** Just by looking at the sequence names can you understand what sequences they are?</br>
**Q10:** How many nucleotides are there in chromosome 1? How many in chromosome 20?