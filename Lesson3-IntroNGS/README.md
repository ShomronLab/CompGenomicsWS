# Lesson 3
[Class Slides](lecture3.pdf)

## Tools to be used today
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
conda install -c bioconda fastqc
```

## Getting the data
Through the next few lessons we will make use of publicly available data deposited at the [SRA](https://www.ncbi.nlm.nih.gov/sra/?term=ERS5326207). We will download two paired end sequence samples, one taken from normal tissue and the second from Cervical Adenocarcinoma tumor biopsy, both from the same patient, a 55 year old female. </br>
Let's grub all files simoutnesly from the SRA ftp server:
```bash
curl -o AD0699_18N_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/007/ERR4833597/ERR4833597_1.fastq.gz &
curl -o AD0699_18N_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/007/ERR4833597/ERR4833597_2.fastq.gz &
curl -o AD0729_18T_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/001/ERR4833621/ERR4833621_1.fastq.gz &
curl -o AD0729_18T_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR483/001/ERR4833621/ERR4833621_2.fastq.gz &
```
As you can see we have R1 & R2 for the normal (18N) and the tumor (18T) </br></br>
**Q1:** The file [Metadata_18.csv](Metadata_18.csv) contain additional meta information about the files we are downloading. What is 'neoplasm' and 'normal tissue adjacent to neoplasm' mean? </br>
**Q2:** The 23rd column (LIBRARY_STRATEGY) is marked as WXS. What is Whole Exome Sequencing? What will be the fastq file sizes if we did WGS instead?</br></br>
Once downloads have completed start by simply looking at each file and making sure the data format is OK. Make sure read titles make sense and match between R1 and R2.</br>
**Q3:** what is the difference between the header of the first read in R1 and the header of the first read in R2?</br>
**Q4:** What is the read length? How many reads are there per file?</br></br>
Run FastQC on all files.</br>
Did you get the number of reads per file correct?</br>
**Q5:** Check sequence quality per position. Does overall base qualities are ok? Do you notice any strange behavior?</br>
**Q6:** Does per base sequence content behaves as expected?</br>
**Q7:** What does deduplicated sequences mean in the sequence duplication levels plot? Can you find a tool that removes duplicated sequences?</br></br>
On the next lesson we will learn all about alignment and mapping. Let's pull the human reference genome in a fasta format from the [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html#human)
```bash
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
**Q8:** What is hg38?</br>
**Q9:** How many sequences are present in the human reference genome you just downloaded?</br>
**Q10:** Just by looking at the sequence names can you understand what sequences they are?</br>
**Q11:** How many nucleotides are there in chromosome 1? How many in chromosome 20?