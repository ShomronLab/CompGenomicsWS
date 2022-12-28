# Lesson 9
[Class Slides](Slides9.pdf)

## Tools to be used today
[Samtools](http://www.htslib.org/)</br>
[IGV](https://software.broadinstitute.org/software/igv/)</br>
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)</br>
[QualiMap](http://qualimap.conesalab.org/)</br>
[STAR](https://github.com/alexdobin/STAR)</br>
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)</br>
```bash
conda create -c conda-forge -c bioconda -n rnaseq trimmomatic samtools qualimap star igv fastqc -y
conda activate rnaseq
```

## Assignment
In this exercise you will perform alignment and QA of a RNA-Seq experiment. We’ll start from RNA-seq reads, preprocess them, then map the reads to a reference genome and QA the results.

The file [S288C_RNA-seq_rep1_subsample.fq](S288C_RNA-seq_rep1_subsample.fq) is a small subsample from one of six RNA-seq experiments (three replicates) of the yeast strain [S288C](https://www.yeastgenome.org/strain/s288c), a widely used and known strain of yeast. Note that single-end (SE) sequencing was used. Start by running FastQC on the file to get a general idea of the quality of data.</br>
**Q1:** What is the read length? How many reads?

Use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), a read trimming tool for Illumina NGS data, and filter low-quality reads. Filter reads below average quality of 30. </br>
Rerun FastQC on the clean data to make sure filtration went fine.</br>
**Q2:** What percentage of the data was discarded?</br>

Next, we’ll map the filtered RNA-seq reads to the reference genome using STAR. First, we need to create an index. Create a directory called `S288C.STAR`, then run STAR to create the indexed genome in this directory. The reference genome ([S288C_reference_sequence_R64-2-1_20150113.fasta](S288C_reference_sequence_R64-2-1_20150113.fasta)) and annotation [S288C_reference_annotation_proc.gtf](S288C_reference_annotation_proc.gtf). You’ll need to include the option `--genomeSAindexNbases 10` to avoid warning/error messages.

Before we perform the actual mapping, we need to determine the maximal *intron* size we want to use. While most yeast genes don’t have introns, some have quite long ones, as can be seen in the [S288C_reference_annotation_proc.gff](S288C_reference_annotation_proc.gff) file. We can use the `intron` features in the annotation file to determine the largest intron size.</br>
**Q3:** what is the length of the longest intron?</br>
Run STAR again, this time to map the clean reads to the reference genome (note that some of the parameters, e.g. runMode, might change). Use the --alignIntronMax <X> option, with <X>  being the largest intron size + 20%. Also use --outSAMtype BAM SortedByCoordinate, to tell STAR to produce a sorted BAM output. The run should take a few minutes.

While you’re waiting for the run to complete, open an IGV session and load the reference genome (via the `Genomes` tab) and annotation (via the `File` tab). Right-click the annotation track and choose the `expanded` view.</br>
When STAR is finished, you should see the output file `Aligned.sortedByCoord.out.bam`, don't forget to index the bam file.</br>
Load the bam file to IGV and take some time for exploration. Particularly pay attention to factors like the depth distribution and uniformity between genes and intergenic regions.</br>
Look for splice alignments and see if they make sense.</br>
**Q4:** Look at a few examples of genes containing introns:</br> 
Are there reads mapped to the intronic regions?</br>
Were the splice junctions detected correctly?</br>
Attach a relevant screenshot of the region you examined.</br>
Here are some genes with introns for example</br>
| YDL083c | YEL003W | YKL180W | YNL162W |
|---------|---------|---------|---------|


You can simply paste the names into the search box in IGV and hit Enter to zoom to  a specific gene.

Finally, run [QualiMap](http://qualimap.conesalab.org/) for a more systematic QA of the mapping.</br>
You’ll find the report under a new directory called [`input bam name`]_rnaseq_qc. Go over the report. Make sure you understand the meaning of the different statistics. Notice the results for the `Strand specificity estimation` stat.</br> 
**Q5:** What can you say about the library prep protocol used, in terms of strand-specificity?</br>
Rerun Qualimap using the -p option with a value matching the library prep protocol.</br> 
You can run qualimap rnaseq -h to find the possible values for the -p option.</br>
**Q6:** What is the percentage of reads mapped to an exonic region?
