# Final Assignment
The final assignment for our DNA sequencing and data analysis course will involve a series of bioinformatics steps and biological reasoning from the long-read sequencing data that was generated in class, with a specific focus on the microbiome. The assignment will include the following tasks:

1. Preparing the long-read sequencing data through trimming and filtering of low-quality reads.

2. Determining the sample of origin using comparative alignment assessment out of a pool of genomes.

3. Locate a functional region in the genome that was sequenced fully or partially and give a biological explanation of its purpose.

4. Conducting microbiome exploration to identify the composition, diversity and richness of the microbial population in the sample, this step will include taxonomy classification and visualization.

5. Choose one highly abundant microbe to discuss in class and assemble its genome. Locate its reference genome and find a metric to assess the correctness of the assembly.

6. Summarize the main findings of the project in a 15-minute presentation to be delivered during the last class meeting. The presentation is to be divided into: technical aspects and methods (5 min), bioinformatics (5 min), and biological findings (5 min).

These tasks will require careful attention to detail and the ability to interpret and communicate the results effectively. The final outcome will be a comprehensive understanding of the bioinformatic and biological concepts related to DNA sequencing, as well as a skill to apply that knowledge on real-world data.


You will be faced with two main challenges:

1. Searching appropriate tools to complete the steps. As oppose to the home assignment, here you won't be served **all** the needed tools for the job and you'll need to find relevant software yourself or write your own code to complete the task. You'll probably find yourself struggling to decide what and how to use software. Don't be afraid to experiment and to try different tools, but reading the docs will save you much time and frustration.
2. Interpretation of the results. Dealing with biological data could be a messy task and the results might be misleading or not clear. Try and avoid convoluted paths and seek simple and clear steps. You should better your data as much as possible, it is preferable (in this assignment) to have less but high quality data. We are here to help if you find yourself stuck in a hurdle you cannot pass.

Our objective (and yours!) is to experiment with real-world sequencing data. There's is **no** right solution here. You are now a bioinformatician dealing with unknown and new dataset.



## Get your data
Links for your raw sequences and run reports.

**You need Reichman username and password to access the links.**

[Team A](https://postidcac-my.sharepoint.com/:f:/g/personal/amit_levon_post_runi_ac_il/EgzEHl0VxtJIrbnlSmzW5GoBe0sDu_T5zn1cRW-d6Lp-pg?e=WwK50C)

[Team B](https://postidcac-my.sharepoint.com/:f:/g/personal/amit_levon_post_runi_ac_il/EhXD92dMOeBDkmLFcbxu2rgB7p2u1Glh_peLgEXyM1i92Q?e=KLSxsG)

[Team C](https://postidcac-my.sharepoint.com/:f:/g/personal/amit_levon_post_runi_ac_il/Eq05Y-nW2AFOuBSsQwOkMiwB35VkkKlQUt8T-Rk4MucS_g?e=uMDuyy)

[Team D](https://postidcac-my.sharepoint.com/:f:/g/personal/amit_levon_post_runi_ac_il/EtrJdqF2hfVCiFr2lTD2VmsB0RRvVjWQ2vgbEALFeAi1SA?e=Sgu9LS)

## QC/A

We will start by examining the FASTQ file produced for your sequencing.

1. Describe general statistics of the sequences and produce visualization for at least 3 metrics.
2. Based on your understanding decide how to trim and/or filter the raw data. You'll need to provide explanation for the logic you implemented during QC.
3. Provide with the same metrics as in 1 and describe the differences.

## **Comparative alignment assessment**

You'll try and decide which animal feces sample you sequenced by aligning the sequences to each one of the below genomes and coming up with a method to evaluate which alignment is "better". Don't worry if you get a low number of alignments, that is to be expected when dealing with fecal samples. The majority of reads are coming from the animal's microbiome.

Present your chosen animal and how and why you decided that.

Here's a list of aligners that you might find useful. The list is not exhaustive, and the popularity of each aligner may vary depending on the specific use case and community

1. [Minimap2](https://github.com/lh3/minimap2) - Minimap2 is a highly popular and versatile sequence aligner that can align long, noisy DNA sequences such as PacBio or Oxford Nanopore reads.
2. [Nanopolish](https://github.com/jts/nanopolish) - Nanopolish is a widely used package for signal-level analysis of Oxford Nanopore data that can be used to improve the quality of consensus sequences and call SNPs and indels with high confidence.
3. [Medaka](https://github.com/nanoporetech/medaka) - Medaka is a widely used neural network-based variant caller that can be used to call SNPs and indels from Oxford Nanopore data.

Genomes:

1. House mouse (GRCm39) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/

2. Red fox (VulVul3) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048418805.1/

3. Jungle cat (FelChav1.0) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019924945.1/

4. Norway rat (GRCr8) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036323735.1/

5. Coyote (Cla-1) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_034620425.1/

6. Little brown bat (Myoluc2.0) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000147115.1/


On the right you have a `Actions` and choose `Download RefSeq` or `Download GenBandk` button. Download the **Genome sequences (FASTA)**

## Functional region that was sequenced

After deciding animal sequences you obtained, you are requested to locate one region with functional capabilities, investigate its role and present it in class. Your animal might have an annotation file you can download, refer to NCBI (the link for your genome) and download a GTF or GFF file. Alternatively, extract the sequences that were aligned to the genome and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) search to check the BLAST database for any functionality.

Functional regions in the genome refer to specific regions of DNA that have a known or predicted function. These regions include:

1. Introns and exons: regions of a gene that are transcribed into RNA. Introns are non-coding regions that are transcribed but removed before the final mRNA is produced, while exons are coding regions that are transcribed and retained in the final mRNA.
2. Regulatory regions: where genes are turned on or off, such as promoter regions and enhancer regions.


## Microbiome exploration

Let's try and determine the origin of the remaining reads that were sequenced and not aligned to our animal of choice. Most of these reads are coming from microbes living alongside our animal. Find a tool that does taxonomy classification (below are two examples). You do not need to build your own database, you can download one compatible with the tool of choice (*read the docs!*). Download a small size database. 

1. [Kraken2](https://ccb.jhu.edu/software/kraken2/) - a software tool that uses a k-mer based classification method to identify the taxonomic origin of DNA sequences. It is designed to work with both short-read and long-read sequencing data, including data from Nanopore sequencing.
2. [Centrifuge](https://ccb.jhu.edu/software/centrifuge/) - a similar tool to Kraken2 that uses a k-mer based classification method and can work with both short-read and long-read sequencing data, including data from Nanopore sequencing.

Present the microbiome diversity. Make a table presenting the top species/genus/etc.. and find a tool to visualize the results. Below are some example tools

1. [Krona](https://github.com/marbl/Krona/wiki) - Interactive metagenomic visualization in a Web browser
2. [MetagenomeScope](https://github.com/marbl/MetagenomeScope) - MetagenomeScope is a web-based platform that supports a wide range of features for metagenomic data analysis, including taxonomic and functional profiling, gene-centric analysis, and interactive visualization.

## Species assembly

Decide on one species that is highly abounded in your sample and assemble its genome. You should pick a species (or a similar one) that has a reference genome deposited in NCBI assembly collection in order to make a comparison with your own results.

Extract all reads associated with your species of choice and assemble them using one of the below tools:

1. [Miniasm](https://github.com/lh3/miniasm) - A lightweight and memory-efficient de novo assembler for long, error-prone reads from Oxford Nanopore Technologies (ONT) platform.
2. [Canu](http://canu.readthedocs.io/en/latest/) - A fork of the Celera Assembler designed for high-noise single-molecule sequencing reads, such as those produced by the Oxford Nanopore and PacBio platforms.
3. [Flye](https://github.com/fenderglass/Flye) - An assembler for long, noisy reads using de Bruijn graphs, and support for both PacBio and Oxford Nanopore data.

Assess your resulting FASTA file using [assembly-stats](https://github.com/sanger-pathogens/assembly-stats)

1. How many contigs are there and what is the length of the longest one?

To compare your assembly to the reference genome use the tool `dnadiff` that is part of the [Mummer package](https://github.com/garviz/MUMmer/blob/master/docs/dnadiff.README). Mummer is a fast aligner that can align complete genomes in relatively short time. (You can choose a different tool if you wish)

Present your results!

### An example presentation from last year -> [Example_presentation.pdf](Example_presentation.pdf)

### Good Luck!



