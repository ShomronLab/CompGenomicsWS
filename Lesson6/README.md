# Lesson 6
[Class Slides](Slides6.pdf)

## Tools to be used today
[MEGAHIT](https://github.com/voutcn/megahit)
```bash
conda install -c bioconda megahit
```

## Assignment
In this assignment we will develop a simplified assembler to assemble the *de novo* genome of a special strain of a galactic virus (not really). The assembler will use the brute force Shortest Common Superstring (SCS) algorithm discussed in class. We will not overcomplicate the problem with sequencing errors or ploidy. Moreover, as opposed to the data you downloaded in Lesson5, the galactic viral genome was sequenced in a *single end* technique.

**A** Write a function that reads all sequence reads from [CosmicVirus.fastq.gz](CosmicVirus.fastq.gz) and stores them in an array of your choosing.</br>
**Q1:** What is the read length?</br>
**Q2:** How many reads are there in the fastq file?</br></br>

**B** Write an additional function that excepts two sequences, a and b, and a minimum overlap length. The function returns length of longest suffix of a matching a prefix of b that is at least the minimum length. return 0 if no match found.</br></br>

**C** Write an implementation of the brute force method of the Shortest Common Superstring. Remember to iterate on all permutations possible with respect to the order of which the sequences are given. Consider a minimum overlap 1 base and concatenate sequences with no overlap. **Run your algo on first five sequence only!**</br>
**Q3:** How many iteration were needed to complete the SCS finding?</br>
**Q4:** What is the SCS of the first five sequences?</br></br>

Run MEGAHIT for comparison on the entire CosmicVirus sequences</br>
**Q5:** What is the length of the virus reported by MEGAHIT?</br>
**Q6:** MEGAHIT reports an *N50* value. What is this metric using for?

