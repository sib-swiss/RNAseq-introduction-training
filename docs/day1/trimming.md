

Following a QC analysis on sequencing results, one could detect stretches of low quality bases along reads, or a contamination by adapter sequence.
Depending on your application and the sofware you use for mapping, you may have to remove these bad quality / spurious sequences out of your data.

During this block, you will learn to :

 * trim your data with trimmomatic


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_03_trimming.pdf){: .md-button }

[Trimmomatic website](http://www.usadellab.org/cms/?page=trimmomatic){: .md-button }


## to trim or not to trim ?

If the data will be used to perform **transcriptome assembly, or variant analysis, then it must be trimmed**.


In contrast, for applications based on **counting reads, such as Differential Expression analysis**, most aligners, such as [STAR](https://github.com/alexdobin/STAR), [HISAT2](http://daehwankimlab.github.io/hisat2/), [salmon](https://salmon.readthedocs.io/en/latest/salmon.html), [kallisto](https://pachterlab.github.io/kallisto/manual) can handle bad quality sequences and adapter content by performing "soft-clipping" on reads, and consequently they **usually do not need trimming.
In fact, hard trimming can be detrimental** to the number of successfully quantified reads \[[William et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2)\].

Nevertheless, they usually recommend to perform some amount of soft trimming (*eg.* [kallisto](https://www.biostars.org/p/389324/), [salmon](https://github.com/COMBINE-lab/salmon/issues/398) ).
If possible, we recommend to perform the mapping for both the raw data and the trimmed one, in order to compare the results for both and choose the best.


## trimming with Trimmomatic


The [trimmomatic website](http://www.usadellab.org/cms/?page=trimmomatic) gives very good examples of their software usage for both paired-end (`PE`) and single-end (`SE`) reads. We recommend you read their quick-start section attentively.


**Task :** Conduct a soft trimming on the Liu2015 data.

**Extra (if you have the time) :** run a QC analysis on your trimmmed reads and compare with the raw ones.

Important notes :

 * Adapter sequences can be found in `/home/SHARED/DATA/adapters/`
 * Trimmomatic RAM requirements : 
 * Trimmomatic time requirements : ~ 10 min/ read file 



??? done "Trimmomatic script"

	The Liu2015 dataset has paired-end reads and we have to take that into account during trimming.
	For a soft-trimming, we chose the following options :

	 * **SLIDINGWINDOW:4:20** Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
	 	- **4**  : windowSize: specifies the number of bases to average across
	 	- **20** : requiredQuality: specifies the average quality required.
	 * **ILLUMINACLIP:/shared/home/SHARED/DATA/adapters/TruSeq3-PE.fa:2:30:10** Cut adapter and other illumina-specific sequences from the read.
	 	 - Cut adapter and other illumina-specific sequences from the read.
	 	 - **2**  : seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
	 	 - **30** : palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
	 	 - **10** : simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.


	Here is a script for a single sample : 

	```sh
	#!/bin/bash
	#SBATCH --job-name=trim
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem-per-cpu=4G
	#SBATCH -o trim.o
	#SBATCH -e trim.e
	
	ml 	trimmomatic
	
	dataDIR=/shared/home/SHARED/DATA/Liu2015
	
	trimmomatic="/easybuild/apps/Trimmomatic/0.36-Java-1.7.0_80/trimmomatic-0.36.jar"
	
	outDIR=Liu2015_trimmed_reads
	
	mkdir -p $outDIR
	
	$trimmomatic PE -threads 4 -phred33 \
	             $dataDIR/SRR1272187_1.fastq.gz \
	             $dataDIR/SRR1272187_2.fastq.gz \
	             $outDIR/SRR1272187_NFLV_trimmed_paired_1.fastq $outDIR/SRR1272187_NFLV_trimmed_unpaired_1.fastq \
	             $outDIR/SRR1272187_NFLV_trimmed_paired_2.fastq $outDIR/SRR1272187_NFLV_trimmed_unpaired_2.fastq \
	             SLIDINGWINDOW:4:20 ILLUMINACLIP:/adapters/TruSeq3-PE.fa:2:30:10 
	
	## compressing the resulting fastq files to save some space.
	gzip $outDIR/SRR1272187_NFLV_trimmed_paired_1.fastq
	gzip $outDIR/SRR1272187_NFLV_trimmed_unpaired_1.fastq
	gzip $outDIR/SRR1272187_NFLV_trimmed_paired_2.fastq
	gzip $outDIR/SRR1272187_NFLV_trimmed_unpaired_2.fastq
	
	```