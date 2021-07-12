Read counting refers to the quantification of an *"expression level"*, or abundance, from reads mapped onto a reference genome/transcriptome.
This *expression level* can take several forms, such as a count, or a fraction (RPKM/TPM), and concern different entities (exon, transcript, genes) depending on your biological application.

During this lesson, you will learn to:

 * differentiate between different levels of counting and their relevance for different questions
 * perform read counting at the gene level for Differential Gene expression


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_05_ReadCounting.pdf){: .md-button }

[featureCounts website](http://subread.sourceforge.net/featureCounts.html){: .md-button }

## Read counting with featureCounts

The [featureCount website](http://subread.sourceforge.net/featureCounts.html) provides several useful command-line examples to get started.
For more details on the algorithm behavior (with multi/overlapping reads for instance), you can refer to the package's [User's guide](http://subread.sourceforge.net/SubreadUsersGuide.pdf) (go to the read summarization chapter).


**Task** :

 * Decide which parameters are appropriate for counting reads from the Ruhland dataset. Assume you are interested in determining which genes are differentially expressed.
 * Count the reads from one of your BAM files using featureCount
 * How do the count compare to the counts from STAR ?

 * featureCount : xxx RAM / bam
 * featureCount : xxx cpu time / bam


??? done "featureCounts script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=htseq
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=4G
	#SBATCH -o htseq-count.%a.o
	#SBATCH -e htseq-count.%a.e

	
	G_GTF=/home/SHARED/DATA/Mus_musculus.GRCm38.99.gtf
	
	inFOLDER=${HOME}/Ruhland2016/STAR_Ruhland2016
	outFOLDER=${HOME}/Ruhland2016/HTSQCOUNT_Ruhland2016
	
	ml featureCounts

	mkdir -p $outFOLDER

	featureCounts -T 8 -a /home/SHARED/DATA/Mus_musculus.GRCm38.99.gtf -t exon -g gene_id -o featureCounts_Ruhland2016.counts.txt \
										inFOLDER/EtOH1_Aligned.sortedByCoord.out.bam \
										inFOLDER/EtOH2_Aligned.sortedByCoord.out.bam \
										inFOLDER/EtOH3_Aligned.sortedByCoord.out.bam \
										inFOLDER/TAM1_Aligned.sortedByCoord.out.bam \
										inFOLDER/TAM2_Aligned.sortedByCoord.out.bam \
										inFOLDER/TAM3_Aligned.sortedByCoord.out.bam

	```

??? done "comparison with STAR counts"

	...