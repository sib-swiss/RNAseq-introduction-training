
Read counting refers to the quantification of an *"expression level"*, or abundance, from reads mapped onto a reference genome/transcriptome.
This *expression level* can take several forms, such as a count, or a fraction (RPKM/TPM), and concern different entities (exon, transcript, genes) depending on your biological application.


**During this lesson, you will learn to:**

 * differentiate between different levels of counting and their relevance for different questions.
 * perform read counting at the gene level for Differential Gene expression.


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_05_ReadCounting.pdf){: .md-button }

[featureCounts website](http://subread.sourceforge.net/featureCounts.html){: .md-button }


## Read counting with featureCounts

<!-- Suggestion: Perhaps add a quite note on why use featureCounts when we've already quantified using STAR. -->

The [featureCount website](http://subread.sourceforge.net/featureCounts.html) provides several useful command-line examples to get started.
For more details on the algorithm behavior (with multi/overlapping reads for instance), you can refer to the package's [User's guide](http://subread.sourceforge.net/SubreadUsersGuide.pdf) (go to the read summarization chapter).


**Task :** 

 * Decide which parameters are appropriate for counting reads from the Ruhland dataset. Assume you are interested in determining which genes are differentially expressed.
 * Count the reads from one of your BAM files using featureCount.
 * How do the featureCount-derived counts compare to those from STAR ?

 * You can find bam files at `/shared/data/Solutions/Liu2015/STAR_Liu2015` and `/shared/data/Solutions/Ruhland2016/STAR_Ruhland2016`
 * featureCount requirements : 400M RAM / BAM file
 * featureCount requirements : 2 min CPU time / BAM file


??? done "featureCounts script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=featurecount
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=4G
	#SBATCH -o count.o
	#SBATCH -e count.e

	
	G_GTF=/shared/data/DATA/Mus_musculus.GRCm39.105.gtf
	
	inFOLDER=/shared/data/Solutions/Ruhland2016/STAR_Ruhland2016
	outFOLDER=featureCOUNT_Ruhland2016
	
	ml subread

	mkdir -p $outFOLDER

	featureCounts -T 8 -a $G_GTF -t exon -g gene_id -o featureCounts_Ruhland2016.counts.txt \
										$inFOLDER/SRR3180535_EtOH1_1.fastq.gzAligned.sortedByCoord.out.bam \
										$inFOLDER/SRR3180536_EtOH2_1.fastq.gzAligned.sortedByCoord.out.bam \
										$inFOLDER/SRR3180537_EtOH3_1.fastq.gzAligned.sortedByCoord.out.bam \
										$inFOLDER/SRR3180538_TAM1_1.fastq.gzAligned.sortedByCoord.out.bam \
										$inFOLDER/SRR3180539_TAM2_1.fastq.gzAligned.sortedByCoord.out.bam \
										$inFOLDER/SRR3180540_TAM3_1.fastq.gzAligned.sortedByCoord.out.bam

	```

??? done "comparison with STAR counts"

	You can use this little R script to check they are the same :

	```
	fc = read.table( "featureCounts_SRR3180535_EtOH1_1.counts.txt" , header =T)
	rownames( fc ) = fc$Geneid
	head( fc )

	star = read.table( "SRR3180535_EtOH1_1.fastq.gzReadsPerGene.out.tab")
	rownames( star ) = star$V1
	head( star )


	star_count = star[ rownames( fc ) , 'V2' ]
	fC_count = fc$STAR_Ruhland2016.SRR3180535_EtOH1_1.fastq.gzAligned.sortedByCoord.out.bam
	plot(log10( star_count + 1),
    	log10(fC_count+1) )

	quantile( star_count  - fC_count)

	```
	