
Once you are happy with your read sequences in your FASTQ files, you can use a mapper software to align the reads to the genome and thereby find where they originated from.


**At the end of this lesson, you will be able to :**

 * identify the differences between a local aligner and a pseudo aligner.
 * perform genome indexing appropriate to your data.
 * map your RNA-seq data onto the genome.



## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_04_Mapping.pdf){: .md-button }

[STAR website](https://github.com/alexdobin/STAR){: .md-button }



## Building a reference genome index

Before any mapping can be achieved, you must first *index* the genome want to map to. 

To do this with STAR, you need two files:
 * a *fasta* file containing the sequences of the chromosome (or genome contigs)
 * a *gtf* file containing annotations (ie. where the genes and exons are)

We will be using the Ensembl references, with their accompanying GTF annotations.

!!! note

	While the data are already on the server here, in practice or if you are following this course without a teacher,
	you can grab the reference genome data from the [Ensembl ftp website](https://www.ensembl.org/info/data/ftp/index.html).

	In particular, you will want a mouse [DNA fasta file](http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/) and [gtf file](http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/) \[release-104 at the time we are linking this. Checking for more recent release is recommended, but may slightly alter the results\].


**Task :** Using STAR, build a genome index for chromosome 19 of *Mus musculus* using the associated GTF

Important notes :

 * the module name for this aligner is `star`.
 * .fasta and .gtf files are in : `/shared/data/DATA/Mouse_chr19/`.
 * refer to the [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) to determine which options to use.
 * the `--genomeDir` parameter is the folder where the indexed genome will be output to.
 * this job should require less than 4Gb and 30min to run.

!!! note
	
	While your indexing job is running, you can read ahead in STAR's manual to prepare the next step : mapping your reads onto the indexed reference genome.


??? done "STAR indexing script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-build
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=4G
	#SBATCH -o star-build.o
	#SBATCH -e star-build.e
		
	G_FASTA=/shared/data/DATA/Mouse_chr19/Mus_musculus.GRCm38.dna.chromosome.19.fa
	G_GTF=/shared/data/DATA/Mouse_chr19/Mus_musculus.GRCm38.101.chr19.gtf
	
	ml star

	mkdir -p STAR_references
	
	STAR --runMode genomeGenerate \
	     --genomeDir STAR_references \
	     --genomeFastaFiles $G_FASTA \
	     --sjdbGTFfile $G_GTF \
	     --runThreadN 4 \
	     --genomeSAindexNbases 11 \
	     --sjdbOverhang 49 

	```


**Extra task :** Determine how you would add an additional feature to your reference, for example for a novel transcript not described by the standard reference.

??? done "Answer"

	Edit the gtf file to add your additional feature(s), following the [proper format](https://www.ensembl.org/info/website/upload/gff.html).


<!-- Suggestion for a note: what to do if you've got multiple FASTA files for your genome, typically 1 per chromosome -->


## Mapping reads onto the reference

**Task :** Using STAR, align ONE of the FASTQ files from the Ruhland2016 study against the mouse genome.

 * Use the full indexed genome at `/shared/data/DATA/Mouse_STAR_index/`, rather than the one we just made.
 * **IMPORTANT**: use the following option in your STAR command: `--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}/`. You can use the manual to look up what this option does. The slurm variables ensure a distinct directory is created in `/tmp/` for each user and for each job.
 * Generate a BAM file sorted by coordinate.
 * Generate a geneCounts file.
 * Mapping reads and generating a sorted BAM from one of the Ruhland2016 et al. FASTQ files should take about 20 minutes.


!!! Note

	Take the time to read the parts of the [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) which concern you : a bit of planning ahead can save you a lot of time-consuming/headache-inducing trial and error on your script.


!!! Warning

	Remember : request a maximum of 30G and 8 CPUs for 1 hour.


??? done "STAR mapping script"


	```
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o star-aln.o
	#SBATCH -e star-aln.e


	ml star
	outDIR=STAR_Ruhland2016

	mkdir -p $outDIR

	dataDIR=/shared/data/DATA/Ruhland2016

	genomeDIR=/shared/data/DATA/Mouse_STAR_index

	STAR --runThreadN 8 --genomeDir $genomeDIR \
                  --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
                  --outFileNamePrefix $outDIR/SRR3180535_EtOH1_1 \
                  --quantMode GeneCounts \
                  --readFilesIn $dataDIR/SRR3180535_EtOH1_1.fastq.gz --readFilesCommand zcat \

	```

	The options of STAR are :

	 * **--runThreadN 8 ** : 8 threads to go faster.
	 * **--genomeDir $genomeDIR** : path of the genome to map to.
     * **--outSAMtype BAM SortedByCoordinate ** : output a coordinate-sorted BAM file.
     * **--outReadsUnmapped Fastx** : output the non-mapping reads (in case we want to analyse them).
     * **--outFileNamePrefix $outDIR/$fastqFILE** : prefix of output files.
     * **--quantMode GeneCounts** : will create a file with counts of reads per gene.
     * **--readFilesIn $dataDIR/$fastqFILE ** : input read file.
     * **--readFilesCommand zcat** : command to unzip the input file.
	 * **--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}** : temporary file folder, for STAR temp files.


??? done "advanced : STAR mapping script with array job"

	The following sets up an array of tasks to align all samples.

	Source file : `Ruhland2016.fastqFiles.txt` :

	```
	SRR3180535_EtOH1_1.fastq.gz
	SRR3180536_EtOH2_1.fastq.gz
	SRR3180537_EtOH3_1.fastq.gz
	SRR3180538_TAM1_1.fastq.gz
	SRR3180539_TAM2_1.fastq.gz
	SRR3180540_TAM3_1.fastq.gz
	```

	sbatch script :

	```
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o star-aln.%a.o
	#SBATCH -e star-aln.%a.e
	#SBATCH --array 1-1%1


	ml star
	outDIR=STAR_Ruhland2016

	mkdir -p $outDIR

	dataDIR=/shared/data/DATA/Ruhland2016

	sourceFILE=Ruhland2016.fastqFiles.txt

	fastqFILE=`sed -n ${SLURM_ARRAY_TASK_ID}p $sourceFILE`

	genomeDIR=/shared/data/DATA/Mouse_STAR_index

	STAR --runThreadN 8 --genomeDir $genomeDIR \
                  --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
                  --outFileNamePrefix $outDIR/$fastqFILE \
                  --quantMode GeneCounts \
                  --readFilesIn $dataDIR/$fastqFILE --readFilesCommand zcat \

	```

	The options of STAR are :

	 * **--runThreadN 8 ** : 8 threads to go faster.
	 * **--genomeDir $genomeDIR** : path of the genome to map to.
     * **--outSAMtype BAM SortedByCoordinate ** : output a coordinate-sorted BAM file.
     * **--outReadsUnmapped Fastx** : output the non-mapping reads (in case we want to analyse them).
     * **--outFileNamePrefix $outDIR/$fastqFILE** : prefix of output files.
     * **--quantMode GeneCounts** : will create a file with counts of reads per gene.
     * **--readFilesIn $dataDIR/$fastqFILE ** : input read file.
     * **--readFilesCommand zcat** : command to unzip the input file.
	 * **--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}** : temporary file folder, for STAR temp files.


## QC on the aligned reads

You can call MultiQC on the STAR output folder to gather a report on the individual alignments.

Here we've aligned a single sample, but usually this would cover all your samples.

**Task :** use `multiqc` to generate a QC report on the results of your mapping.

 * Evaluate the alignment statistics. Do you consider this to be a good alignment?
 * How many unmapped reads are there? Where might this come from, and how would you determine this?
 * What could you say about library strandedness ? 

??? done "script and answers"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o multiqc_star_Ruhland2016.o
	#SBATCH -e multiqc_star_Ruhland2016.e
	
	
	mkdir -p STAR_MULTIQC_Ruhland2016/
	
	multiqc -o STAR_MULTIQC_Ruhland2016/ STAR_Ruhland2016/
	
	```

	Result : 

	[ Download the report ](../assets/html/multiqc_report.SRR3180535_EtOH1_1.html){: .md-button }


## ADDITIONNAL : STAR 2-Pass

Genome annotations are incomplete, particularly for complex eukaryotes : there are many as-of-yet unannotated splice junctions.

The first pass of STAR can create a splice junction database, containing both known and novel junctions.
This splice junction database can, in turn, be used to guide an improved second round of alignment, using a command like:

```sh
STAR <1st round options> --sjdbFileChrStartEnd sample_SJ.out.tab
```

**Task :** run STAR in this STAR-2pass mode on the same sample as before and evaluate the results.


??? done "script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln2
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o star-aln-2pass.%a.o
	#SBATCH -e star-aln-2pass.%a.e
	#SBATCH --array 1-1%1

	ml star

	outDIR=STAR_Ruhland2016
	
	mkdir -p $outDIR
	
	dataDIR=/shared/data/DATA/Ruhland2016
	
	sourceFILE=Ruhland2016.fastqFiles.txt
	
	fastqFILE=`sed -n ${SLURM_ARRAY_TASK_ID}p $sourceFILE`

	genomeDIR=/shared/data/DATA/Mouse_STAR_index
	
	STAR --runThreadN 8 --genomeDir $genomeDIR \
	                  --outSAMtype BAM SortedByCoordinate \
	                  --outFileNamePrefix $outDIR/$fastqFILE.2Pass. \
	                  --outReadsUnmapped Fastx --quantMode GeneCounts \
	                  --sjdbFileChrStartEnd $outDIR/${fastqFILE}SJ.out.tab \
	                  --readFilesIn $dataDIR/$fastqFILE --readFilesCommand zcat \

	```


## ADDITIONAL : pseudo-aligning with salmon

[salmon website](https://salmon.readthedocs.io/en/latest/salmon.html){: .md-button }

salmon can allow you to quantify transcript expression without explicitly aligning the sequenced reads onto the reference genome with its gene and splice junction annotations, but to a simplification of the corresponding transcriptome, thus saving computational resources.

We refer you to the tool's documentation in order to see [how the reference index is computed](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).


**Task :** run salmon to quantify the expression of either the Ruhland or Liu dataset. 
 
 * Use the tool documentation to craft your command line.
 * precomputed indices can be found in `/shared/data/Mouse_salmon_index` and `/shared/data/Human_salmon_index`.


??? done "script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=salmonRuhland
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o salmon_ruhland2016.%a.o
	#SBATCH -e salmon_ruhland2016.%a.e
	#SBATCH --array 1-6%1

	ml salmon

	outDIR=salmon_Ruhland2016
	
	mkdir -p $outDIR
	
	dataDIR=/shared/data/DATA/Ruhland2016
	
	sourceFILE=Ruhland2016.fastqFiles.txt
	
	fastqFILE=`sed -n ${SLURM_ARRAY_TASK_ID}p $sourceFILE`

	genomeDIR=/shared/data/DATA/Mouse_salmon_index
	
	salmon quant -i $genomeDIR -l A \
				-r $dataDIR/$fastqFILE \
				-p 8 --validateMappings --gcBias --seqBias \
				-o $outDIR
				
	```



<!--
## ADDITIONNAL : Assessing read coverage for biases

The [RSeQC](http://rseqc.sourceforge.net/) package includes a function for evaluating [“gene body coverage”](http://rseqc.sourceforge.net/#genebody-coverage-py), which
can be used to assess 5’ or 3’ bias, which might happen if your RNA is degraded or otherwise biased

Requirements:

 * Genome annotations in the 12-column BED format : [you can grab one for mouse here](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/), or at `/data/GRCm38/Mus_musculus.GRCm38.89.bed12` on the server `CHECK LINK. gtf2bed instead?`
 * *Indexed* and sorted BAM file, which can be generated from a sorted BAM using the SAMtools package: `samtools index sample1_sorted.bam`

Example command :
```sh
geneBody_coverage.py -r /data/GRCm38/Mus_musculus.GRCm38.89.bed12 \
                     -i sample1_sorted.bam \
                     -f pdf \
                     -o output_prefix
```

**Task** : Evaluate the gene body coverage for the sample you aligned.

 * expected RAM : 
 * expected time : 

??? done "geneBody_coverage script"

	TBD	
-->