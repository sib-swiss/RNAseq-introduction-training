


At the end of this lesson, you will be able to :

 * identify the differences between a local aligner and a pseudo aligner
 * perform genome indexing appropriate to your data
 * map your RNAseq data onto the genome


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_04_Mapping.pdf){: .md-button }

[STAR website](https://github.com/alexdobin/STAR){: .md-button }

## Building a reference genome index

Before any mapping can be achieved, you must first *index* the genome want to map to. 

We will be using the Ensembl versions of iGenome references, with their accompanying GTF annotations.

!!! note

	While the data are already on the server here, in practice, and/or if you are following this course without a teacher
	you can grab the reference genome data from [Ensembl ftp website](https://www.ensembl.org/info/data/ftp/index.html).

	In particular, you will want a mouse [DNA fasta file](http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/) and [gtf file](http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/) \[release-104 at the time we are linking this. Checking for more recent release is recommended\].

**Task :** Using STAR, build a genome index for chromosome 19 of *Mus musculus* using the associated GTF

Important notes :

 * .fasta and .gtf files are in : `/shared/data/DATA/Mouse_chr19/`
 * refer to the [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) to determine which options to use.
 * the `--genomeDir` parameter is the output folder
 * this job should require less than 4Gb and 30min to run.

!!! note
	
	While your indexing job is running, you can read ahead in STAR's manual to prepare the next step : mapping your reads onto the reference.


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
	
	$singularity_exec STAR --runMode genomeGenerate \
	     --genomeDir STAR_references \
	     --genomeFastaFiles $G_FASTA \
	     --sjdbGTFfile $G_GTF \
	     --outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID} \
	     --runThreadN 4
	```

**Extra task :** Determine how you would add an additional feature to your reference, for example for a novel transcript not described by the standard reference.

??? done "Answer"

	Edit the gtf file to add your additionnal feature(s), following the [proper format](https://www.ensembl.org/info/website/upload/gff.html)


!!! Warning

	Remember : request a maximum of 30G and 8 cpus.

## Mapping reads onto the reference

**Task :** Using STAR, align ONE of the FASTQ files from the Ruhland2016 study against the mouse genome.

 * Use the full indexed genome at `/shared/data/DATA/Mouse_STAR_index/`
 * **IMPORTANT**: on the server use the following option in your STAR commands: `--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}/`
 * Generate a BAM file sorted by coordinate
 * Generate a geneCounts file
 * Mapping reads and generating a sorted BAM from one of the Ruhland2016 et al. FASTQ files should take about 20 minutes


!!! Note

	Take the time to read the parts of the STAR [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) which concern you : a bit of planning ahead can save you a lot of time-consuming/headache-inducing trial and error on your script.
	
!!! Warning

	Remember : request a maximum of 30G and 8 cpus.

??? done "STAR mapping script"

	The following sets-up an array of tasks to align all samples.

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

	genomeDIR=/shared/data/Mouse_STAR_index

	STAR --runThreadN 8 --genomeDir $genomeDIR \
                  --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
                  --outFileNamePrefix $outDIR/$fastqFILE \
                  --quantMode GeneCounts \
                  --readFilesIn $dataDIR/$fastqFILE --readFilesCommand zcat \
	              --outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
	```

	The options of STAR are :

	 * **--runThreadN 8 ** : 8 threads to go faster
	 * **--genomeDir $genomeDIR** : path of the genome to map to
     * **--outSAMtype BAM SortedByCoordinate ** : output a sorted BAM file
     * **--outReadsUnmapped Fastx** : output the non-mapping reads (in case we want to analyse them)
     * **--outFileNamePrefix $outDIR/$fastqFILE** : prefix of output files
     * **--quantMode GeneCounts** : will create a file with counts of reads per gene
     * **--readFilesIn $dataDIR/$fastqFILE ** : input read file 
     * **--readFilesCommand zcat** : command to unzip the input file
	 * **--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}** : temporary file folder, for STAR temp files


## QC on the aligned reads

You can call MultiQC on the STAR output folder to gather a report on the alignment.

Here this concern a single sample but usually this would cover all your samples.

**Task :** use `multiqc` to generate a QC report on the result of your mapping.

 * Evaluate the alignment statistics. Do you consider this to be a good alignment?
 * How many unmapped reads are there? Where might this come from, and how would you determine this?
 * what could you say about library strandedness ? 

??? done "script and answers"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o multiqc_star_Ruhland2016.o
	#SBATCH -e multiqc_star_Ruhland2016.e
	
	ml multiqc
	
	mkdir -p STAR_MULTIQC_Ruhland2016/
	
	multiqc -o STAR_MULTIQC_Ruhland2016/ STAR_Ruhland2016/
	```

	result : 

	[ Download the report ](../assets/html/multiqc_report.SRR3180535_EtOH1_1.html){: .md-button }

## ADDITIONNAL : STAR 2-Pass

Genome annotations are incomplete, particularly for complex eukaryotes : ther are many missing splice junctions.

The first pass of STAR can create a splice junction database, known and novel.
This splice junction database can, in turn, be used to guide an improved second round of alignment using a command like:

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
	
	$singularity_exec STAR --runThreadN 8 --genomeDir $genomeDIR \
	                  --outSAMtype BAM SortedByCoordinate \
	                  --outFileNamePrefix $outDIR/$fastqFILE.2Pass. \
	                  --outReadsUnmapped Fastx --quantMode GeneCounts \
	                  --sjdbFileChrStartEnd $outDIR/${fastqFILE}SJ.out.tab \
	                  --readFilesIn $dataDIR/$fastqFILE --readFilesCommand zcat \
	                  --outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}
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