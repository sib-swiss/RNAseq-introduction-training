
Once you are happy with your read sequences in your FASTQ files, you can use a mapper software to align the reads to the genome and thereby find where they originated from.


**At the end of this lesson, you will be able to :**

 * identify the differences between a local aligner and a pseudo aligner.
 * perform genome indexing appropriate to your data.
 * map your RNA-seq data onto a genome.



## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_04_Mapping.pdf){target=_blank : .md-button }

[STAR website](https://github.com/alexdobin/STAR){target=_blank : .md-button }



## Building a reference genome index

Before any mapping can be achieved, you must first *index* the genome want to map to. 

To do this with STAR, you need two files:

 * a *fasta* file containing the sequences of the chromosome (or genome contigs)
 * a *gtf* file containing annotations (ie. where the genes and exons are)

We will be using the Ensembl references, with their accompanying GTF annotations.

!!! note

	While the data are already on the server here, in practice or if you are following this course without a teacher,
	you can grab the reference genome data from the [Ensembl ftp website](https://www.ensembl.org/info/data/ftp/index.html).

	In particular, you will want a mouse [DNA fasta file](http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/) and [gtf file](http://ftp.ensembl.org/pub/current_gtf/mus_musculus/).

  Take note of the genome sequence and annotation versions, you will need this in your paper's methods section!


**Task :** Using STAR, build a genome index for the mouse mitochondrial chromosome.

 * .fasta and .gtf files are in : `/shared/data/DATA/Mouse_MT_genome/`.
 * create the index in the folder `041_d_STAR_mouseMT_reference`
 * the module name for this aligner is `star`.
 * this job should require less than 4Gb and 10min to run. 

!!! info "STAR basic parameter for genome index generation"

	From the [manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf). Refer to it for more details

	 * `--runMode genomeGenerate` : running STAR in index generation mode
	 * `--genomeDir </path/to/genomeDir>` : output folder for the index
	 * `--genomeFastaFiles </path/to/genome/fasta1>` : chromosome sequences fasta file (can be several files)
	 * `--sjdbGTFfile </path/to/annotations.gtf>` : annotation gtf file
	 * `--runThreadN <NumberOfThreads>` : number of threads to run on 
	 * `--sjdbOverhang <ReadLength-1>` : length of the genomic sequence around the annotated junctions to be used in constructing the splice junctions database. Ideally : read length - 1.

	 Additionally, because the genome is so small here (we only use the mitochondrial chromosome after all), you will need the following advanced option:

	 * `--genomeSAindexNbases 5` : must be scaled to `min(14, log2(GenomeLength)/2 - 1)`, so 5 in our case


!!! note
	
	While your indexing job is running, you can read ahead in STAR's manual to prepare the next step : mapping your reads onto the indexed reference genome.


??? success "STAR indexing script"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-build
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=2
	#SBATCH --mem=3G
	#SBATCH -o 041_l_star_index.o

		
	G_FASTA=/shared/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.dna.chromosome.MT.fa
	G_GTF=/shared/data/DATA/Mouse_MT_genome/Mus_musculus.GRCm39.MT.gtf
	
	ml star

	mkdir -p 041_d_STAR_mouseMT_reference
	
	STAR --runMode genomeGenerate \
	     --genomeDir 041_d_STAR_mouseMT_reference \
	     --genomeFastaFiles $G_FASTA \
	     --sjdbGTFfile $G_GTF \
	     --runThreadN 4 \
	     --genomeSAindexNbases 5 \
	     --sjdbOverhang 99 

	```

	It can be found on the cluster at `/shared/data/Solutions/mouseMT/041_s_star_index.sh`

**Extra task :** Determine how you would add an additional feature to your reference, for example for a novel transcript not described by the standard reference.

??? success "Answer"

	You would edit the gtf file to add your additional feature(s), following the [proper format](https://www.ensembl.org/info/website/upload/gff.html).



!!! note "Note"

	In case you've got multiple FASTA files for your genome (eg, 1 per chromosome), you may just list them with the `genomeFastaFiles` option as follow:

	`--genomeFastaFiles /path/to/genome/fasta1.fa /path/to/genome/fasta2.fa /path/to/genome/fasta3.fa ...`


## Mapping reads onto the reference


**Task :** Using STAR, align the raw FASTQ files of the mouseMT dataset against the mouse mitochondrial reference you just created.

 * if were not able to complete the previous task, you can use the index in `/shared/data/Solutions/mouseMT/041_d_STAR_mouseMT_reference` .
 * search the STAR manual for the option to output a BAM file sorted by coordinate.
 * search the STAR manual for the option to output a geneCounts file.
 * put the results in folder `042_d_STAR_map_raw/` .


!!! info "STAR basic parameters for mapping"

	Taken again from the manual:

	- `--genomeDir </path/to/genomeDir>` : folder where you have put the genome index
	- `--readFilesIn </path/to/read1> ` : path to a fastq file. If the reads are paired, then also include the path to the second fastq file
	- `--runThreadN <NumberOfThreads>`: number of threads to run on.
	- `--outFileNamePrefix  <prefix> ` : prefix of the output files, typically something like `output_directory/sampleName` . 



!!! Note

	Take the time to read the parts of the [STAR manual](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) which concern you: a bit of planning ahead can save you a lot of time-consuming/headache-inducing trial-and-error on your script.


!!! Warning

	Mapping reads and generating a sorted BAM from one of the mouseMT FASTQ file will take less than a minute and very little RAM, but on a real dataset it should take from 15 minutes to an hour per sample and require at least 30GB of RAM.



??? success "STAR mapping script"

	We will be using a job array to map each file in different job that will run at the same time.

	First create a file named `sampleNames.txt`, containing the sample names:

	```
	sample_a1
	sample_a2
	sample_a3
	sample_a4
	sample_b1
	sample_b2
	sample_b3
	sample_b4
	```
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/sampleNames.txt`

	Then for our script:

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=00:10:00
	#SBATCH --cpus-per-task=2
	#SBATCH --mem=1G
	#SBATCH -o 042_l_STAR_map_raw.%a.o
	#SBATCH --array 1-8%8

	ml star

	mkdir -p 042_d_STAR_map_raw

	SAMPLE=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt`

	FASTQ_NAME=/shared/data/DATA/mouseMT/${SAMPLE}.fastq

	STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
                  --outSAMtype BAM SortedByCoordinate \
                  --outFileNamePrefix  042_d_STAR_map_raw/${SAMPLE}. \
                  --quantMode GeneCounts \
                  --readFilesIn $FASTQ_NAME

	```
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/042_s_STAR_map_raw.sh`

	and its results can be found at `/shared/data/Solutions/mouseMT/042_d_STAR_map_raw/`


	The options of STAR are :

	 * `--runThreadN 4 ` : 4 threads to go faster.
	 * `--genomeDir 041_STAR_reference` : path of the genome to map to.
     * `--outSAMtype BAM SortedByCoordinate ` : output a coordinate-sorted BAM file.
     * `--outFileNamePrefix 042_STAR_map_raw/${SAMPLE}.` : prefix of output files.
     * `--quantMode GeneCounts` : will create a file with counts of reads per gene.
     * `--readFilesIn $FASTQ_NAME` : input read file.



## QC on the aligned reads

You can call MultiQC on the STAR output folder to gather a report on the individual alignments.


**Task :** use `multiqc` to generate a QC report on the results of your mapping.

 * Evaluate the alignment statistics. Do you consider this to be a good alignment?
 * How many unmapped reads are there? Where might this come from, and how would you determine this?
 * What could you say about library strandedness ? 

??? success "script and answers"

	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=map-multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o 043_l_multiqc_map_raw.o
	
	multiqc -n 043_r_multiqc_mouseMT_mapped_raw.html -f --title mapped_raw 042_d_STAR_map_raw/
	```
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/043_s_multiqc_map_raw.sh`


	[ Download the report ](../assets/html/043_multiqc_mouseMT_mapped_raw.html){target=_blank : .md-button }


## Comparison of mapping the trimmed reads

After having mapped the raw reads, we also map the trimmed reads and then compare the results to decide which one we want to use for the rest of our analysis.

We will spare you the mapping of the trimmed reads, and let you directly download the mapping multiqc report:


[ trimmed reads mapping  report ](../assets/html/045_multiqc_mouseMT_mapped_trimmed.html){target=_blank : .md-button }



??? note "For the curious: scripts for the mapping of trimmed reads"
	
	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln
	#SBATCH --time=00:10:00
	#SBATCH --cpus-per-task=2
	#SBATCH --mem=1G
	#SBATCH -o 044_l_STAR_map_trimmed.%a.o
	#SBATCH --array 1-8%8
	
	ml star
	
	mkdir -p 044_d_STAR_map_trimmed
	
	SAMPLE=`sed -n ${SLURM_ARRAY_TASK_ID}p sampleNames.txt`
	
	FASTQ_NAME=030_d_trim/${SAMPLE}.trimmed.fastq
	
	STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
	                 --outSAMtype BAM SortedByCoordinate \
	                 --outFileNamePrefix  044_d_STAR_map_trimmed/${SAMPLE}_trimmed. \
	                 --quantMode GeneCounts \
	                 --readFilesIn $FASTQ_NAME
	```
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/044_s_STAR_map_trimmed.sh`


	```sh
	#!/usr/bin/bash
	#SBATCH --job-name=map-trim-multiqc
	#SBATCH --time=00:30:00
	#SBATCH --cpus-per-task=1
	#SBATCH --mem=1G
	#SBATCH -o 045_l_multiqc_mouseMT_mapped_trimmed.o
	
	multiqc -n 045_r_multiqc_mouseMT_mapped_trimmed.html -f --title mapped_trimmed 044_d_STAR_map_trimmed/
	```
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/045_s_multiqc_mouseMT_mapped_trimmed.sh`


## QC report of mapping for the Liu2015 and Ruhland2016 dataset

**Liu2015**

Take the time to look at the following reports: 

[ Liu2015 raw reads mapping  report ](../assets/html/042_r_star_aln_raw_QC_Liu2015.html){target=_blank : .md-button }
[ Liu2015 trimmed reads mapping  report ](../assets/html/044_r_star_map_trim_QC_Liu2015.html){target=_blank : .md-button }

Which one would you choose? 


**Ruhland**

[ Ruhland2016 raw reads mapping  report ](../assets/html/034_r_STAR_multiqc_Ruhland2016.html){target=_blank : .md-button }




## ADDITIONAL : pseudo-aligning with salmon

[salmon website](https://salmon.readthedocs.io/en/latest/salmon.html){target=_blank : .md-button }

Salmon can allow you to quantify transcript expression without explicitly aligning the sequenced reads onto the reference genome with its gene and splice junction annotations, but instead to a simplification of the corresponding transcriptome, thus saving computational resources.

We refer you to the tool's documentation in order to see [how the reference index is computed](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).


**Task :** run salmon to quantify the expression of either the Ruhland or Liu dataset. 
 
 * Use the tool documentation to craft your command line.
 * precomputed indices can be found in `/shared/data/DATA/Mouse_salmon_index` and `/shared/data/DATA/Human_salmon_index`.


??? success "script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=salmonRuhland
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o 033_l_salmon_ruhland2016.%a.o
	#SBATCH --array 1-6%1

	ml salmon
	
	dataDIR=/shared/data/DATA/Ruhland2016
	
	sourceFILE=Ruhland2016.fastqFiles.txt
	
	fastqFILE=`sed -n ${SLURM_ARRAY_TASK_ID}p $sourceFILE`

	genomeDIR=/shared/data/DATA/Mouse_salmon_index

	outDIR=033_d_salmon_Ruhland2016_${fastqFILE%.*}
	
	mkdir -p $outDIR
	
	salmon quant -i $genomeDIR -l A \
				-r $dataDIR/$fastqFILE \
				-p 8 --validateMappings --gcBias --seqBias \
				-o $outDIR
				
	```
	it can also be found in the cluster at `/shared/data/Solutions/Ruhland2016/033_s_salmon_Ruhland2016.sh`


## ADDITIONAL Mapping reads from Ruhland2016 on the reference 

**Task :** Using STAR, align the raw FASTQ files of the Ruhland2016 dataset against thed mouse mitochondrial reference you just created

 * Mapping reads and generating a sorted BAM from one of the Ruhland2016 et al. FASTQ files should take about 20 minutes.
 * Use the full indexed genome at `/shared/data/DATA/Mouse_STAR_index/`, rather than the one we just made.
 * **IMPORTANT**: use the following option in your STAR command: `--outTmpDir /tmp/${SLURM_JOB_USER}_${SLURM_JOB_ID}/`. You can use the manual to look up what this option does. The slurm variables ensure a distinct directory is created in `/tmp/` for each user and for each job.


??? success "STAR mapping script of the Ruhland2016 data"

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
	#SBATCH --job-name=star-aln-Ruhland2016
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o 031_l_STAR_aln_Ruhland2016.%a.o
	#SBATCH --array 1-1%1


	ml star
	outDIR=031_d_STAR_aln_Ruhland2016

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
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/031_s_STAR_aln_Ruhland2016.sh`

	The options of STAR are :

	 * **--runThreadN 8 ** : 8 threads to go faster.
	 * **--genomeDir $genomeDIR** : path of the genome to map to.
     * **--outSAMtype BAM SortedByCoordinate ** : output a coordinate-sorted BAM file.
     * **--outReadsUnmapped Fastx** : output the non-mapping reads (in case we want to analyse them).
     * **--outFileNamePrefix $outDIR/$fastqFILE** : prefix of output files.
     * **--quantMode GeneCounts** : will create a file with counts of reads per gene.
     * **--readFilesIn $dataDIR/$fastqFILE ** : input read file.
     * **--readFilesCommand zcat** : command to unzip the input file.




## ADDITIONNAL : STAR 2-Pass

Genome annotations are incomplete, particularly for complex eukaryotes : there are many as-of-yet unannotated splice junctions.

The first pass of STAR can create a splice junction database, containing both known and novel junctions.
This splice junction database can, in turn, be used to guide an improved second round of alignment, using a command like:

```sh
STAR <1st round options> --sjdbFileChrStartEnd sample_SJ.out.tab
```

**Task :** run STAR in this STAR-2pass mode on the same sample as before and evaluate the results.


??? success "script"

	```
	#!/usr/bin/bash
	#SBATCH --job-name=star-aln2-Ruhland2016
	#SBATCH --time=01:00:00
	#SBATCH --cpus-per-task=8
	#SBATCH --mem=30G
	#SBATCH -o 032_l_STAR_2PASS_Ruhland2016.%a.o
	#SBATCH --array 1-1%1

	ml star

	outDIR=032_d_STAR_Ruhland2016
	
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
	it can also be found in the cluster at `/shared/data/Solutions/mouseMT/032_s_STAR_2PASS_Ruhland2016.sh`


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