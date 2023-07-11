# Slides notes

This document contains notes about the course slides.

This is a somewhat internal document, so expect a fairly draft-ish and concise style.


## 01 Overview


slides 3 - 7

 * the central dogma of molecular biology is know to be not so simple
 * RNA is not only a messenger but may have an effect
 * most of these elements interacts and regulates one-another 
 * alternative splicing in eukaryots adds a layer of possibilities to all this
 * main takeaway maybe : measuring RNA is a proxy for protein levels, which is a proxy for protein activity , which is a proxy for the physiological state of the cell

slides 8 - 9

 * non-exhaustive list of sequencing possibilities 
 * https://liorpachter.wordpress.com/seq/

slides 10-12

 * RNAseq, the challenges
 * slide 10 : from human gff, includes ncRNAs
 * slide 11 : 
 	* data source: https://gtexportal.org/home/datasets (V8 gene TPMs)
 	* left: 1 sample -> from 1 to $10^5$ TPM 
 	* right: 50 random samples -> 10% of genes contribute 90% of the transcripts
 	* NB: mammalian cell 10-30pg RNA/cell , around 360 000 mRNA molecules (source: https://www.qiagen.com/us/resources/faq?id=06a192c2-e72d-42e8-9b40-3171e1eb4cb8&lang=en )

 * slide 12 : important considerations as well

slide 13

 * Illumina : market leader 50-600bp (generally 50-100), 0.1 (nextseq) to 3 (Hiseq) billion reads
 * Ion torrent : 600bp, 260M reads
 * Pacbio : 10-30kb N50 , 4M CCS reads 
 * nanopore : theory single molecule, practice: variable (N50 >100kb on ultra-long kit, up to 4.2Mb )

slides 14-22 : describe different technologies


Ion Torrent :

	- cell sequentially flooded with A T G C 

PacBio SMRT

	- DNApol at bottom of Zero-Mode-Waveguide 
	- fluorescent dye on dNTPs

Illumina seq :

	- formation of clusters with the same sequence 
	- SBS : labelled nucleotides have reversible terminators, so only 1 base is incorporated at a time.


slide 23: paired end sequencing https://www.france-genomique.org/technological-expertises/whole-genome/sequencage-a-courtes-lectures-par-clusterisation/?lang=en

slide 24: stranded sequencing https://www.ecseq.com/support/ngs/how-do-strand-specific-sequencing-protocols-work

slide 25-26: RIN , RNA purification

slide 27-34: sequencing depth and replicates

 * slide 33-34 : this pattern applies to low- mid- and high- expressors (see their supp doc)

slide 35-42: schematic analysis

 * slide 35 : before the sequencing
 * slide 36 : basic analysis
 * slide 37 : basic analysis with trimmed reads
 * slide 38 : main QC steps
 * slide 39 : QC steps provide feedback on previous steps
 * slide 40 : analysis for variant calling / isoform descriptions / ...
 * slide 41 : when no reference genome: de novo assembly
 * slide 42 : the analysis which we'll do during this course
