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


## 02 Quality control


slide 04 : https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm

 * control bit: 0 when none of the control bits are on, otherwise it is an even number. On HiSeq X and NextSeq systems, control specification is not performed and this number is always 0.

slide 16-... : interpretation of fastQC report


## 03 trimming

slide 02-05 : enumerating reasons we may want to trim

 * slide 04 : adapter can be present if, for instance, insert size is shorter than sequenced length. https://www.ecseq.com/support/ngs/trimming-adapter-sequences-is-it-necessary

slide 06-11 : different cases where we may trim

## 04 mapping

slides 4-8:

 * bowtie2 : https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
 * STAR : https://github.com/alexdobin/STAR
 * tophat2 : https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/tophat_tophat2/
 * RSEM : https://github.com/deweylab/RSEM
 * cufflinks : http://cole-trapnell-lab.github.io/cufflinks/
 * salmon : https://salmon.readthedocs.io/en/latest/salmon.html
 * kallisto : https://pachterlab.github.io/kallisto/
 * tximport : https://bioconductor.org/packages/release/bioc/html/tximport.html
 * featurecount : https://subread.sourceforge.net/featureCounts.html
 * stringtie : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 * GATK variant calling pipeline https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-


## 05 DE

slide 3-9 : challenges for RNAseq 

 * slide 4: sequencing depth varies accross libraries
     * left : 100 samples from the Gtex V8 dataset: https://gtexportal.org/home/datasets
     * right : samples from a random binomial with and without a library size factor applied

 * slide 5: most of the expression is taken by very few genes + a lot of genes have 0 reads  (plots: data from 100 samples from the Gtex V8 dataset: https://gtexportal.org/home/datasets )

 * slide 6: small number of samples. 10k simulation of negative binomial draws

 * slide 7-9 : xkcd.com/882

slide 10-11 : input

slide 12-14 : https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

slide 15: filtering: 
    * https://academic.oup.com/bioinformatics/article/29/17/2146/240530 max and mean refer to CPM thresholds ; CPM 1: genes with a CPM less than one in more than half the samples are filtered
    * http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt
    * section 2.7 of https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

slide 16: normalization

https://www.biostars.org/p/284775/
EdgeR: Trimmed Mean of M-values (TMM): 
        
        also based on the hypothesis that most genes are not DE. 

        The TMM factor is computed for each lane, with one lane being considered as a reference sample and the others as test samples. 

        For each test sample, TMM is computed as the weighted mean of log ratios between this test and the reference, after exclusion of the most expressed genes and the genes with the largest log ratios. 

        According to the hypothesis of low DE, this TMM should be close to 1. If it is not, its value provides an estimate of the correction factor that must be applied to the library sizes (and not the raw counts) in order to fulfill the hypothesis. 

        [source: https://www.ncbi.nlm.nih.gov/pubmed/22988256]
DESeq2
    DESeq:  is based on the hypothesis that most genes are not DE. 
            
        the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. 

        The underlying idea is that non-DE genes should have similar read counts across samples, leading to a ratio of 1. 

        Assuming most genes are not DE, the median of this ratio for the lane provides an estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis. 

            [source: https://www.ncbi.nlm.nih.gov/pubmed/22988256]


slide 17 : https://pubmed.ncbi.nlm.nih.gov/22988256/
    * TC: Total count (CPM) - UQ: Upper Quartile - Med: median - Q: quantile
    * top left: coef of variation in housekeeping genes in H. sapiens data
    * top right: average false-positive rate over 10 independent datasets simulated with varying proportions of differentially expressed genes (from 0% to 30% for each normalization method). 
    * bottom:
        * distribution: distribution inter samples look the same
        * Intra-variance: intra group variance 
        * Housekeeping : coef of variation in 30 housekeeping genes, which are supposed similarly expressed
        * clustering: similarity of DE genes with other methods
        * false positive rate : see above

slide 19: NB model
    * https://doi.org/10.1186/gb-2010-11-10-r106
    * orange line is the fit w(q)
    * purple line show the variance implied by the Poisson distribution 
    * dashed orange line is the variance estimate used by edgeR. 


## 06 Enrichment

08: geneontology.org
09: reactome.org
10: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
11: https://www.genome.jp/kegg/

15 to 21 : GSEA 

GSEA is used a lot, but it comes in many flavour, with a large number of options
and it is not always easy to understand what is happening.

"The enrichment score (ES) represents the degree to which a set S is over-represented at the top or bottom of the ranked list L. The score is calculated by walking down the list L, increasing a running-sum statistic when we encounter a gene in S and decreasing when it is not encountered. The magnitude of the increment depends on the gene statistics (e.g., correlation of the gene with phenotype). The ES is the maximum deviation from zero encountered in the random walk; it corresponds to a weighted Kolmogorov-Smirnov(KS)-like statistic (Subramanian et al. 2005)." (from https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#gsea-algorithm)

Then the ES is normalized (NES) in order to compute a p-value.

 * in the base method [Subramanian 2005](https://www.pnas.org/doi/full/10.1073/pnas.0506580102) for each gene sets they create a number of permutated dataset for which they compute an ES, and they then compare the ES on the original data to the distribution of permutated ES for that set.
    2 flavours of permutation are described : "sample permutation" and "gene permutation". 
        * The sample permutation is only possible if the expression data for each sample is given to the method. It is recommended to have at least 7 samples per condition for it to make sense. 
        * The gene permutation is performed from the ranking metric directly. Hence it is sometimes called "preranked-GSEA" and to my knowledge this is the most often used permutation scheme of the 2.

 * in fGSEA[Korotkevich et al. 2021](https://www.biorxiv.org/content/10.1101/060012v3) several tricks are used to make the p-value computation of "preranked-GSEA" faster and more accurate for low p-values. Consequently it has become the default GSEA method in some libraries such as clusterProfiler.

Finally, a parameter to consider when performing (preranked-)(f)GSEA is which metric to use to rank genes. 
The most commons are logFC , -log10(pvalues) * logFC , or some signed test statistics (eg. t-test or Wald statistic).

