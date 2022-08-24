# Precourse preparations

On top of a thirst for knowledge, and a working Internet connection, here is what you will need for the course : 

## NGS

As announced in the [course registration webpage](https://www.sib.swiss/training/course/20220901_IRNAS),
we expect participants to already have a basic knowledge in Next Generation Sequencing (NGS) techniques. 

## UNIX

Practical knowledge of the UNIX command line is also required to be able to follow this course, given that that the tools used to process sequenced reads use this interface.

If you are unsure about your capabilities or feel a bit rusty, we strongly recommend you spend some time practicing before the course : in our experience, the more comfortable you are with UNIX, the more you will be able to focus on the RNA-seq during the course, and the more you will gain from it.

You may refer to the [SIB's UNIX e-learning module](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html)

## R 

A basic knowledge of the [R language](https://www.r-project.org/) is required to perform most analytical steps after reads have been mapped and quantified : differential gene expression, gene set enrichment, over-representation analysis.
<!-- Can we point them towards an intro to R course? EG https://github.com/sib-swiss/first-steps-with-R-training -->

## Software

To replicate the technical condition of today's real-life data analysis, we will perform our computations on a distant HPC cluster.
To access it: 

 * macOS / Linux : you can use your pre-installed terminal.
 * Windows : you should install a terminal which lets you do ssh (for instance [mobaXterm](https://mobaxterm.mobatek.net/)). 

Additionally, a graphical client for file transfer to and from the distant server can be useful. MobaXterm integrates this functionality, so if you use it there is no need for additional software. 
Otherwise, we recommend [FileZilla](https://filezilla-project.org/download.php?type=client).
