

Following a QC analysis on sequencing results, one may detect stretches of low quality bases along reads, or a contamination by adapter sequence.
Depending on your research question and the software you use for mapping, you may have to remove these bad quality / spurious sequences out of your data.


**During this block, you will learn to :**

 * trim your data with trimmomatic


## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/RNA-Seq_03_trimming.pdf){target=_blank : .md-button }

[Trimmomatic website](http://www.usadellab.org/cms/?page=trimmomatic){target=_blank : .md-button }


## to trim or not to trim ?

There are several ways to deal with poor quality bases or adapter contamination in reads, and several terms are used in the field, sometimes very loosely. We can talk about:

 * **Trimming**: to remove a part of, or the entirety of, a read (for quality reasons).
 	* **Hard trimming**: trim with a high threshold (eg. remove everything with QUAL<30).
 	* **Soft trimming**: trim with a low threshold (eg. remove everything with QUAL<10).
 * **Clipping**: to remove the end part of a read (typically because of adapter content).
 	* **Hard clipping**: actually removing the end of the read from the file (ie. with trimmomatic).
 	* **Soft clipping**: ignoring the end of the read at mapping time (ie. what STAR does).


If the data will be used to perform **transcriptome assembly, or variant analysis, then it MUST be trimmed**.


In contrast, for applications based on **counting reads**, such as **Differential Expression analysis**, most aligners, such as [STAR](https://github.com/alexdobin/STAR), [HISAT2](http://daehwankimlab.github.io/hisat2/), [salmon](https://salmon.readthedocs.io/en/latest/salmon.html), and [kallisto](https://pachterlab.github.io/kallisto/manual), can handle bad quality sequences and adapter content by soft-clipping, and consequently they _usually_ do not need trimming.
In fact, **trimming can be detrimental** to the number of successfully quantified reads \[[William et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2)\].

Nevertheless, it is usually recommended to perform some amount of soft trimming (*eg.* [kallisto](https://www.biostars.org/p/389324/), [salmon](https://github.com/COMBINE-lab/salmon/issues/398) ).

If possible, we recommend to perform the mapping for both the raw data and the trimmed one, in order to compare the results for both, and choose the best.

**Question:** what could be a good metric to choose the best between the trimmed and untrimmed?

??? success "Answer"

	The number of uniquely mapped reads is generally what would matter in differential expression analysis. Of course, this means that you can only choose after you have mapped both the trimmed and the untrimmed reads.




## trimming with Trimmomatic


The [trimmomatic website](http://www.usadellab.org/cms/?page=trimmomatic) gives very good examples of their software usage for both paired-end (`PE`) and single-end (`SE`) reads. We recommend you read their quick-start section attentively.



**Task 1:** 

 * Conduct a soft trimming on the mouseMT data

     - name the output folder : `030_d_trim/`.
     - Adapter sequences can be found in `/data/DATA/adapters/TruSeq3-PE.fa` ([adapter sequence source](https://github.com/usadellab/Trimmomatic/tree/main/adapters)).
     - unlike fastqc, you will have to launch trimmomatic for each sample separately
     - to facilitate QC afterward, add the following at the end of your trimmomatic command (substituting `<sample name>`):
		 		`2> 030_d_trim/trim_out.<sample name>.log`
		 		This will send part of the output of trimmomatic to a file in the same folder as the trimmed reads, which multiQC will be able to use afterward.


!!! note
	On real data Trimmomatic requires ~0.5G / CPU of RAM, and ~10 min/read file.

!!! warning

	trimmomatic is a Java-based program, and thus must be run by passing its .jar file to the Java interpreter:

	```{sh}
	trimmomatic
	```


??? success "trimmomatics script"

    We chose the following option:

     * **SLIDINGWINDOW:3:25** Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
       - **3**  : windowSize: specifies the number of bases to average across
       - **25** : requiredQuality: specifies the average quality required.
     * **ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10** Cut adapter and other Illumina-specific sequences from the read.
       - Cut adapter and other illumina-specific sequences from the read.
       - **2**  : seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
       - **30** : palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
       - **10** : simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.


    ```sh
    #!/usr/bin/bash
    # trimming the mouseMT fastq files


    ## creating output folder, in case it does not exists
    mkdir -p 030_d_trim

    ## we store the input folder in a variable, 
    ## to access its value, we will write $INPUT_FOLDER
    INPUT_FOLDER=/data/DATA/mouseMT

    ## by ending a line with \ we can continue the same command on the line below

    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_a1.fastq \
             030_d_trim/sample_a1.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_a1.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_a2.fastq \
             030_d_trim/sample_a2.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_a2.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_a3.fastq \
             030_d_trim/sample_a3.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_a3.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_a4.fastq \
             030_d_trim/sample_a4.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_a4.log
    
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_b1.fastq \
             030_d_trim/sample_b1.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_b1.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_b2.fastq \
             030_d_trim/sample_b2.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_b2.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_b3.fastq \
             030_d_trim/sample_b3.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_b3.log
    
    trimmomatic SE -phred33 \
             $INPUT_FOLDER/sample_b4.fastq \
             030_d_trim/sample_b4.trimmed.fastq \
             ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
             SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.sample_b4.log

    ```

    On the server, this script is also in `/data/Solutions/mouseMT/030_s_trim.sh`


??? success "alternative script using a for loop"

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
	it can also be found in the server at `/data/Solutions/mouseMT/sampleNames.txt`

	```
	#!/usr/bin/bash
    # trimming the mouseMT fastq files
	

	## creating output folder, in case it does not exists
	mkdir -p 030_d_trim
	
	INPUT_FOLDER=/data/DATA/mouseMT
	
	## each job grab a specific line from sampleNames.txt
	for SAMPLE in `cat sampleNames.txt`
	do
	 trimmomatic SE -phred33 \
                  $INPUT_FOLDER/${SAMPLE}.fastq \
                  030_d_trim/${SAMPLE}.trimmed.fastq \
                  ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
                  SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.${SAMPLE}.log
	done
	```

	On the server, this script is also in	`/data/Solutions/mouseMT/030bis_s_trim_loop.sh`


**Task 2:** 

 * Use the the following script to run a QC analysis on your trimmmed reads and compare with the raw ones.


```sh
#!/usr/bin/bash
# QC abalysis of the mouseMT trimmed reads

## fastQC on trimmed fastq files
fastqc 030_d_trim/*.fastq -o 030_d_trim


## multiqc on the fastQC reports AND the trimmomatic logs
multiqc -n 032_r_multiqc_mouseMT_trimmed.html -f --title trimmed_fastq 030_d_trim/
```
On the server, you can find this script in : `/data/Solutions/mouseMT/032_s_multiqc_trimmed.sh`


!!! note 
    The script above presumes that you have successfully trimmed the reads. 

    If not, you can grab them on the server in `/data/Solutions/mouseMT/030_d_trim/`



??? success "trimmed reads multiqc report"

    [:fontawesome-solid-file: Download the multiqc report](../assets/html/032_multiqc_mouseMT_trimmed.html){target=_blank : .md-button }

    **Note the second section, Trimmomatic, which lets you know the number/percentage of reads dropped**







**Extra Task: if you have the time**

We have trimmed the Liu2015 data with trimmomatic with the following option:

 * SLIDINGWINDOW:4:20 
 * MINLEN:36 
 * ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 

Download the multiQC report:

[:fontawesome-solid-file: Liu2015 trim multiqc report](../assets/html/032_r_trimQC_Liu2015.html){target=_blank : .md-button }

Analyze it : what do you think? Did the trimming do its job?
