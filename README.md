# Introduction to Mammalian Sequence Analysis
Sergio Llaneza-Lago

# Setting up

## Bash terminal

#### Accessing Your Workshop Environment

1.  **Follow [this
    link](https://github.com/UEA-Cancer-Genetics-Lab/MHC_RNA-seq_Workshop)
    to GitHub.**

2.  **Log in to GitHub** if you are not already.

3.  Find and click the green **`< > Code`** button.

4.  In the pop-up window, select the **`Codespaces`** tab.

5.  Click on **`Create codespace on main`**.

6.  A new browser window will open and begin loading your terminal
    environment.

#### Setting Up Your Workspace

Once the Codespace has fully loaded and you see the terminal in the
bottom pane:

1.  Type the following command into the terminal:

``` bash
bash download_reference.sh
```

2.  Press enter.

## Script to install required libraries in Rstudio

1.  Open Rstudio.
2.  Copy and paste the following and run it:

``` r
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[!req]

  if (length(need) > 0) {
    # Install BiocManager if missing
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    # Try installing from CRAN first, then Bioconductor
    for (pkg in need) {
      tryCatch({
        install.packages(pkg)
      }, error = function(e) {
        message(paste("Trying Bioconductor for:", pkg))
        BiocManager::install(pkg, ask = FALSE)
      })
      library(pkg, character.only = TRUE)
    }
  }
}

# Load/install the required packages
using("tidyverse", "tximport")
```

# RNA-seq analysis

The goal of RNA-seq is often to perform differential expression testing
to determine which genes or transcripts are expressed at different
levels between conditions. These findings can offer biological insight
into the processes affected by the condition(s) of interest. Below is an
overview of the typical analysis workflow that is followed for
differential gene expression analysis with bulk RNA-seq data.

![RNA-seq workflow and examples of software for each
step.](img/workflow.png)

This session is divided into two parts. In the first part, we will cover
the process of transforming our RNA sequences into counts using a Bash
terminal. During the second part, we will use R to identify
differentially expressed genes and perform functional analysis.

## How to use Bash?

### 1. Make a Folder

Use `mkdir` to make a folder (called a “directory”) to store your files.

``` bash
mkdir rnaseq_project
```

### 2. Move Into That Folder

Use `cd` to change directory and go into it.

``` bash
cd rnaseq_project
```

Use `cd ..` to go back one level.

### 3. See What’s Inside

Use ls to list files and folders.

``` bash
ls
```

If you want more detail, try:

``` bash
ls -l
```

### 4. Open compressed files

Use `zcat` to see the content of a file

``` bash
zcat data/file.fq.gz
```

If the file is too big, you can combine `zcat` with `head` to get the
first ten lines:

``` bash

zcat data/file.fq.gz | head
```

## 1. From Sequence to Counts Using Bash

### 1.1 The FASTQ file

The first step in the RNA-Seq workflow is to take the FASTQ files
received from the sequencing facility (e.g., Novogene) and assess the
quality of the reads.

The FASTQ file format is the defacto file format for sequence reads
generated from next-generation sequencing technologies. This file
contains sequence data, but also contains quality information. For a
single record (sequence read) there are four lines, each of which are
described below:

| Line | Description                                                                                            |
|------|--------------------------------------------------------------------------------------------------------|
| 1    | The header line. Always begins with ‘@’ and then shows information about the read                      |
| 2    | The actual DNA sequence                                                                                |
| 3    | Always begins with a ‘+’ and sometimes repeats info from line 1                                        |
| 4    | String of characters which represent the quality scores; must have same number of characters as line 2 |

Let’s use the following read as an example:

    @HWI-ST330:304:H045HADXX:1:1101:1111:61397
    CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
    +
    @?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################

As mentioned previously, line 4 has characters encoding the quality of
each nucleotide in the read. The legend below provides the mapping of
quality scores (Phred-33) to the quality encoding characters.

    Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                      |         |         |         |         |
       Quality score: 0........10........20........30........40                                

### 1.2 Assessing quality with FastQC

Now we understand what information is stored in a FASTQ file, the next
step is to examine quality metrics for our data.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
provides a simple way to do some quality control checks on raw sequence
data coming from high throughput sequencing pipelines which you can use
to give a quick impression of whether your data has any problems of
which you should be aware before doing any further analysis.

#### 1.2.a Run FastQC

Below there is an example of how to use FastQC for paired reads from one
sample:

``` bash
fastqc mySample_1.fq.gz mySample_2.fq.gz -o outputDirectory/
```

FASTQ files are usually formatted as `fq` and are compressed into the
`gzip` format (denoted as `gz`) to save space. FastQC can be applied
directly to the compessed FASTQ files. The `-o` flag denotes where to
store the output from the FastQC analysis.

#### 1.2.b Understanding FastQC results

Let’s take a closer look at the report that FastQC generates. One thing
to keep in mind is that ***FastQC is just an indicator of what’s going
on with your data, don’t take the “PASS”es and “FAIL”s too seriously.***

FastQC has a really well documented [manual
page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with
[more
details](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)
about all the plots in the report. We recommend looking at [this
post](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010)
for more information on what bad plots look like and what they mean for
your data.

Below are two of the most important analysis modules in FastQC, the
**“Per base sequence quality”** plot and the **“Overrepresented
sequences”** table.

The **“Per base sequence quality”** plot provides the distribution of
quality scores across all bases at each position in the reads.

![Sequence quality per base in a good quality
sample.](img/fastqc_good.png)

![Sequence quality per base in a poor quality
sample.](img/fastqc_poor.png)

The **“Overrepresented sequences”** table displays the sequences (at
least 20 bp) that occur in more than 0.1% of the total number of
sequences. This table aids in identifying contamination, such as vector
or adapter sequences.

![FastQC Overrepresented sequences example](img/over_seq.png)

### Exercises \#1

1.  For this first part of the session, we will use a paired-end
    chondrosarcoma sample located in `data/`. Examine the FASTQ files
    with `zcat`.
2.  Perform FastQC on both reads of the sample and store the output in
    `results/fastqc/`.
3.  Explore the FastQC report, are these good quality samples?

### 1.3 To count or to pseudocount?

Traditional RNA-seq pipelines follow a specific path: After quality
control (QC), sequencing reads are trimmed to remove adapters and
low-quality bases. These cleaned reads are then aligned base-by-base to
a reference genome to pinpoint their original genomic locations.

<img src="img/RNAseqWorkflowTrad.png" width="404"
alt="Standard RNA-seq workflow performing trimming and alignment" />

This alignment-based method is computationally demanding and takes a lot
of time. **Another strategy for quantification which has more recently
been introduced involves transcriptome mapping**. Tools that fall in
this category include
[Kallisto](https://pachterlab.github.io/kallisto/about),
[Sailfish](http://www.nature.com/nbt/journal/v32/n5/full/nbt.2862.html)
and [Salmon](https://combine-lab.github.io/salmon/); each working
slightly different from one another.

Common to all of these tools is that **base-to-base alignment of the
reads is avoided**, and these tools **provide quantification estimates
much faster than do standard approaches** (typically more than 20 times
faster) with **improvements in accuracy** at **the transcript level**.
These transcript expression estimates, often referred to as
‘pseudocounts’, can be converted for use with DGE tools like DESeq2 or
the estimates can be used directly for isoform-level differential
expression using a tool like
[Sleuth](http://www.biorxiv.org/content/biorxiv/early/2016/06/10/058164.full.pdf).

![Aligment free alternative using Salmon and
pseudocounts.](img/alignmentfree_workflow_aug2017.png)

In this workshop we will explore Salmon in more detail.

### 1.4 What is Salmon?

[Salmon](http://salmon.readthedocs.io/en/latest/salmon.html#using-salmon)
is a tool that takes a **reference transcriptome** (in FASTA format) and
raw **sequencing reads** (in FASTQ format) as input. Unlike traditional
methods, it does *not* align full reads. Instead, Salmon performs both
mapping and quantification efficiently. It is known for being extremely
fast at “mapping” reads to the transcriptome and often provides more
accurate results than standard alignment-based approaches.

Let’s break down how Salmon works:

![](img/salmon_workflow_subset.png)

#### **1.4.a Indexing**

This step involves creating an index to evaluate the sequences for all
possible unique sequences of length k (kmer) in the **transcriptome**
(genes/transcripts).

**The index helps creates a signature for each transcript in our
reference transcriptome.** The Salmon index has two components:

- a suffix array (SA) of the reference transcriptome
- a hash table to map each transcript in the reference transcriptome to
  it’s location in the SA (is not required, but improves the speed of
  mapping drastically)

#### **1.4.b Quasi-mapping and quantification**

The quasi-mapping approach estimates the numbers of reads mapping to
each transcript, then generates final transcript abundance estimates
after modeling sample-specific parameters and biases.

> **NOTE:** If there are sequences in the reads that are not in the
> index, they are not counted. As such, trimming is not required when
> using this method. Accordingly, if there are reads from transcripts
> not present in the reference transcriptome, they will not be
> quantified. Quantification of the reads is only as good as the quality
> of the reference transcriptome.

![Salmon quasi-mapping process](img/salmon_quasialignment.png)

### 1.5 Running Salmon

As you can imagine from the description above, when running Salmon there
are also two steps.

**Step 1: Indexing** “Index” the transcriptome using the `index`
command:

``` bash
salmon index -t transcripts.fa -i transcripts_index --type quasi -k 31
```

> **NOTE:** Default for salmon is –type quasi and -k 31, so we do not
> need to include these parameters in the index command. The kmer
> default of 31 is optimized for 75bp or longer reads, so if your reads
> are shorter, you may want a smaller kmer to use with shorter reads
> (kmer size needs to be an odd number).

**We are not going to run this in class, but it only takes a few
minutes.** We will be using an index we have generated from transcript
sequences (all known transcripts/ splice isoforms with multiples for
some genes) for human (GRCh37). The transcriptome data (FASTA) was
obtained from the [GENCODE
website.](https://www.gencodegenes.org/human/release_19.html)

**Step 2: Quantification:** Get the transcript abundance estimates using
the `quant` command and the parameters described below (more information
on parameters can be found
[here](http://salmon.readthedocs.io/en/latest/salmon.html#id5)):

- **`-i`:** specify the location of the index directory
  (`reference_files/salmon_index/`)
- **`-l`:** library type (single-end reverse reads `SR`, but we can use
  `A` to automatically infer the library type) - more information is
  available
  [here](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype))
- **`-o`:** output quantification file name (`Mov10_oe_1.subset.salmon`)
- **`--writeMappings`:** instead of printing to screen, write to a file
  (`--writeMappings=salmon.out`)
- `-p`: Number of threads to use (how much computational power to
  dedicate for this task)

To run the quantification step on a single sample we have the command
provided below. Let’s try running it on our subset sample for
`Mov10_oe_1.subset.fq`:

``` bash
% salmon quant -i reference_files/salmon_index/ \
 -l A \
 -1 mySample_1.fq.gz \
 -2 mySample_2.fq.gz \
 -o results/salmonQuant \
 --writeMappings \
 --gcBias \
 -p 4
```

> **NOTE:** To correct for the various sample-specific biases you could
> add the following parameters to the Salmon command:
>
> - `--seqBias` will enable it to learn and correct for
>   sequence-specific biases in the input data.
> - `--gcBias` to learn and correct for fragment-level GC biases in the
>   input data
> - `--posBias` will enable modeling of a position-specific fragment
>   start distribution

### Exercise \#2

1.  Perform `salmon quant` on the chondrosarcoma samples. Use the same
    arguments than the example when appropriate.

2.  After Salmon finishes running, you should see a new directory
    created with the name you provided in the `-o` parameter. **Take a
    look inside this directory.**

3.  Have you noticed the file named `quant.sf`? This is the primary
    quantification file. Each row in `quant.sf` corresponds to a
    transcript (identified by its Ensembl ID), and the columns list
    various quantification metrics for each transcript. **Explore this
    file to understand its contents.**

### 1.6 Salmon output

Salmon creates a directory containing a logs directory, which contains
all of the text that was printed to screen as Salmon was running.
Additionally, there is a file called `quant.sf`.

This is the **quantification file** in which each row corresponds to a
transcript, listed by Ensembl ID, and the columns correspond to metrics
for each transcript:

``` bash
Name    Length  EffectiveLength TPM     NumReads
ENST00000632684.1       12      3.00168 0       0
ENST00000434970.2       9       2.81792 0       0
ENST00000448914.1       13      3.04008 0       0
ENST00000415118.1       8       2.72193 0       0
ENST00000631435.1       12      3.00168 0       0
ENST00000390567.1       20      3.18453 0       0
ENST00000439842.1       11      2.95387 0       0

....
```

- The first two columns are self-explanatory, the **name** of the
  transcript and the **length of the transcript** in base pairs (bp).
- The **effective length** represents the the various factors that
  effect the length of transcript due to technical limitations of the
  sequencing platform.
- Salmon outputs ‘pseudocounts’ which predict the relative abundance of
  different isoforms in the form of three possible metrics (KPKM, RPKM,
  and TPM). **TPM (transcripts per million)** is a commonly used
  normalization method and is computed based on the effective length of
  the transcript.
- Estimated **number of reads** (an estimate of the number of reads
  drawn from this transcript given the transcript’s relative abundance
  and length).

### Exercise \#3

1.  Download the `quant.sf` to your computer. We will not need to use
    the bash terminal anymore in this session, so feel free to close it.

2.  In Rstudio, follow these steps:

``` r
library(tximport)

tx2knowngene <- read_csv(...)

txi <- tximport("quant.sf",
                type = "salmon",
                tx2gene = tx2knowngene,
                ignoreTxVersion = TRUE)
```

> **NOTE:** Make sure that `quant.sf` is in your current working
> directory.

3.  The previous command will import the quantification files and map it
    to genes. Which gene has the highest number of counts?

------------------------------------------------------------------------

*Part of the material of this workshop was adapted from Dr Sarah Boswell
from Laboratory of Systems Pharmacology, Single Cell Core; and Dr
Radhika Khetani ,Training Director at the Harvard Chan Bioinformatics
Core RNA-seq workshop available at
https://hbctraining.github.io/rnaseq-cb321/.*
