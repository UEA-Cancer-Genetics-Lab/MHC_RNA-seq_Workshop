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
using("tidyverse", "tximport", "DESeq2", "EnhancedVolcano", "gprofiler2")
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

## 1. How to use Bash?

### 1.1 Make a Folder

Use `mkdir` to make a folder (called a “directory”) to store your files.

``` bash
mkdir rnaseq_project
```

### 1.2 Move Into That Folder

Use `cd` to change directory and go into it.

``` bash
cd rnaseq_project
```

Use `cd ..` to go back one level.

### 1.3 See What’s Inside

Use `ls` to list files and folders.

``` bash
ls
```

If you want more detail, try:

``` bash
ls -l
```

### 1.4 Open compressed files

Use `zcat` to see the content of a file

``` bash
zcat data/file.fq.gz
```

If the file is too big, you can combine `zcat` with `head` to get the
first ten lines:

``` bash

zcat data/file.fq.gz | head
```

## 2. From Sequence to Counts Using Bash

### 2.1 The FASTQ file

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

### 2.2 Assessing quality with FastQC

Now we understand what information is stored in a FASTQ file, the next
step is to examine quality metrics for our data.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
provides a simple way to do some quality control checks on raw sequence
data coming from high throughput sequencing pipelines which you can use
to give a quick impression of whether your data has any problems of
which you should be aware before doing any further analysis.

#### 2.2.a Run FastQC

Below there is an example of how to use FastQC for paired reads from one
sample:

``` bash
fastqc mySample_2.fq.gz mySample_2.fq.gz -o outputDirectory/
```

FASTQ files are usually formatted as `fq` and are compressed into the
`gzip` format (denoted as `gz`) to save space. FastQC can be applied
directly to the compessed FASTQ files. The `-o` flag denotes where to
store the output from the FastQC analysis.

#### 2.2.b Understanding FastQC results

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

### 2.3 To count or to pseudocount?

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

### 2.4 What is Salmon?

[Salmon](http://salmon.readthedocs.io/en/latest/salmon.html#using-salmon)
is a tool that takes a **reference transcriptome** (in FASTA format) and
raw **sequencing reads** (in FASTQ format) as input. Unlike traditional
methods, it does *not* align full reads. Instead, Salmon performs both
mapping and quantification efficiently. It is known for being extremely
fast at “mapping” reads to the transcriptome and often provides more
accurate results than standard alignment-based approaches.

Let’s break down how Salmon works:

![](img/salmon_workflow_subset.png)

#### **2.4.a Indexing**

This step involves creating an index to evaluate the sequences for all
possible unique sequences of length k (kmer) in the **transcriptome**
(genes/transcripts).

**The index helps creates a signature for each transcript in our
reference transcriptome.** The Salmon index has two components:

- a suffix array (SA) of the reference transcriptome
- a hash table to map each transcript in the reference transcriptome to
  its location in the SA (is not required, but improves the speed of
  mapping drastically)

#### **2.4.b Quasi-mapping and quantification**

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

### 2.5 Running Salmon

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
- **`-o`:** output quantification file name
- **`--writeMappings`:** instead of printing to screen, write to a file
  (`--writeMappings=salmon.out`)
- `-p`: Number of threads to use (how much computational power to
  dedicate for this task)

To run the quantification step on a single sample we have the command
provided below.

``` bash
salmon quant -i reference_files/salmon_index/ \
 -l A \
 -1 mySample_2.fq.gz \
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

### 2.6 Salmon output

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
library(tidyverse)

tx2knowngene <- read_csv("https://raw.githubusercontent.com/UEA-Cancer-Genetics-Lab/MHC_RNA-seq_Workshop/refs/heads/main/files_for_R/tx2knownGene.csv", col_select = c(2,3))

txi <- tximport("quant.sf",
                type = "salmon",
                tx2gene = tx2knowngene,
                ignoreTxVersion = TRUE)
```

> **NOTE:** Make sure that `quant.sf` is in your current working
> directory.

3.  The previous command will import the quantification files and map it
    to genes. Which gene has the highest number of counts?

## 3. From Counts to Differential Expression and Enrichment

In this second part of the workshop, we will systematically walk through
an RNA-seq differential expression workflow using various R packages.

Our starting point will be gene count data from both normal and tumor
tissue samples derived from prostate cancer patients. We will then:

1.  **Perform differential expression analysis** to identify genes that
    are significantly up- or down-regulated between the two conditions.

2.  **Create a volcano plot** to visually represent our differentially
    expressed genes.

3.  **Perform enrichment analysis** on our identified genes of interest
    to uncover their associated biological pathways or functions.

4.  **Visualize our enrichment results** for better interpretation.

### 3.1 The Prostate Cancer (PRAD) dataset

We will be using publicly available prostate cancer expression data from
**The Cancer Genome Atlas (TCGA)**. This dataset is included in the
workshop repository and can be loaded directly into R.

#### Step 1: Load Expression Data

To load the gene expression data, run the following command in your R
environment:

``` r
pradExpression <- read_csv("https://raw.githubusercontent.com/UEA-Cancer-Genetics-Lab/MHC_RNA-seq_Workshop/refs/heads/main/files_for_R/PRAD_expression.csv")
```

This `pradExpression` dataset contains:

- **108 samples** (as columns).

- **20,828 genes** (as rows).

- Gene names are provided in the column named `SYMBOL`.

#### Step 2: Load Sample Information

In addition to the expression data, we also have detailed information
about each sample. Download this by running:

``` r
pradSampleInfo <- read_csv("https://raw.githubusercontent.com/UEA-Cancer-Genetics-Lab/MHC_RNA-seq_Workshop/refs/heads/main/files_for_R/PRAD_SampleInfo.csv")
```

This `pradSampleInfo` dataset contains:

- **108 rows**, with each row representing a unique sample.

- **Three columns:**

  - `PatientID`: A unique identifier for each patient.

  - `SampleID`: A unique identifier for each sample.

  - `SampleType`: Indicates if the sample is from **tumoural** tissue
    (denoted as ‘Tumour’) or **normal** tissue.

You will observe that for all patients, we have paired samples: one
normal and one tumoural.

### 3.2 Differential expression analysis with DESeq2

**DESeq2** is a widely used R package for gene-level differential
expression analysis. It employs the negative binomial distribution and
is known for balancing sensitivity and specificity, helping to reduce
both false positives and false negatives in your results. For more
details and helpful tips, refer to the [DESeq2
vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

#### 3.2.a DESeq2 input

We can use the function *`DESeqDataSetFromMatrix`* for our already
prepared read count matrix. For this, we require three inputs:

1.  A **count matrix** (as an R `matrix` object).
2.  A **sample information data frame** (where rows correspond to
    samples and columns contain metadata).
3.  A **design formula**.

**Important:** The row names of your sample information data frame *must
exactly match* the column names of your count matrix.

A design formula tells the statistical software which sources of
variation to test for. This includes your primary factor of interest
(e.g., treatment vs. control) and any other known covariates that
contribute significantly to variation in your data (e.g., sex, age,
batch effects). **The design formula should have all of the factors in
your metadata that account for major sources of variation in your
data.**

For example, suppose we have the following sample info:

![](img/meta_example.png)

If you want to examine the expression differences between treatments,
your design formula would be:

`design = ~ treatment`

However, if you know that `sex` and `age` are also major sources of
variation, your formula should include them to account for their
effects:

`design = ~ sex + age + treatment`

The tilde (`~`) should always precede your factors and tells DESeq2 to
model the counts using the following formula. Note the **factors
included in the design formula need to match the column names in the
metadata**.

#### 3.2.b Preparing our data

Currently, both our `pradExpression` (counts) and `pradSampleInfo`
(sample information) are data frames, and neither has appropriate row
names for DESeq2. Let’s fix this.

##### Exercise \# 4

1.  Complete and run the following code:

    ``` r
    expressionMatrix <- [Type here your count matrix file] |>
      column_to_rownames([Type here the column to convert into row names, use quotation marks]) |>
      as.matrix()
    ```

2.  Modify `pradSampleInfo` so the column `SampleID` becomes the row
    names. Keep it as a data frame.

3.  All data is ready, time to create our `DESeq2` object. Complete and
    run the following code using the sample type as design:

``` r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = [Type here your count matrix],
                              colData = [Type here your sample info dataframe],
                              design = ~ [Add here the factor to compare by])
```

#### 3.2.c Running DESeq2

In the previous section, we created a `DESeq2` object in our R
environment. To run the analysis just run this code:

``` r
des <- DESeq(dds)
```

`DESeq2` will start to analyse your data and print the progress to your
console. You’ll notice that the analysis performs several steps. We will
be taking a detailed look at each of these steps to better understand
how DESeq2 is performing the statistical analysis and what metrics we
should examine to explore the quality of our analysis.

![](img/deseq2_workflow_separate.png)

##### Step 1: Estimate size factors

The first step in the differential expression analysis is to estimate
the size factors. In RNA-seq experiments, the raw counts we obtained
depend on how much the gene is expressed and how much sequencing was
done. This means that some samples may have more reads purely for
technical reasons, not because of biological differences.

To correct this, DESeq2 calculates a size factor for each sample to
adjust for differences in sequencing depth and scales them accordingly.
This way we can ensure that we compare biologically meaningful
differences.

In the example below, imagine the sequencing depths are similar between
Sample A and Sample B, and every gene except for gene DE presents
similar expression level between samples. The counts in Sample B would
be greatly skewed by the DE gene, which takes up most of the counts.
Other genes for Sample B would therefore appear to be less expressed
than those same genes in Sample A.

![](img/normalization_methods_composition.png)

##### Step 2: Estimate gene-wise dispersion

Next, DESeq2 estimates the dispersion for each gene. Dispersion measures
how much a gene’s expression varies between replicates, after correcting
for sequencing depth.

A gene with stable counts across samples has low dispersion. A gene with
more fluctuation has high dispersion. This information is used to decide
how much confidence to place in any differences observed between groups.

For example, if we imagine two genes:

- Gene A: Counts are always around 100 across samples. It has low
  dispersion.

- Gene B: Sometimes the counts are 10, sometimes is 500. It has high
  dispersion.

Even if both genes have similar average counts, Gene B is less
consistent and needs stronger evidence to be considered significantly
differentially expressed.

##### Step 3: Fit curve to gene-wise dispersion estimates

The next step in the workflow is to fit a curve to the gene-wise
dispersion estimates. The idea here is that when dispersion is plotted
against mean expression for all genes, a clear trend usually appears:
genes with lower expression tend to show higher variability.

DESeq2 fits a smooth curve through these points to model this
relationship. The curve reflects the expected dispersion for a gene
based on its average expression level.

<img src="img/deseq_dispersion1.png" width="400"
alt="In this plot we have dispersion on the y-axis and mean normalized counts on the x-axis. Each black dot represents a gene. Simply looking at the trend of black dots, we observe an inverse relationship between mean and dispersion." />

##### Step 4: Shrink gene-wise dispersion estimates toward the values predicted by the curve

The raw dispersion estimates, especially for lowly expressed genes, are
often unreliable. DESeq2 applies a method called shrinkage, which
adjusts these estimates towards the fitted curve.

This improves the accuracy of the analysis. Noisy or extreme values are
moderated, while more stable genes remain largely unaffected. The extent
of shrinkage depends on how far a gene’s dispersion is from the expected
trend and how many samples are in the dataset.

Shrinkage is important for reducing false positives, particularly when
working with small sample sizes.

![Example of genes (black dots) being pulled to the curve. Genes with a
dispersion very away from the curve are not pulled as they are assumed
that their values hold biological value.](img/deseq_dispersion2.png)

##### Step 5: Comparing groups

Once dispersions have been adjusted, DESeq2 compares gene expression
between groups.

- If comparing two groups (e.g. High vs Low risk), DESeq2 uses the Wald
  test.

- If comparing more than two groups, it uses the Likelihood Ratio Test
  (LRT).

Each gene is assigned a p-value. However, when testing thousands of
genes, some may appear significant by chance. This is known as the
multiple testing problem.

DESeq2 corrects for this by adjusting the p-values using the False
Discovery Rate (FDR), based on the Benjamini-Hochberg method. By
default, it uses a threshold of 0.05. This means that among the genes
called significant, no more than 5% are expected to be false positives.

For instance, if 500 genes are declared differentially expressed at FDR
\< 0.05, we expect around 25 of those to be false positives.

#### 3.2.c. DESeq2 results

After running the DESeq2 analysis, we’ll use the `results()` function to
extract the outcomes.

Run the following command to get your differential expression results:

``` r
deResults <- as.data.frame(results(des, name = "SampleType_Tumour_vs_Normal")) |>
  rownames_to_column('SYMBOL')
```

> **Note on `name` parameter:** In this analysis, we are only comparing
> two groups (`Tumour` vs. `Normal` within `SampleType`), so DESeq2 will
> by default return this specific comparison. Therefore, explicitly
> setting `name = "SampleType_Tumour_vs_Normal"` is good practice for
> clarity, but not strictly necessary here. However, if you were
> comparing three or more groups, or if your design included additional
> cofactors, you **would need to specify** which comparison you want to
> extract.

This command will produce a data frame named `deResults`, where each row
corresponds to one of the genes we tested for differential expression.
The columns provide key metrics for each gene:

- **`baseMean`**: The average of the normalised count values, taken over
  all samples. This represents the gene’s overall expression level.

- **`log2FoldChange`**: This value indicates the magnitude and direction
  of expression change.

  - A **positive** `log2FoldChange` means the gene is overexpressed
    (up-regulated) in tumour samples compared to normal samples.

  - A **negative** `log2FoldChange` means the gene is underexpressed
    (down-regulated) in tumour samples compared to normal samples.

  - **Interpretation:** A `log2FoldChange` of 1 means the gene is 2
    times more expressed. A value of 2 means it’s 4 times more
    expressed, and so on.

- **`lfcSE`**: The standard error of the `log2FoldChange` estimate.

- **`stat`**: The Wald test statistic, which is the `log2FoldChange`
  divided by its `lfcSE`.

- **`pvalue`**: The p-value obtained from the Wald test, indicating the
  probability of observing such a `log2FoldChange` by chance if there
  were no true difference.

- **`padj`**: The adjusted p-value (using the Benjamini-Hochberg (BH)
  method).

### Exercise \#5

1.  Plot the gene-wise dispersion curve by running:

    ``` r
    plotDsipEsts(des)
    ```

    Does the dispersion follow the curve?

2.  Create a new data frame containing only the significantly
    differentially expressed genes. This can be done by **filtering** in
    the genes in `deResults` with `padj < 0.05`.

3.  Find the 5 genes significantly more overexpressed and the five more
    underexpressed. Search a couple on Google to find what is their role
    in prostate cancer. Do they make biological sense?

### 3.3 Visualising our DESeq2 output with `EnhancedVolcano`

When we are working with large amounts of data it can be useful to
display that information graphically to gain more insight. Visualization
deserves an entire course of its own, but during this lesson we will get
you started with volcano plots.

Volcano plots show the log transformed adjusted p-values plotted on the
y-axis and log2 fold change values on the x-axis. They are useful
because they give an overall distribution of the expression of our
genes. There is no built-in function for the volcano plot in DESeq2, so
instead we will use the `EnhancedVolcano` package. This package uses
`ggplot2` as a base to create the plots so it allows for a lot of
customisation, more information can be [found
here](https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html).

It can be used this way:

``` r
library(EnhancedVolcano)

EnhancedVolcano([Add here the results dataframe we got from DEseq2],
                lab = [Column from the results containing gene names as results$geneNames],
                x = [Name of the column from the results containing the fold change values, needs to be within quotation marks, e.g. 'fold'],
                y = [Name of the column from the results containing the adjusted p-values, needs to be within quotation marks],
                title = 'PRAD Tumour vs Normal',
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                pointSize = 1.5,
                labSize = 3.0)
```

### Exercise \#6

Use `EnhancedVolcano` to visualise your differentially expressed genes.

### 3.4 Functional enrichment analysis with `gProfiler2`

The output of RNA-seq differential expression analysis is a list of
significant differentially expressed genes (DEGs). To gain greater
biological insight on the differentially expressed genes there are
various analyses that can be done:

- determine whether there is enrichment of known biological functions,
  interactions, or pathways

- identify genes’ involvement in novel pathways or networks by grouping
  genes together based on similar trends

- use global changes in gene expression by visualizing all genes being
  significantly up- or down-regulated in the context of external
  interaction data

Generally for any differential expression analysis, it is useful to
interpret the resulting gene lists using freely available web-and
R-based tools. In this session, we will use the R package `gProfiler2`
to perform functional enrichment through the Gene Ontology (GO)
database, the Kyoto Encyclopedia of Genes and Genomes (KEGG) database,
and the Reactome Pathway Database.

***Note that all tools described below are great tools to validate
experimental results and to make hypotheses. These tools suggest
genes/pathways that may be involved with your condition of interest;
however, you should NOT use these tools to make conclusions about the
pathways involved in your experimental process. You will need to perform
experimental validation of any suggested pathways.***

#### 3.4.a How does Functional Enrichment work?

These three databases (KEGG, GO, and Reactome) categorise genes into
groups (gene sets) based on a shared function, or involvement in a
pathway, or presence in a specific cellular location, or other
categorisations such as functional pathways.

For each gene set, `gProfiler2` calculates the probability that the
observed number of your significantly differentially expressed genes
within that set is higher than expected by chance. This calculation
compares your gene list to a “background” set (in our case, the entire
human transcriptome). The statistical significance of this
over-representation is determined using the **hypergeometric test**.

![](img/go_proportions.png)

![](img/go_proportions_table3.png)

Each of the three databases you will use today have their own system to
categorise genes.

- GO is helpful for understanding the biological processes that
  genes/products are involved in. It is subdivided into three
  ‘ontologies’:

  - Biological processes (BP): refers to the biological role involving
    the gene or gene product, and could include “transcription”, “signal
    transduction”, and “apoptosis”.

  - Cellular component (CC): refers to the location in the cell of the
    gene product. Cellular components could include “nucleus”,
    “lysosome”, and “plasma membrane”.

  - Molecular function (MF): represents the biochemical activity of the
    gene product, such activities could include “ligand”, “GTPase”, and
    “transporter”.

- KEGG is more appropriate for understanding the precise pathway where
  the gene/product is active. It also contains various biological
  systems and diseases but requires licensing for certain access.

- Reactome is a manually curated database primarily focused on human
  biological processes. It is completely open-source and contains
  cross-links to other databases.

#### 3.4.b How to use `gprofiler2`?

To use `gprofiler2`, the first step is to query for gene sets associated
with your differentially expressed genes. It’s a good practice to
separate your genes into overexpressed and underexpressed lists, as this
makes interpreting the enrichment results much easier.

Additionally, while optional, providing your genes ranked from most to
least differentially expressed (`log2FoldChange`) significantly improves
the accuracy and biological relevance of the enrichment results. More
details about `gprofiler2` can be found in its
[documentation](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html).

Here’s an example of how to perform a `gprofiler2` query:

``` r
library(gprofiler2)

gostres <- gost(query = deResults$SYMBOL, 
                   organism = "hsapiens",
                   significant = TRUE,
                   sources = c("REAC", "KEGG", "GO:BP"),
                   correction_method = "fdr",
                   ordered_query = TRUE)
```

Let’s break down the parameters used in this `gost()` function:

- **`query`**: This is your list of gene identifiers (e.g., gene
  symbols). For best results, this list should be sorted by a measure of
  differential expression, such as `log2FoldChange` (from most to least
  differentially expressed).

- **`organism`**: Specifies the organism being analysed. For human, use
  `"hsapiens"`.

- **`significant`**: If set to `TRUE`, only significant enrichment
  results (based on the adjusted p-value threshold) will be returned.

- **`sources`**: A character vector indicating which databases to query.
  For this session, we will use:

  - `"REAC"` for Reactome pathways.

  - `"KEGG"` for KEGG pathways.

  - `"GO:BP"` for Gene Ontology Biological Process terms.

- **`correction_method`**: Defines the method used to adjust the
  p-values obtained from the hypergeometric test for multiple testing.
  `"fdr"` (False Discovery Rate) is a common and robust choice.

- **`ordered_query`**: Set to `TRUE` if your `query` gene list is sorted
  (e.g., by `log2FoldChange`). This allows `gprofiler2` to perform a
  more sensitive and specific ranked enrichment analysis.

You can access the full results of your query by inspecting the
`gostres$results` object:

``` r
gostres$results
```

#### 3.4.c Visualising `gprofiler2` results

We can use the function `gostplot` to visualise all our enrichment
results in a Manhattan plot by running this:

``` r
gostplot(gostres, capped = TRUE, interactive = TRUE)
```

> **NOTE:** You will need to set `interactive` to `FALSE` if you want to
> edit this plot in the future.

![Manhattan plot generated by \`gostplot\`. Each circle corresponds to a
significantly enriched gene set. The y-axis denotes how much enriched
they are. In this image, one biological process is higlighted by
clicking on it.](img/Screenshot%202025-07-10%20at%2011.56.33.png)

It is possible to pick specific gene sets to highlight on the plot using
the function `publish_gostplot`. Here is an example highlighting sets
that we found related to muscle contraction:

``` r
gplot <- gostplot(gostres, capped = TRUE, interactive = FALSE)

publish_gostplot(gplot,
                 highlight_terms = c("KEGG:04260", "REAC:R-HSA-390522", "REAC:R-HSA-397014"))
```

![Manhattan plot highlighting the selected gene sets. At the bottom, it
includes a table with the selected gene sets, the number of genes
included in the set (term_size) and the p-value from the geometric
test.](img/Screenshot%202025-07-10%20at%2012.05.50.png)

##### Visualising using `ggplot2`

Another alternative is to create our own plots using the `ggplot2`
package. A common way to represent gene sets is using barplots, for
example:

``` r
ggplot(gostres$result, aes(x = intersection_size, y = reorder(term_name, intersection_size), fill = p_value)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low = "blue", high = "red")
```

> **NOTE:** The `reorder` function ensures that the bars are sorted from
> higher to smaller. It is optional but it makes the plot to look much
> better.

![](img/Screenshot%202025-07-10%20at%2014.41.07.png)

### Exercise \#7

1.  Create two new data frames, one containing only the significant
    upregulated genes from `deResult` and another only the significant
    downregulated genes. Sort them by the log2 Fold Change.

2.  Use the `gost` function on both sets of data, and explore the
    outcome.

3.  Use `gostplot` to visualise your enrichment results.

4.  Use `publish_gostplot` and highlight the ten upregulated Reactome
    pathways with the biggest `intersection_size`.

5.  Using `ggplot2` create a barplot of the upregulated Reactome
    pathways.

6.  Using your previous barplot as a starting point:

    - Change the barplot to a scatterplot.

    - Make the dots on the scatterplot to be coloured by the `p_value`.

    - Make the size of the dots at the scatterplot to be dependant on
      the `intersection_size`.

------------------------------------------------------------------------

*Part of the material of this workshop was adapted from Dr Sarah Boswell
from Laboratory of Systems Pharmacology, Single Cell Core; and Dr
Radhika Khetani ,Training Director at the Harvard Chan Bioinformatics
Core RNA-seq workshop available at
https://hbctraining.github.io/rnaseq-cb321/.*
