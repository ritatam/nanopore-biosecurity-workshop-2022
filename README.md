# Nanopore Sequencing for Biosecurity Workshop 2022 

### High-accuracy long-read amplicon consensus sequences reconstruction with Nanopore sequencing data

This is a tutorial as part of bioinformatics component for the course [Nanopore Sequencing for Biosecurity](https://cba.anu.edu.au/news-events/events/nanopore-sequencing-biosecurity) given on 28-30 November 2022 at the Australian National University.

## Introduction
*Amplicon sequencing* is a highly targeted next-generation sequencing approach that enables analyses of genetic variations in specific genomic regions. Short-read sequencing has been commonly used in amplicon-related research for its high-throughput nature and low error rate, yet its read length limits the maximum amplicon size. Long-read amplicon sequencing allows for capturing long-range genetic information, which offers better solutions for in-depth analyses, such as resolving longer structural variations and phasing.

This tutorial will guide you through bioinformatic workflows to reconstruct high-accuracy amplicon consensus sequences with Nanopore sequencing data, via the two case studies introduced in the presentation.

## **Case study 1 - Fungal identification by generating consensus sequence of full-ITS region for a species**

### Fungal barcoding with ITS region
Fungal identification using DNA barcoding based on internal transcribed spacer (ITS) located among nuclear rRNA genes has been routinely used by a lot of diversity studies. Full-ITS region contains two hypervariable spacers, ITS-1 and ITS-2, separated by the highly conserved 5.8S rRNA gene. The variations in ITS-1 and ITS-2 are species-specific, enabling the discriminatory power for taxonomy classification.

![fungal ITS](./fig/fuganl_ITS.png)

**Here, we will reconstruct a consensus sequence for the full-ITS region of a fungal species using long-read amplicons, which we will use to identify the fungus later.** My colleagues and I amplified the ITS region from the fungal gDNA, sequenced and basecalled the ampilcon long reads, then filtered the reads by read quality (Q>15) and length (2.5-3.5kbp). We will go through the following bioinformatic workflow for consensus reconstruction using these reads.

'insert workflow'

### 1. Generating a draft sequence

We will use the filtered reads in a fastq file to generate a draft sequence. This will be done using a clustering algorithm called USEARCH, which clusters reads based on a user-defined sequence similarity threshold (e.g. 0.75). The centroid (representative) read sequence of the most populated cluster will be used as the draft sequence. This way, we can remove potential contamination or sequencing artifacts, and make sure that only similar reads from the amplicon under question will be considered in downstream analyses.

### Define USEARCH file path

Since USEARCH is installed as one executable file, we will first define a variable called <code>USEARCH</code> as its file path. 

    course_user> USEARCH=/usearch11.0.667_i86linux32

Now, run the USEARCH program by simplying adding a $ sign to the variable. 

    course_user> $USEARCH
    
... which will print the following message. It means USEARCH is all set!

    usearch v11.0.667_i86linux32, 4.0Gb RAM (132Gb total), 32 cores
    (C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
    https://drive5.com/usearch
    License: personal use only

We can now proceed to clustering our reads!

### Run USEARCH to cluster the reads

First, change into the <code>fungal_ITS</code> directory, where the input <code>reads.fastq</code> file is located. You can use command <code>less</code> to view its content in a full-screen display mode. 

    course_user> cd ~/amplicon_prac/fungal_ITS
    course_user> less reads.fastq
    
To quit the display mode, simply press the "Q" key. 

Create an output directory called <code>draft</code> for USEARCH.

    course_user> mkdir draft

Now let's run a USEARCH command to cluster our reads. (Remember to use $ to call variable)

    $USEARCH -cluster_fast reads.fastq -id 0.75 -strand both -centroids draft/centroids.id75.fasta -sizeout

Option explanations from USEARCH documentation:
- <code>-cluster_fast</code>: Cluster sequences using UCLUST
- <code>-id 0.75</code>: Minimum sequence identity threshold. Identity is defined as the fraction of columns in an alignment with matching letters.
- <code>-strand both</code>: Search for hits of both forward and reverse complement strands.
- <code>-centroids</code>: write the centroid sequences to a FASTA file.
- <code>-sizeout</code>: annotate cluster size in the header of every centroid sequence in the FASTA file.

This will generate a fasta file storing the centroid/representative sequences from all clusters. They are sorted in descending order of cluster size, which means the first sequence on top represents the highest number of similar reads. You can inspect the fasta file using <code>less</code>. Note the cluster size annotated at the end of each header. 

    course_user> cd draft
    course_user> less centroids.id75.fasta

 *Tip: You can scroll up and down in the <code>less</code> display mode to look at other centroid sequences. (What does size=1 indicate? How do their centroid sequences look like?)*

Let's extract the first centroid sequence in the fasta file as the draft sequence.

    awk "/^>/ {n++} n>1 {exit} {print}" centroids.id75.fasta > draft.fasta

Note: the draft sequence is a Nanopore **read** that is highly representative, so it might still contain a lot of errors that need to be polished up downstream!

### Generate the consensus sequence and polish it

Now we can generate consensus sequence from the draft. The centroid sequence 