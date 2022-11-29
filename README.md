# Nanopore Sequencing for Biosecurity Workshop 2022 

This tutorial will take you through steps of **high-accuracy long-read amplicon consensus sequences reconstruction with Nanopore sequencing data**.

## Introduction
*Amplicon sequencing* is a highly targeted next-generation sequencing approach that enables analyses of genetic variations in specific genomic regions. Short-read amplicon sequencing has been commonly adopted for its high-throughput nature and low error rate, yet its read length limits the maximum amplicon size. Long-read amplicon sequencing allows for capturing long-range genetic information, which offers better solutions for in-depth analyses such as resolving structural variations and phasing.

This tutorial will guide you through a simple workflow to reconstruct high-accuracy amplicon consensus sequences with Nanopore long-read sequencing data, via the two introduced case studies.

## **Case study 1 - Fungal identification by generating consensus sequence of full-ITS region for a species**


#### Fungal barcoding with ITS region
Fungal identification using DNA barcoding based on internal transcribed spacer (ITS) regions has been routinely used by a lot of diversity studies. ITS region is located among nuclear rRNA genes. which contains two hypervariable spacers, ITS-1 and ITS-2, separated by the highly conserved 5.8S rRNA gene. The variations in ITS-1 and ITS-2 are species-specific, thus enabling the discriminatory power for taxonomy classification.

'insert ITS image'

Here, we will reconstruct a consensus sequence for the full-ITS region of a fungal species using long-read amplicons, which we will later use to identify the fungus. We have amplified the ITS region from fungal gDNA, sequenced and basecalled the ampilcon long reads, and filtered by read quality (Q>15) and length (2.5-3.5kbp). 

'insert workflow'


We will use the filtered reads in a fastq file to generate a draft sequence. This will be done using a clustering algorithm called USEARCH, which clusters reads based on a user-defined sequence similarity threshold (e.g. 0.75). The centroid (representative) read sequence of the most populated cluster will be used as the draft. This way, we can remove potential contamination or sequencing artifacts, and make sure that only similar reads from amplicon under question will be considered downstream.

#### Define USEARCH file path

Since USEARCH is installed as one executable file, we will first define a variable called <code>USEARCH</code> as its file path. 

    course_user> USEARCH=/usearch11.0.667_i86linux32

Now, call USEARCH by simplying adding a $ sign to the variable. 

    course_user> $USEARCH
    
... which will print the following message. It means USEARCH is all set!

    usearch v11.0.667_i86linux32, 4.0Gb RAM (132Gb total), 32 cores
    (C) Copyright 2013-18 Robert C. Edgar, all rights reserved.
    https://drive5.com/usearch
    License: personal use only

We can now proceed to clustering our reads!

#### Run USEARCH to cluster the reads

First, change into the <code>fungal_ITS</code> directory. The input read fastq file is located in the *fungal_ITS* directory. You can use command <code>less</code> to view its content in a full-screen display mode. 

    course_user> cd ~/amplicon_prac/fungal_ITS
    course_user> less reads.fastq
    
To quit the display mode, simply press the "Q" key. 

Create a output directory called <code>draft</code> for USEARCH.

    course_user> mkdir draft

Now let's run a USEARCH command to cluster our reads. (Remember to add the $ sign to call the variable!)

    $USEARCH -cluster_fast reads.fastq -id 0.75 -strand both -centroids draft/centroids.id75.fasta -sizeout

Option definitions from USEARCH documentation:
<code>-cluster_fast</code>: Cluster sequences using UCLUST
<code>-id 0.75</code>: Minimum sequence identity threshold. Identity is the fraction of columns in an alignment with matching letters.
<code>-strand both</code>: Search for hits of both forward and reverse complement strands.
<code>-centroids</code>: write the centroid sequences to a FASTA file.
<code>-sizeout</code>: annotate cluster size in the header of every centroid sequence in the FASTA file.
