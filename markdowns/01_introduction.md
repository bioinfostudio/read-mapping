## Introduction

The goal of this session is to perform an unspliced alignment
for a small subset of raw reads. We will align raw sequencing data to
the mouse genome using Bowtie2 and then we will manipulate the SAM output in order to visualize the alignment on the IGV browser.

### Prepare the Environment

We will use one data set in this practical, which can be found in the
`mapping` directory.

First, go to the right folder, where the data are stored.

```bash
$ cd /home/bioinfo/mapping
```

The `.fastq` file that we will align is called `Oct4.fastq`. This file
is based on Oct4 ChIP-seq data published by [Chen *et al.* (2008)](https://www.cell.com/fulltext/S0092-86740800617-X). For
the sake of time, we will align these reads to a single mouse chromosome.
