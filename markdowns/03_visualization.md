## Manipulate SAM output

SAM files are rather big and when dealing with a high volume of NGS
data, storage space can become an issue. As we have already seen, we can
convert SAM to BAM files (their binary equivalent that are not human
readable) that occupy much less space.

Convert SAM to BAM using `samtools view` and store the output in the
file `Oct4.bam`. You have to instruct `samtools view` that the input is
in SAM format (`-S`), the output should be in BAM format (`-b`) and that
you want the output to be stored in the file specified by the `-o`
option:

```bash
mapping$ samtools view -bSo Oct4.bam Oct4.sam
```

Compute summary stats for the Flag values associated with the alignments
using:

```bash
mapping$ samtools flagstat Oct4.bam
```

## Visualize alignments in IGV

IGV is a stand-alone genome browser. Please check their [website](<http://www.broadinstitute.org/igv/>) for all the formats that IGV can
display. For our visualization purposes we will use the BAM and bigWig formats.

When uploading a BAM file into the genome browser, the browser will look
for the index of the BAM file in the same folder where the BAM files is.
The index file should have the same name as the BAM file and the suffix
`.bai`. Finally, to create the index of a BAM file you need to make sure
that the file is sorted according to chromosomal coordinates.

Sort alignments according to chromosomal position and store the result
in the file with the prefix `Oct4.sorted`:

```bash
mapping$ samtools sort Oct4.bam -o Oct4.sorted.bam
```

Index the sorted file.

```bash
mapping$ samtools index Oct4.sorted.bam
```

The indexing will create a file called `Oct4.sorted.bam.bai`. Note that
you don’t have to specify the name of the index file when running
`samtools index`, it simply appends a `.bai` suffix to the input BAM
file.

Another way to visualize the alignments is to convert the BAM file into
a bigWig file. The bigWig format is for display of dense, continuous
data and the data will be displayed as a graph. The resulting bigWig
files are in an indexed binary format.

The BAM to bigWig conversion takes place in two steps. Firstly, we
convert the BAM file into a bedgraph, called `Oct4.bedgraph`, using the
tool `genomeCoverageBed` from BEDTools. Then we convert the bedgraph
into a bigWig binary file called `Oct4.bw`, using `bedGraphToBigWig`
from the UCSC tools:

```bash
mapping$ genomeCoverageBed -bg -ibam Oct4.sorted.bam -g bowtie_index/mouse.mm10.genome > Oct4.bedgraph
mapping$ bedGraphToBigWig Oct4.bedgraph bowtie_index/mouse.mm10.genome Oct4.bw
```

Both of the commands above take as input a file called
`mouse.mm10.genome` that is stored under the subdirectory
`bowtie_index`. These genome files are tab-delimited and describe the
size of the chromosomes for the organism of interest. When using the
UCSC Genome Browser, Ensembl, or Galaxy, you typically indicate which
species/genome build you are working with. The way you do this for
BEDTools is to create a “genome” file, which simply lists the names of
the chromosomes (or scaffolds, etc.) and their size (in basepairs).

BEDTools includes pre-defined genome files for human and mouse in the
`genomes` subdirectory included in the BEDTools distribution.

Now we will load the data into the IGV browser for visualization. In
order to launch IGV double click on the `IGV 2.3` icon on your Desktop.
Ignore any warnings and when it opens you have to load the genome of
interest.

On the top left of your screen choose from the drop down menu
`Mouse (mm10)`. If it doesn’t appear in list, click `More ..`, type
`mm10` in the Filter section, choose the mouse genome and press OK. Then
in order to load the desire files go to:

    File > Load from File

On the pop up window navigate to Desktop -> chipseq folder and select
the file `Oct4.sorted.bam`.

Repeat these steps in order to load `Oct4.bw` as well.

Select `chr1` from the drop down menu on the top left. Right click on
the name of `Oct4.bw` and choose Maximum under the Windowing Function.
Right click again and select Autoscale.

In order to see the aligned reads of the BAM file, you need to zoom in
to a specific region. For example, look for gene `Lemd1` in the search
box.

### The main difference between the visualization of BAM and bigWig files

The actual alignment of reads that stack to a particular region can be displayed using the information stored in a BAM format. The bigWig format is for display of dense, continuous data that will be displayed in the Genome Browser as a graph.
