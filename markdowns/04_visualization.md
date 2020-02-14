## Visualize alignments in IGV

IGV is a genome browser. Please check their [website](https://software.broadinstitute.org/software/igv/) for all the formats that IGV can
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
tool `genomeCoverageBed` from [BEDTools](https://bedtools.readthedocs.io/en/latest/). Then we convert the bedgraph
into a bigWig binary file called `Oct4.bw`, using `bedGraphToBigWig`
from [the UCSC tools](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads):

```bash
mapping$ genomeCoverageBed -bg -ibam Oct4.sorted.bam > Oct4.bedgraph
mapping$ bedGraphToBigWig Oct4.bedgraph bowtie_index/mouse.mm10.genome Oct4.bw
```

The commands above take as input a file called
`mouse.mm10.genome` that is stored under the subdirectory
`bowtie_index`. These genome files are tab-delimited and describe the
size of the chromosomes for the organism of interest. When using the
[UCSC Genome Browser](https://genome.ucsc.edu/index.html), [Ensembl](https://www.ensembl.org/), or [Galaxy](https://usegalaxy.org/), you typically indicate which
species/genome build you are working with. The way you do this for
BEDTools is to create a “genome” file, which simply lists the names of
the chromosomes (or scaffolds, etc.) and their size (in basepairs).

BEDTools includes pre-defined genome files for human and mouse in the
`genomes` subdirectory included in the BEDTools distribution.

In order to see the aligned reads of the BAM file, you need to zoom in
to a specific region. For example, look for gene `Lemd1` in the search
box.

<div id="igv-div"></div>

### The main difference between the visualization of BAM and bigWig files

The actual alignment of reads that stack to a particular region can be displayed using the information stored in a BAM format. The bigWig format is for display of dense, continuous data that will be displayed in the Genome Browser as a graph.

<script type="text/javascript">
  var igvDiv = document.getElementById("igv-div");
  var options =
    {
        genome: "mm10",
        locus: "Chr1",
        tracks: [
            {
                type: "alignment",
                format: "bam",
                name: "BAM",
                url: "gs://bioinfostudio/mapping/data/Oct4.sorted.bam",
                indexURL: "gs://bioinfostudio/mapping/data/Oct4.sorted.bam.bai",
            },
            {
                type: "wig",
                format: "bigwig",
                name: "BigWig",
                url: "gs://bioinfostudio/mapping/data/Oct4.bw",
                indexed: false,
            },
        ]
    };
    igv.createBrowser(igvDiv, options)
</script>
