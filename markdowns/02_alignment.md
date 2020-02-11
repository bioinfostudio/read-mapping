## Alignment

You already know that there are a number of competing tools for short
read alignment, each with its own set of strengths, weaknesses, and
caveats. Here we will try Bowtie2, a widely used ultrafast, memory
efficient short read aligner.

Bowtie2 has a number of parameters in order to perform the alignment. To view them all type

```bash
$ bowtie2 --help
```

Bowtie2 uses indexed genome for the alignment in order to keep its
memory footprint small. Because of time constraints we will build the
index only for one chromosome of the mouse genome. For this we need the
chromosome sequence in FASTA format. This is stored in a file named
`mm10`, under the subdirectory `bowtie_index`.

The indexed chromosome is generated using the command:
DO NOT run this command. This has already been run for you.

**bowtie2-build bowtie_index/mm10.fa bowtie_index/mm10**

This command will output 6 files that constitute the index. These files
that have the prefix `mm10` are stored in the `bowtie_index`
subdirectory. To view if they files have been successfully created type:

```bash
mapping$ ls -l bowtie_index
```

Now that the genome is indexed we can move on to the actual alignment.
The first argument for `bowtie2` is the basename of the index for the
genome to be searched; in our case this is `mm10`. We also want to make
sure that the output is in SAM format using the `-S` parameter. The last
argument is the name of the FASTQ file.

Align the Oct4 reads using Bowtie2:

```bash
mapping$ bowtie2 -x bowtie_index/mm10 -q Oct4.fastq > Oct4.sam
```

The above command outputs the alignment in SAM format and stores them in
the file `Oct4.sam`.

In general before you run Bowtie2, you have to know what quality
encoding your FASTQ files are in. The available FASTQ encodings for
bowtie are:

- phred33-quals: Input qualities are Phred+33 (default).
- phred64-quals: Input qualities are Phred+64 (same as `–solexa1.3-quals`).
- solexa-quals: Input qualities are from GA Pipeline ver. < 1.3.
- solexa1.3-quals: Input qualities are from GA Pipeline ver. >= 1.3.
- integer-quals: Qualities are given as space-separated integers (not ASCII).
    
The FASTQ files we are working with are Sanger encoded (Phred+33), which
is the default for Bowtie2.

Bowtie2 will take 2-3 minutes to align the file. This is fast compared
to other aligners which sacrifice some speed to obtain higher
sensitivity.

Look at the top 10 lines of the SAM file using head (record lines are
wrapped). Then try the second command, note use arrow navigation and to
exit type ‘q’.

```bash
$ head Oct4.sam
$ less -S Oct4.sam
```

**Can you distinguish between the header of the SAM format and the actual alignments?**

The header line starts with the letter ‘@’, i.e.:

  ----- ------------ -------------- ---------- ----------------------------------------------------------------------------------------------------------
  @HD   VN:1.0       SO:unsorted               
  @SQ   SN:chr1      LN:195471971              
  @PG   ID:Bowtie2   PN:bowtie2     VN:2.2.4   CL:“/tools/bowtie2/bowtie2-default/bowtie2-align-s –wrapper basic-0 -x bowtie\_index/mm10 -q Oct4.fastq”
  ----- ------------ -------------- ---------- ----------------------------------------------------------------------------------------------------------

While, the actual alignments start with read id, i.e.:

  -------------- ---- ------ -----
  SRR002012.45   0    etc    
  SRR002012.48   16   chr1   etc
  -------------- ---- ------ -----

**What kind of information does the header provide?**

- @HD: Header line; VN: Format version; SO: the sort order of alignments.
- @SQ: Reference sequence information; SN: reference sequence name; LN: reference sequence length.
- @PG: Read group information; ID: Read group identifier; VN: Program version; CL: the command line that produces the alignment.
