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

