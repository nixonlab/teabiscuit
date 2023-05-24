# Telescope

## Alignment

### SAM/BAM files

The [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) is a standardized file format for representing sequence alignments to a reference genome. It is able to include all the data needed to describe an alignment and the original input FASTQ files can be reconstructed (lossless). SAM files are tab-delimited text files and can be viewed with common command-line tools. BAM files are compressed (zipped) binary files that may be sorted and indexed for fast retrieval of alignments overlapping a genomic region.

The [samtools](http://www.htslib.org/) package is useful for working with SAM/BAM format files.

SAM/BAM files can also be used to store unaligned data (uBAMs)

CRAM files are able to achieve better compression than BAM files by storing reads aligned to a reference sequence and only storing the bases that differ from the reference. One may also opt for controlled loss of data for even better compression, i.e. eliminating read names or qualities.


#### Viewing BAM files

Activate teabiscuit environment:

```bash
conda activate teabiscuit
```

```bash
samtools view ENCSR693KOP/Aligned.out.bam | less -S
```

This python code will print out alignments for all the reads with 500 alignments:

```python
import pysam
af = pysam.AlignmentFile('ENCSR693KOP/Aligned.out.bam')
iter = af.fetch(until_eof=True)
for aseg in iter:
    if aseg.has_tag('NH') and aseg.get_tag('NH')==500:
        print(aseg.to_string())
```

## Visualizing alignments

Tools for visualizing alignments:

+ [Intergrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)
+ [UCSC Genome Browser](https://genome.ucsc.edu/)

### Download and install IGV

Go to [https://software.broadinstitute.org/software/igv/download](https://software.broadinstitute.org/software/igv/download)

Follow the directions for your OS. I recommend the "Java included" versions unless you have a good reason not to.


PTPRC (CD45, microglia marker)
