# Alignment with STAR

See the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
for details about running STAR.

## Index genome

Example index command

```
STAR\
 --runThreadN 16\
 --runMode genomeGenerate\
 --genomeDir data/STAR_gdc38_gencode38\
 --outFileNamePrefix data/STAR_gdc38_gencode38\
 --genomeFastaFiles data/refs/GRCh38.d1.vd1.fa.gz\
 --sjdbGTFfile data/refs/gencode.v38/gencode.v38.REF.annotation.gtf.gz\
 --sjdbOverhang 74
```

## Run aligner

```
cd teabiscuit/01-rnaseq_workflow
conda activate teabiscuit
STAR\
 --runThreadN 1\
 --genomeDir STAR_gdc38_gencode38\
 --readFilesIn SRR19243439_1.short.fastq SRR19243439_2.short.fastq\
 --outSAMattributes NH HI NM MD AS XS\
 --outSAMtype BAM Unsorted\
 --quantMode GeneCounts\
 --outSAMstrandField intronMotif\
 --outFilterMultimapNmax 500\
 --outFilterMultimapScoreRange 5\
 --outSAMunmapped Within KeepPairs
```

### Important multimapping parameters

+ **outFilterMultimapNmax** max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
+ **outFilterMultimapScoreRange** the score range below the maximum score for multimapping alignments

## Output

1. What are the output files produced by STAR?
2. 