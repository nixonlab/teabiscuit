# Reference Genome Practical

In your teabiscuit directory, there are two reference "human" genomes.

1. Are the genomes the same? How are they different?
2. Where did these reference genomes come from?
3. Which genome would you use and why?

**HINT**: The `*.fai` files alongside each genome are plain-text index files used by `samtools faidx` and contain information about each sequence in the genome FASTA file. See [faidx manual page](http://www.htslib.org/doc/faidx.html)

#### Some useful commands:

##### `head`

Unzip gzipped file and view the first 10 lines:

```bash
zcat GRCh38.d1.vd1.fa.gz | head
```

##### `less`

View gzipped file contents one page at a time:

```bash
zcat GRCh38.d1.vd1.fa.gz | less
```

+ Press up and down arrow keys to move up and down
+ Press "F" to advance a whole screen
+ Press "Q" to quit

##### `grep`

Print all the lines beginning with ">"

```bash
zcat GRCh38.d1.vd1.fa.gz | grep '^>'
```



