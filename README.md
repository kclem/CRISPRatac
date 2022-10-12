# CRISPRatac

CRISPRatac aligns ATAC-sequencing reads to a genome and to regions produced by combination of CRISPR cut sites

Program requirements:
samtools
bowtie2
matplotlib
CRISPResso2

Parameters may be passed in via the command line, or a settings file.

To run, create the file 'settings.txt'

\--- settings.txt ---
```
fastq_r1	read1.fq.gz
fastq_r2	read2.fq.gz
genome	/data/pinello/COMMON_DATA/REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa
cuts	chr1:1111111 chr2:222222
```

Then run 'python CRISPRatac.py settings.txt'
