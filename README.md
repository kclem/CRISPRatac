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


## Parameters:
```
CRISPRatac: tools for aligning ATAC-seq reads to CRISPR-cleaved genomes

positional arguments:
  settings_file         Tab-separated settings file (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --debug               Tab-separated settings file (default: False)
  --root ROOT           Output directory file root (default: None)
  --keep_intermediate   If true, intermediate files are not deleted (default:
                        False)
  --cuts [CUTS ...], --cut_sites [CUTS ...]
                        Cut sites in the form chr1:234 (multiple cuts are
                        separated by spaces) (default: [])
  --genome GENOME       Genome sequence file for alignment. This should point
                        to a file ending in ".fa", and the accompanying index
                        file (".fai") should exist. (default: None)
  --bowtie2_genome BOWTIE2_GENOME
                        Bowtie2-indexed genome file. (default: None)
  --fastq_r1 FASTQ_R1   Input fastq r1 file (default: None)
  --fastq_r2 FASTQ_R2   Input fastq r2 file (default: None)

Custom target settings:
  --alignment_extension ALIGNMENT_EXTENSION
                        Number of bp to extend beyond av read length around
                        cut site for custom index (default: 50)

Alignment cutoff parameters:
  --arm_min_seen_bases ARM_MIN_SEEN_BASES
                        Number of bases that are required to be seen on each
                        "side" of translocated reads. E.g. if a artificial
                        target represents a translocation between chr1 and
                        chr2, arm_min_seen_bases would have to be seen on chr1
                        as well as on chr2 for the read to be counted.
                        (default: 15)
  --arm_min_matched_start_bases ARM_MIN_MATCHED_START_BASES
                        Number of bases that are required to be matching (no
                        indels or mismatches) at the beginning of the read on
                        each "side" of artifical targets. E.g. if a artificial
                        target represents a translocation between chr1 and
                        chr2, the first arm_min_matched_start_bases of the
                        read would have to match exactly to chr1 and the last
                        arm_min_matched_start_bases of the read would have to
                        match exactly to chr2 (default: 10)
  --ignore_n IGNORE_N   If set, "N" bases will be ignored. By default (False)
                        N bases will count as mismatches in the number of
                        bases required to match at each arm/side of the read
                        (default: False)

CRISPResso settings:
  --crispresso_min_count CRISPRESSO_MIN_COUNT
                        Min number of reads required to be seen at a site for
                        it to be analyzed by CRISPResso (default: 50)
  --crispresso_min_aln_score CRISPRESSO_MIN_ALN_SCORE
                        Min alignment score to reference sequence for
                        quantification by CRISPResso (default: 20)
  --crispresso_quant_window_size CRISPRESSO_QUANT_WINDOW_SIZE
                        Number of bp on each side of a cut to consider for
                        edits (default: 1)

Pipeline parameters:
  --samtools_command SAMTOOLS_COMMAND
                        Command to run samtools (default: samtools)
  --bowtie2_command BOWTIE2_COMMAND
                        Command to run bowtie2 (default: bowtie2)
  --crispresso_command CRISPRESSO_COMMAND
                        Command to run CRISPResso2 (default: CRISPResso)
  --n_processes N_PROCESSES
                        Number of processes to run on (may be set to "max")
                        (default: 1)
```
