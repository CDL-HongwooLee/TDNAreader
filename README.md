# T-DNAreader: Fast and precise identification of T-DNA insertion sites in plant genomes using RNA-sequencing data
  
![Pipeline](https://github.com/user-attachments/assets/03d6f6aa-8fde-4780-9d17-79505598d499)

For any question about TDNAreader, please contact miso5103@snu.ac.kr

## Dependencies
python3 >= 3.8

Bowtie2 >= 2.5.1 (https://github.com/BenLangmead/bowtie2)

Samtools >= 1.17 (http://www.htslib.org)

Trim_galore >= 0.6.10 (https://github.com/FelixKrueger/TrimGalore)

Deeptools >= 3.5.1 (https://deeptools.readthedocs.io/en/develop/content/installation.html)

GNU parallel (https://www.gnu.org/software/parallel/)

Required packages: sys, os, subprocess, argparse, time, multiprocessing, functools, pandas, numpy

## Installation
Github install

```
git clone https://github.com/CDL-HongwooLee/TDNAreader
```

# Quick start

- TDNAreader requires pre-aligned BAM (SAM) files and TDNA.fa file as inputs.
- The BAM (SAM) files generated by STAR, Bowtie2 (with --local), BWA or other local sequencing aligners are acceptable.
- **Note: For paired-end sequencing, forward and reverse FASTQ files are separately aligned to the plant genome (treating them as single-end sequencing)
- Bowtie2 index for the provided TDNA.fa is required. (can be generated by ```bowtie2-build```)



### Quick run with the default parameters

For single-end sequencing,
```
python3 TDNAreader.py -b sample_rep1.bam,sample_rep2.bam,sample_rep3.bam -tdna tdna.fa -o ./ -p sample_name -t 3
```

For paired-end sequencing,
```
python3 TDNAreader.py -b sample_rep1_1.bam:sample_rep1_2.bam,sample_rep2_1.bam:sample_rep2_2.bam,sample_rep3_1.bam:sample_rep3_2.bam --paired -tdna tdna.fa -o ./ -p sample_name -t 6
```

For example run
```
python3 TDNAreader.py -b BAM/SRR13765633_1Aligned.sortedByCoord.out.bam:BAM/SRR13765633_2Aligned.sortedByCoord.out.bam --paired -tdna pROK2.fas -o Results/ -p test -t 2
```

## TDNAreader.py

### Usage

```
usage: TDNAreader.py [-h] -b BAM [--paired] [-tdna TDNA] [-o OUTDIR] -p PREFIX [-bL BOWTIEL] [-d D] [-l1 L1] [-l2 L2] [-m M] [-t THREAD] [-bed BED] [-bl BLACKLIST] [--tmp]

optional arguments:
  -h, --help       show this help message and exit
  -b BAM         Input BAM or SAM file(s).
                 Multiple files should be separated by ',' (e.g., rep1.bam,rep2.bam,rep3.bam).
                 For paired-end reads, pairs of BAM files should be separated by ':' (e.g., rep1_1.bam:rep1_2.bam,rep2_1.bam:rep2_2.bam).
                 Strongly recommend using all biological replicates.

  --paired       The input BAM files are from paired-end sequencing.
                 **Note: it is required to align paired-end reads separately to the plant genome (treating them as single-end reads).
                 Paired BAM files should be delimited by ':' in the input.

  -tdna TDNA     Path to the T-DNA fasta file containing sequences from LB to RB.
                 A Bowtie2 index for this fasta file is required and must be created prior to running T-DNAreader.
                 Bowtie2 index can be generated using 'bowtie2-build tdna.fa'
                 Sequences outside the LB-RB regions (e.g., vector backbone sequences) should be excluded.

  -o OUTDIR      Path to the output directory where all output files will be written.
                 Default = './'

  -p PREFIX      Prefix for naming output files.

  -bL BOWTIEL    Length of seed substring (in bp) used in Bowtie2 alignment.
                 Must be >3 and <=l1.
                 Default = 10

  -d D           Maximum distance from T-DNA borders to distinguish between border regions and intergnal regions.
                 Reads mapped within this distance from T-DNA borders used the threshold specified by -l1.
                 Reads mapped beyond this distance use the threshold specified by l2.
                 Default = 200 (100<=d<=300 is recommended).

  -l1 L1         Minimum length of sequences matched to T-DNA junctions (reads mapped close to T-DNA borders).
                 Higher values increase stringecy, potentially reducing false positives but may miss some true TISs.
                 Default = 18 (15<=l1<=300 is recommended).

  -l2 L2         Minimum length of sequences matched with T-DNA internal regions (reads mapped far from T-DNA borders).
                 Higher values increase stringecy, potentially reducing false positives but may miss some true TISs.
                 Default = 30 (30<=l2<=50 is recommended)

  -m M           Maximum number of mismatches allowed in T-DNA alignment.
                 Default = 1 (bp). (0<=m<=3 is recommended)

  -t THREAD      Number of threads to use.
                 Default = 1. (t = {number of samples} is recommended more most efficient analysis)

  -bed BED       (Optional) A bed file containing genomic information.
                 BED6 format is required (6 columns).
                 Gene names or identifiers should be in column 4.
                 If specified, overlapping genes with TISs will be identified.

  -bl BLACKLIST  (Optional) A bed file specifying regions to exclude from analysis.
                 If the inserted T-DNA of transgenic plants contain endogenous genomic sequences,
                 it is recommended to specify the positions of those sequences to avoid misinterpretation.

  --tmp          (Optional) Keep intermediate files generated by T-DNAreader

```



## Input and output data

### Input

#### (-b) .bam(.sam) file 

- The BAM (SAM) files generated by STAR, Bowtie2 (with --local), BWA or other local sequencing aligners are acceptable.
  - For STAR, use {name}.sortedByCoord.bam files.
- **Note: For paired-end sequencing, forward and reverse FASTQ files are separately aligned to the plant genome (treating them as single-end sequencing).
- (Optional) The modification of STAR alignment options may improve the sensitivity of TDNAreader. You can specify ```--outFilterMatchNminOverLread``` and ```--outFilterScoreMinOverLread``` to an lower value (default=0.66) to identify more reads. The range from 0.33 to 0.5 is recommended depending on your read length.

#### tdna.fa

- FASTA file containing T-DNA sequences inserted into the plant genome: from left-border (LB) to right-border (RB)
- Sequences outside LB-RB regions should be excluded
- Bowtie2 index for this file is required prior to running TDNAreader (use ```bowtie2-build```)

### Output

Output file is generated as ```{outdir}/{prefix}.TDNA.txt```

Example 
```
Sample  Position(REF)   Geneid  Readid  Flag(REF)       Cigar(REF)      Sequence        TDNA    Flag(TDNA)      Position(TDNA)  Cigar(TDNA)
1_floe1-1_2-R   Chr4:14016003-TDNA      AT4G28300       SRR13765633.1338395     0       131M19S TCGTCGTCTCACTCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTA>>> pROK2(SALK)     0       29      19M
1_floe1-1_3-R   Chr4:14016003-TDNA      AT4G28300       SRR13765634.8319149     16      131M19S TCGTCGTCTCACTCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTA>>> pROK2(SALK)     0       29      19M
1_floe1-1_2-R   Chr4:14016003-TDNA      AT4G28300       SRR13765633.11143970    0       130M20S CGTCGTCTCACTCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTAA>>> pROK2(SALK)     0       29      20M
1_floe1-1_2-F   Chr4:14016003-TDNA      AT4G28300       SRR13765633.16534429    16      129M21S GTCGTCTCACTCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTAAT>>> pROK2(SALK)     0       29      21M
1_floe1-1_3-R   Chr4:14016003-TDNA      AT4G28300       SRR13765634.13245767    16      129M21S GTCGTCTCACTCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTAAT>>> pROK2(SALK)     0       29      21M
1_floe1-1_2-F   Chr4:14016003-TDNA      AT4G28300       SRR13765633.17743684    16      119M31S TCTCAACATGGTGAGGACCGTGTCGCCACTCCTGTTCCAGAGCCTAAGAAGAGCGAGAACACCTCTGATGCACACAACCAGCAGCTTGCACTTGCTTTGCCTCACCAAATAGCCCCACA>>>TTGACGCTTAGACAACTTAATAACACATTGC>>> pROK2(SALK)     0       29      31M
```
- Sample: {prefix}_{replicate_n}-{F/R}
- Position(REF): the identified genomic positions of T-DNA insertion (=TISs).
- Geneid: the name of genes containing the identified TISs.
- Readid: the readID of chimeric reads identified by TDNAreader, which contain both genomic and T-DNA sequences.
- Flag(REF): samtools flag of reads in the initial mapping to the plant genome.
- Cigar(REF): CIGAR string of reads in the initial mapping to the plant genome.
- Sequence: sequences of identified reads. The genome:T-DNA junction is indicated by '<<<' or '>>>' depending on the its orientation.
- TDNA: name of TDNA (derived from tdna.fa)
- Flag(TDNA): samtools flag of cleaved reads in the second mapping to the T-DNA sequences.
- Position(TDNA): the aligned position of cleaved reads to the T-DNA sequences.
- Cigar(TDNA): CIGAR string of reads in the second mapping to the T-DNA sequences.


## Example run

Run TDNAreader with example data (available at zenodo, download Examples)

```
python3 TDNAreader.py -b Examples/BAM/SRR13765632_1Aligned.sortedByCoord.out.bam:Examples/BAM/SRR13765632_2Aligned.sortedByCoord.out.bam,Examples/BAM/SRR13765633_1Aligned.sortedByCoord.out.bam:Examples/BAM/SRR13765633_2Aligned.sortedByCoord.out.bam,Examples/BAM/SRR13765634_1Aligned.sortedByCoord.out.bam:Examples/BAM/SRR13765634_2Aligned.sortedByCoord.out.bam -o Examples/Results/ -bed Examples/Reference/Araport11_GTF_genes.Apr2022.rmMC.bed -t 12 -tdna Examples/TDNA/pROK2.fas --paired -p 1_floe1-1
```

