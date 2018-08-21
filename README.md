# insiM: *in silico* Mutator software for bioinformatics pipeline validation of clinical Next Generation Sequencing (NGS) assays

## insiM_v1.0

### Contact

Sushant Patil: spatil6@bsd.uchicago.edu

Sabah Kadri: skadri@bsd.uchicago.edu

### Publication (please cite if using software):

Currently in press at Journal of Molecular Diagnostics. 


### Introduction 

Lack of reliable reference samples containing different mutations of interest across large sets of disease-relevant loci limits the extensive validation clinical next generation sequencing (NGS) assays and their associated bioinformatics pipelines. Here, we have created a highly flexible tool, insiM (in silico Mutator) to introduce point mutations, insertions, deletions, and duplications of any size into real datasets of amplicon-based or hybrid-capture NGS assay. insiM accepts an alignment file along with target territory and produces paired-end FASTQ files containing specified mutations via modification of original sequencing reads.  Mutant signal is thus created within the context of existing real-world data to most closely mimic assay performance.  Resulting files may then be passed through the assay’s bioinformatics pipeline to assist with assay/bioinformatics validation and to identify performance gaps in detection. 

### How it works

The software extracts uniquely mapped reads at each specified mutation locus randomly, based on the specified variant frequency, and outputs all (mutated or otherwise) reads to new paired FASTQ files. insiM also outputs the exact fraction of mutated reads and variant sequence (for insertion and SNV) at each mutation locus to a VCF file for downstream performance evaluation of the bioinformatics pipeline. 
insiM works by evaluating reads in the context of their pairs and mutating one or both reads based on the overlap at the mutation position (highlighted in blue and green in figure, respectively). For example, if both read mates of a pair overlap the mutation locus, then both will be altered. The new FASTQ files should be nearly identical to the original FASTQ files other than the mutated reads and any hard clipping that the original alignment produced in the input BAM. 

The specific logic of mutation generation in NGS data varies depending on the type of anomaly as described below.

(i) For SNVs, only the specified base is mutated within the read. 

(ii) For insertions and duplications, new sequence is added at the specified locus and the end of the existing reads are then truncated to the original read lengths. 

(iii) In case of deletion simulations in amplicon-based assays, the reads are extended to the specified read length by first appending the genomic sequence until the amplicon end and then non-specific nucleotides may be added. In future updates, library-specific primers may also be used to more accurately extend the reads, past the amplicon end.  

(iv) For deletions in hybrid-capture assays, insiM first creates a library fragment from the read pair, introduces the deletion in the fragment and then re-constructs the read pair. In cases where the fragment is longer than the sum of the read lengths, the read mates do not overlap each other. For such non-overlapping reads, the region of the fragment that is not covered by reads is derived from the reference genome, whereas for overlapping reads (shorter fragments), the fragment sequence in the overlapped region is obtained from the forward read. The deletion is introduced in the fragment at the specified position, and downstream genomic sequence is used to extend the fragment by the deletion length, thus preserving the fragment length/insert size. Finally, read sequences of the original length are extracted from either end of mutated fragment. 

![figure1-1](https://user-images.githubusercontent.com/9405995/44405375-b3dce400-a51e-11e8-960d-1fcf3216c7e0.jpg)

### Installation and usage

insiM has been tested with python v2.7.6. It requires ‘pysam’ module in the python PATH. 'Input data' section below lists all the mandatory and optional parameters required by the software in the form of a plain text configuration file. insiM can be executed as:

``` python insiM_v1.0.py <configuration file> ```

### Input data

insiM input parameters in configuration file. 

| Argument		| Value           |
| ------------------- |---------------|
| -assay (m)	| Assay type: 'amplicon' or 'capture' (Hybrid Capture) |
| -bam (m)      | BAM file     |
| -target (m) | BED file containing target loci. Mutations are introduced at the mid-points of the BED coordinates. For simulating mutation at a specific position, start and end BED coordinates should denote that exact same position.       |
| -mutation (m)		| Mutation type: 'snv' (Single Nucleotide Variants), 'ins' (Insertion), 'del' (Deletion), 'dup' (Duplication) or 'mix' (Multi-type) |
| -genome (m)		| Genome FASTA file |
| -ampliconbed (m)		| Amplicon BED file. Required only when -assay == 'amplicon’ |
| -read (m)	| Read length. Required only when -assay == 'amplicon’ |
| -vaf (m)	| Variant allele frequency. For multi-type mutations, a separate value can be specified for each locus in BED file.  |
| -len (o)	| Size of insert. If not specified, random sequence of length 10 is used for ins/del/dup |
| -seq (o)	| Sequence of insert. For insertion, this sequence is inserted; if not specified, random sequence of length specified by -len is inserted |
| -out (o)	| Output FASTQ basename |

m = mandatory, o = optional

### FAQs
