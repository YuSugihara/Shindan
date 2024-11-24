# ViiR User Guide
#### Version 0.0.1

## Table of contents
- [What is ViiR?](#what-is-viir)
- [Installation](#installation)
  + [Dependencies](#dependencies)
  + [Installation](#installation)
- [Usage](#usage)
  + [Example 1 : Run ViiR with default settings](#example-1--run-viir-with-default-settings)
  + [Example 2 : Run ViiR with more threads and CPU memories](#example-2--run-viir-with-more-threads-and-cpu-memories)
  + [Example 3 : Run ViiR with strand specific library](#example-3--run-viir-with-strand-specific-library)


## What is ViiR?

ViiR is a software for 'Virus identification independent of Reference sequence'.

## Installation
### Dependencies
#### Softwares
- [Python3](https://www.python.org/downloads/)
- [wget](https://www.gnu.org/software/wget/)
- [Trinity package](https://github.com/trinityrnaseq/trinityrnaseq)
  + [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)
  + [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  + [Samtools](http://www.htslib.org/doc/samtools.html)
- [HMMER](http://hmmer.org/)
- [R](https://www.r-project.org/)
- [DESeq2](https://bioconductor.org/packages/3.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [RSEM](https://deweylab.github.io/RSEM/)
- [Biopython](https://biopython.org)


### Installation
ViiR and its dependencies are easily installed via [bioconda](https://bioconda.github.io/index.html) like below:

```
conda create -n viir -c bioconda trinity hmmer wget samtools biopython barrnap
conda activate viir
conda install -c bioconda rsem bioconductor-deseq2 bioconductor-edger r bowtie
git clone https://github.com/YuSugihara/ViiR.git
cd ViiR
pip install . 
```

**If you install RSEM and DESeq2 with other dependencies at the same time, anaconda will take too long time to solve the environment or cannot solve it.** Therefore, we highly recommend the users to install them separately via bioconda.


If the error ```samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory``` appeared, the symbolic link might solve the error. 

```
ln -s ~/miniconda3/envs/viir/lib/libcrypto.so.3 ~/miniconda3/envs/viir/lib/libcrypto.so.1.0.0
```

Please change the path ```~/miniconda3/envs/viir/lib``` to your environment.



## Usage

```
usage: viir -l <FASTQ_LIST> -o <OUT_DIR> [-t <INT>]

ViiR version 0.0.1

optional arguments:
  -h, --help          show this help message and exit
  -l , --fastq-list   Fastq list.
  -o , --out          Output directory. Specified name must not
                      exist.
  -t , --threads      Number of threads. [16]
  -a , --adapter      FASTA of adapter sequences. If you don't
                      specify this option, the defaul adapter set
                      will be used.
  --pfam              List of Pfam IDs. If you don't specify
                      this option, the defaul list will be used.
  --SS-lib-type       Type of strand specific library (No/FR/RF). [No]
  --pvalue            Threshold of pvalue in DESeq2. [0.01]
  --max-memory        Max memory used in Trinity. [32G]
  -v, --version       show program's version number and exit
```

### Example 1 : Run ViiR with default settings

```
viir -l sample_list.txt \
     -o result \
```

`-l` : Sample list describing the paired-end FASTQ files.

`-o` : Name of the output directory. Specified name should not exist.

### Example 2 : Run ViiR with more threads and CPU memories

```
viir -l sample_list.txt \
     -o result \
     -t 40 \
     --max-memory 1000G
```

`-l` : Sample list describing the paired-end FASTQ files.

`-o` : Name of the output directory. Specified name should not exist.

`-t` : Number of threads.

`--max-memory` : Maximum memory used for Trinity.

### Example 3 : Run ViiR with strand specific library

```
viir -l sample_list.txt \
     -o result \
     -t 40 \
     --max-memory 1000G \
     --SS-lib-type FR
```

`-l` : Sample list describing the paired-end FASTQ files.

`-o` : Name of the output directory. Specified name should not exist.

`-t` : Number of threads.

`--max-memory` : Maximum memory used for Trinity.

`--SS-lib-type` : Type of strand specific library (FR or RF).
