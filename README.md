# ViiR User Guide
#### Version 0.0.1

## Table of contents
- [What is ViiR?](#What-is-ViiR)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation via bioconda](#Installation-via-bioconda)
- [Usage](#Usage)


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


### Installation via bioconda
ViiR and its dependencies are easily installed via [bioconda](https://bioconda.github.io/index.html) like below:

```
conda create -n viir -c bioconda trinity hmmer wget samtools>=1.9
conda activate viir
conda install -c bioconda rsem bioconductor-deseq2 bioconductor-edger r>=4.1
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
  -t , --threads      Number of threads.
  -a , --adapter      FASTA of adapter sequences. If you don't
                      specify this option, the defaul adapter set
                      will be used.
  --pfam              List of Pfam IDs. If you don't specify
                      this option, the defaul list will be used.
  --SS-lib-type       Type of strand specific library (No/FR/RF) [No].
  --pvalue            Threshold of pvalue in DESeq2. [0.01]
  --max-memory        Max memory to use by Trinity. [32G]
  -v, --version       show program's version number and exit
```

### Example 1 : Run ViiR with default settings

```
viir -l sample_list.txt \
     -o result \
```
`-l` : Sample list describing the paired-end FASTQ files.

`-o` : Name of the output directory. Specified name should not exist.

### Example 2 : Run ViiR with default settings

```
viir -l sample_list.txt \
     -o result \
     -t 40 \
     --max-memory 1000G
```
