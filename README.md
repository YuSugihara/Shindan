# Shindan User Guide
#### Version 0.0.1

## Table of contents
- [What is Shindan?](#What-is-Shindan)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation using bioconda](#Installation-using-bioconda)


## What is Shindan?

Shindan is a pipeline to quickly identify a causal plant virus in the field. 

## Installation
### Dependencies
#### Softwares
- [Python3](https://www.python.org/downloads/)
- [wget](https://www.gnu.org/software/wget/)
- [Trinity package](https://github.com/trinityrnaseq/trinityrnaseq)
  + [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)
  + [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  + [Samtools](http://www.htslib.org/doc/samtools.html)
  + [R](https://www.r-project.org/)
  + [DESeq2](https://bioconductor.org/packages/3.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [HMMER](http://hmmer.org/)

We highly recommend the users to download Trinity via [bioconda](https://bioconda.github.io/).
Other softwares requered in upstream or downstream analysis can be downloaded with Trinity via bioconda. 
