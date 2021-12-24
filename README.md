# ViiR User Guide
#### Version 0.0.1

## Table of contents
- [What is ViiR?](#What-is-ViiR)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation via bioconda](#Installation-via-bioconda)


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
  + [R](https://www.r-project.org/)
  + [DESeq2](https://bioconductor.org/packages/3.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [HMMER](http://hmmer.org/)

We highly recommend the users to download Trinity via [bioconda](https://bioconda.github.io/).
Other softwares requered in upstream or downstream analyses can be downloaded with Trinity via bioconda.


### Installation via bioconda
You can easily install the dependencies of ViiR via [bioconda](https://bioconda.github.io/index.html).
```
conda install -c bioconda trinity hmmer wget
```

If you want to create ViiR specific environment:
```
conda create -n viir -c bioconda trinity hmmer wget
```

ViiR can be installed using the following commands.
```
git clone https://github.com/YuSugihara/ViiR.git
cd ViiR
pip install . 
```

