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
- [R](https://www.r-project.org/)
- [DESeq2](https://bioconductor.org/packages/3.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [RSEM](https://deweylab.github.io/RSEM/)
- [HMMER](http://hmmer.org/)

We highly recommend the users to download Trinity via [bioconda](https://bioconda.github.io/).
Other softwares requered in upstream or downstream analyses can be downloaded with Trinity via bioconda.


### Installation via bioconda
You can easily install the dependencies of ViiR via [bioconda](https://bioconda.github.io/index.html).

```
conda create -n viir -c bioconda trinity hmmer wget samtools>=1.9
conda activate viir
conda install -c bioconda rsem bioconductor-deseq2 bioconductor-edger r>=4.1
```

If you install RSEM and DESeq2 with other dependencies at the same time, anaconda will take a long time to solve the environment or cannot solve the environment.

ViiR can be installed using the following commands.
```
git clone https://github.com/YuSugihara/ViiR.git
cd ViiR
pip install . 
```


If you will face an error like ```samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory```, you can solve the problem using the symbolic link like below. Please change the PATH to the environment of ViiR in anaconda.
```
ln -s ~/miniconda3/envs/viir/lib/libcrypto.so.3 ~/miniconda3/envs/viir/lib/libcrypto.so.1.0.0
```

