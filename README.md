# `RNAsik-pipe` easy and quick RNA-seq(uencing) pipeline

> RNAsik-pipe is written in [BigDataScrip](http://pcingola.github.io/BigDataScript/) (BDS) language
> There are many advantages in using BDS, but main ones are:
>
>  1. One script works on your local machine, on your remote server and on your cluster
>  2. Checkpoints - never need to run your job from the start in case it stopped. Start from where you left off
>  3. Remote access to your data. Don't need to worry about `scp`ing you data to the server, just point RNAsik-pipe to your cloud store if you like
>
> As for `RNAsik-pipe` itself it makes your life easy. One single command line can give you list of differrentially expressed genes.
> Simply give your FASTQ files to RNAsik-pipe with reference files and press go ! [Enter]

## Content

- [Introduction](#introduction)
  - [FASTQ files explained](#fastq-files-explained)
  - [Get RNAseq metrics](#get-rnaseq-metrics)
  - [get your counts](#get-your-counts)
- [Installation](#installation)
  - [The easy way](#the-easy-way)
  - [The other way](#the-other-way)
- [Prerequisites](#prerequisites)
- [Quick start](#quick-start)
- [User manual](#user-manual)
- [Release notes](#release-notes)
  - [Version 1.2](#version-1.2)

## Introduction

Your standard [RNA-seq](https://en.wikipedia.org/wiki/RNA-Seq) workflow as follows:

- get [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) reads
- align them to your reference genome (your reference genome is given in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file). [STAR aligner](https://github.com/alexdobin/STAR) is used in `RNAsik-pipe` for that
- count how many reads have mapped to the gene feature. You need [SAM/BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) for this. [featureCounts](http://subread.sourceforge.net/) is used in `RNAsik-pipe` for that
- Then you can choose to use [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) in [R](https://en.wikipedia.org/wiki/R_(programming_language) to do differential gene expression analysis Or you can simply upload your counts file to [Degust](http://victorian-bioinformatics-consortium.github.io/degust/), which is interactive and user friendly web tool

**OR you can just use** `RNAsik-pipe`

`RNAsik-pipe` is a simple tools that wraps several other tools together to obstruct away the same repetitive work of retyping the same commands over the same number of tools just to get your differentially expressed genes. Simply provide `RNAsik-pipe` with [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files directory and your reference files such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format) and [GTF](http://mblab.wustl.edu/GTF22.html) and press Go !

Even though most tools used inside `RNAsik-pipe` can handle [GFF](https://en.wikipedia.org/wiki/General_feature_format) at this stage only [GTF](http://mblab.wustl.edu/GTF22.html) file is supported by the pipeline.

One important note about [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, `RNAsik-pipe` assumes that if the data is paired-end than it will have upper case **R** following by the number 1 or 2. `RNAsik-pipe` further assumes that both **R1** and **R2** have and underscore prefix - `_R1`, `_R2`. If your data is single-end and just has `_R1` in the file name, this is fine and `RNAsik-pipe` will handle those files correctly. `RNAsik-pipe` will also handle single-end data without any of **R's** in the file name correctly, but you might have to specify your own `-fqRegex` for more refer to [User manual](supplementary/docs.md).

### FASTQ files explained

  Your raw data will always come in FASTQ format. The number of FASTQ files will depend on many things
  including:

  - Number of samples 
  - Number of replicates 
  - Your sample was split into different lanes
  - Your are sequencing paired-end data

  Your FASTQ files might reside in one directory e.g directory per experiment

  ![fqDir](supplementary/rawDataDir.png)

  In this example sample 14-09157, which might WT is split across two lanes `L001` and `L002`, which is 
  typical of Illumina sequencing. That is two files, one for each lane. We can also see that this is
  paire end data. That means for each file there has to be a pair file, that is R1 and R2. In summary
  single sample, e.g WT is covered by four FASTQ files:
     - 14-09157_L001_R1.fastq.gz
     - 14-09157_L001_R2.fastq.gz
                 AND 
     - 14-09157_L002_R1.fastq.gz
     - 14-09157_L002_R2.fastq.gz
  you can use `cat` command to concatenate files across different lanes
  e.g `cat 14-09157_L001_R1.fastq.gz 14-09157_L002_R1.fastq.gz > 14-09157_merged.fastq.gz` or you can merge
  BAM files with `samtools` later. `STAR` aligner can merger on the fly and what RNAsik-pipe is using.

  OR 

  Each sample might be put into its own subdirectory e.g

  ![test](supplementary/rawDataDirs.png)
  
  In this example each sample is placed into its own directory. Now we can see directory `Sample_14-09157`, 
  which will hold four files for described above.

  `RNAsik-pipe` can work with either of those two options. Specify your "root" directory either with `-fqDir`
  or `-fqDirs` options.

### Get RNAseq metrics

  [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc) requires you BAM to be sorted,
  reordered and have duplicates marked. Here is detailed [instructions](supplementary/RNAseQC-manual.pdf)
  for how to prepare your BAMs files for RNAseQC, BUT the good news is you don't even need to worry about this!
  `RNAsik-pipe` takes care of long and laborious BAM file manipulation for RNA-SeQC tools, just flag 
  `-prePro` to get your BAMS in the right shape and `-RNAseQC` to get the actual report

### Get your read counts

  Do you want to do differential gene expression analysis..? just flag `-count` and you will get your counts

## Prerequisites

- [BigDataScript](http://pcingola.github.io/BigDataScript/download.html)
- [STAR aligner](https://github.com/alexdobin/STAR/releases)
- [featureCounts](http://subread.sourceforge.net/)
- [samtools](http://www.htslib.org/download/)
- [Picard tools](http://broadinstitute.github.io/picard/)
- [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc)

### System requirements 

- At least 100 Gb of disk space. [STAR aligner](https://github.com/alexdobin/STAR/releases) you will need 100 Gb of disk space during making an index step. The index itself isn't 100 Gb, but that much space is required. 
- At least 32 Gb of RAM

## Installation

### The easy way

**It is recommended that you use ansible to install RNAsik-pipe with all dependencies**

The easiest way to get `RNAsik-pipe` installed together with all dependencies is to use [ansible RNAsik-playbook](https://github.com/MonashBioinformaticsPlatform/ansible-RNAsik-REST).
All you will need to do is install ansible `sudo pip install ansible`, clone the [ansible RNAsik-playbook](https://github.com/MonashBioinformaticsPlatform/ansible-RNAsik-REST)
`git clone https://github.com/MonashBioinformaticsPlatform/ansible-RNAsik-REST.git` and run ansible `sudo ansible-playbook -i local_server_hosts site.yml`. [Look here for more help](https://github.com/MonashBioinformaticsPlatform/ansible-RNAsik-REST).

### The other way

Make sure to install [BigDataScript](http://pcingola.github.io/BigDataScript/) first. Follow [BDS installation instructions](http://pcingola.github.io/BigDataScript/download.html)

Get [latest stable release](https://github.com/MonashBioinformaticsPlatform/RNAsik-pipe/releases) by 
downloading `*tar.gz` file.

1. Locate your `*tar.gz` file
2. `tar zxvf *tar.gz file` 
3. You `RNAsik-pipe` executable file is located in `src` directory
 
**If you like to get developing version**

`git clone https://github.com/MonashBioinformaticsPlatform/RNAsik-pipe.git`

## Quick start

- If you have a `module` system on you server/cluster make sure to `module load` all required tools e.g `module load STAR`

- run `RNAsik-pipe` by pointing to the executable file. e.g if you downloaded `*tar.gz` file into `Downlaods`
directory and unpacked there, then run `RNAsik-pipe` from anywhere as such `~/Downloads/RNAsik-pipe/src/RNAsik-pipe`

- To align `RNAsik-pipe -star -fqRegex A -genomeIndex path/to/yourIndexDirectory` 
- To get counts `RNAsik-pipe -count -gtfFile path/to/yourGTFfile`
- To get RNA-SeQC report `RNAsik-pipe -prePro -fastaRef path/to/yourFASTAreference-file -RNAseQC`

**You can simply specify all of the options at the start and let `RNAsik-pipe` to do everything for your**

- `RNAsik-pipe` will guide you through. `RNAsik-pipe` will let you know if you have forgotten any files needed for your run. 

## User manual

- [User manual](supplementary/docs.md)

## Caveats 

- At the moment `STAR aligner` has fixed options, here is the default:

```BASH
STAR --runThreadN 26 \
     --genomeDir $genomeIndex \
     --outSAMtype BAM Unsorted \
     --outSAMattrRGline ID:001 CN:AGRF DS:RNA-seq PL:ILLUMINA PM:MiSeq SM:$uniqueName \
     --outSAMunmapped Within \
     --readFilesCommand zcat \
     --readFilesIn $read1 $read2 \
     --outFileNamePrefix $preFix
```

- Right now only `featureCounts` is supported for read counting
- At the moment there isn't an option to choose the strand direction for read counts. `RNAsik-pipe` simply
counts using both stranded NO and stranded REVERSE options

## Release notes

- [Version 1.2](supplementary/releaseNotes1.2.md)
