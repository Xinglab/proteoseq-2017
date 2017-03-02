##Proteotranscriptome pipeline
RNA-seq and proteomics integrated pipeline for identification of peptide evidence of Alu exon or Human specific exons
<br />
contact: ybwang@hust.edu.cn<br />
Please leave message to me via email or report here https://github.com/Xinglab/proteoseq/issues if any questions <br />
update: 2017-03-02

## Requirements:
    1. python 2.7.x or above, samtools
    2. add samtools directory to the $PATH environment variable

## Installation:
    # download using the url: https://github.com/Xinglab/proteoseq/archive/master.zip 
    unzip unzip proteoseq-master.zip
    cd proteoseq-master
    python install.py --homedir /u/home/y/ybwang --install directorToInstall

## USAGE
    Usage: ptp.py -b Aligned.out.sorted.bam -j SJ.tab.out -p proteomicsdir -e HSExonfile -o outdir
    Parameter:
        -b bam file from STAR
        -j SJ.tab.out file from STAR
        -p proteomics dir
        -e exons file (bed format) or use 'None' to search all junctions peptide [default: None]
        -d database file to perform the MS comet search [default: data/Ensembl_Alu_25bp_0.5.unique.sorted.bed]
        -g genomic fasta directory [default: /u/home/f/frankwoe/nobackup/hg19/hg19_by_chrom/]
        -o output directory [default: outdir]
        --l Extend flanking junction ends by this number of bp
        --min-junc-reads Minimum number of reads required spanning the junction

## Example
    # search peptides from Alu exons
    cd ~/scratch/ptptest
    ln -s /u/home/f/frankwoe/nobackup/AST/Yoav_Gilad_YRI/RNA/star_output/GM18486.rna/ RNA
    ln -s /u/home/y/ybwang/nobackup-yxing/data/Yoav_Gilad_proteom/GM18486/ PRO
    python ptp.py -b RNA/Aligned.out.sorted.bam -j RNA/SJ.out.tab -p PRO/ -e data/Ensembl_Alu_25bp_0.5.unique.sorted.bed -d data/UP000005640_9606_additional_cdhit1.fasta -g /u/home/f/frankwoe/nobackup/hg19/hg19_by_chrom/

    # or search for all junctions peptides:
    python ptp.py -b RNA/Aligned.out.sorted.bam -j RNA/SJ.out.tab -p PRO/

    # or submit the job using qsub
    qsub -cwd -V -N PTP -l h_data=30G,h_rt=3:00:00 -M eplau -m bea ./example_submit.sh

## Update (2017-02-23):
    1. modify 'percolator_test2parellel.py' to 'percolator_triesearch.py' to speed up the re-mapping step
    1. output the result to file
    2. modify the output on console



