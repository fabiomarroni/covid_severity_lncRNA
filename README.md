# covid_severity_lncRNA

This is an R function that enables you parsing results of the meta-analsyis of the COVID-19 Host Genetics Initiative (https://www.covid19hg.org/) and annotate then SNPs for which the closest gene is a lncRNA.

# Instructions
1) Choose a working directory, move there using the`cd` command and save the function `close2lncRNA.r` in the directory. 

1) Download the Human Reference Genome gff file for version 38 by typing:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
```
2)  Unzip file (if you are on windows you may just download mobaXterm which has an integrated gzip utility). 

```
gunzip GCF_000001405.39_GRCh38.p13_genomic.gff.gz
```

3) Download the result file from COVID-19 Host Genetics Initiative (here we are downloading the B1_ALL table, a GWAS of Hospitalized COVID-19 vs non-hospitalize COVID-19; you can choose any of the "filtered" file):
```
wget https://storage.googleapis.com/covid19-hg-public/20200915/results/20201020/COVID19_HGI_B1_ALL_20201020.txt.gz_1.0E-5.txt
```
4) Open R (or Rstudio) and move to the working directory using the setwd() command

5) Run the command source("close2lncRNA.r") 

6) Enjoy! The results are in the file `ncRNA_COVID19_HGI_B1_ALL_20201020.txt.gz_1.0E-5.txt`. The file is the same as the one you downloaded from the COVID-19 Host Genetics Initiative, with a last column (lncRNA). The column will show the name of the lncRNA that represent the closest gene to the SNP of the row you are vieweing. If the closest gene to the SNP is not a lncRNA (e.g. is a "normal" mRNA), the column will show the value "No".
