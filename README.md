
# Matsuda2024lifespan

In this repository you find codes used in Matsuda and Makino (2024) Comparative genomics reveals convergent signals associated with the high metabolism and longevity in birds and bats.

## Requirements

### General

- Unix-like environment
- python3 (>= 3.9.0)
- R (>= 4.3.0)

### Bioinformatics tools

- blast (>= 2.11.0+)
- [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
- [mafft](https://mafft.cbrc.jp/alignment/software/) (>= 7.487)
- [trimAl](http://trimal.cgenomics.org/) (>= 1.4.1)
- [iqtree2](http://www.iqtree.org/)

### Python libraries

- [biopython](https://biopython.org/)?
- [CSUBST](https://github.com/kfuku52/csubst)
- [bithon](https://github.com/ymat2/bithon)
- [ete3](https://github.com/etetoolkit/ete)

### R packages

- [RERconverge](https://github.com/nclark-lab/RERconverge)
- tidyverse
- phytools


## Pipeline

### Get sequences, detect orthologs, and alignment

```sh
## execute select_species.R & prep_master_tree.R
Rscript src/c0_select_species/select_species.R
Rscript src/c0_select_species/prep_master_tree.R

## get sequences from NCBI
qsub src/c1_get_seq/get_seq.sh out/accession.list

## run vs_reference_species blastp
qsub src/c2_blast/make_blast_db.sh accession.list
qsub src/c2_blast/run_blast.sh accession.list GCF_000001405.40
bzip2 blst/*.pep.fa.*

## identify SCO from blastp
qsub src/c3_detect_sco/detect_sco.sh

## alignment and trimming
qsub src/c4_align_trim/align.sh out/sco.list
qsub src/c4_align_trim/trim.sh out/sco.list

## filter out trimed alignments (optional)
# qsub src/c4_align_trim/filter.sh
```

### Run PAML

```sh
qsub src/c5_paml/run_paml.sh out/cds.list
echo -e "symbol\tp-value\tforeground-omega\tbackground-omega" > out/paml_aves_result.tsv
ls cds/*/aves_result.txt | while read x; do [ -s $x ] && cat $x >> out/paml_aves_result.tsv; done
echo -e "symbol\tp-value\tforeground-omega\tbackground-omega" > out/paml_chiroptera_result.tsv
ls cds/*/chiroptera_result.txt | while read x; do [ -s $x ] && cat $x >> out/paml_chiroptera_result.tsv; done
```

### Run CSUBST

```sh
qsub src/c6_csubst/run_csubst.sh out/cds.list
python3 src/c6_csubst/summarize_csubst_result.py -d csubst -o out/csubst_result.tsv
```
