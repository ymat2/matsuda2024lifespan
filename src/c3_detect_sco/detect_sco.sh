#$ -S /bin/bash
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e ./job/detect_sco.log.stderr
#$ -cwd

[ ! -d ./sco/pep ] && mkdir -p ./sco/pep
[ ! -d ./sco/cds ] && mkdir -p ./sco/cds

python3 src/c3_detect_sco/blast2sco.py \
  --directory blst \
  --reference GCF_000001405.40 \
  --min_sp 3 \
  --output out/single_copy_orthologs.tsv

python3 src/c3_detect_sco/separate_seq.py \
  -f seq \
  -st pep \
  -l out/single_copy_orthologs.tsv \
  -o sco/pep

python3 src/c3_detect_sco/separate_seq.py \
  -f seq \
  -st cds \
  -l out/single_copy_orthologs.tsv \
  -o sco/cds

ls sco/pep/ | awk -F '.' '{print $1}' > out/sco.list
