#$ -S /bin/bash
#$ -t 1-132:1
#$ -cwd
#$ -o /dev/null
#$ -e ./job
#$ -tc 100

seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

module load singularity

[ ! -d "blst" ] && mkdir ./blst

singularity exec /usr/local/biotools/b/blast:2.14.0--h7d5a4b4_1 makeblastdb \
  -in seq/${seq_id}.pep.fa \
  -out blst/${seq_id}.pep.fa \
  -dbtype prot \
  -hash_index
