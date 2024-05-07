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

singularity exec /usr/local/biotools/b/blast:2.14.0--h7d5a4b4_1 blastp \
  -outfmt 6 \
  -evalue 1e-4 \
  -db blst/${seq_id}.pep.fa \
  -query seq/$2.pep.fa \
  -out blst/ref_to_${seq_id}.blastp \
  -max_target_seqs 1

singularity exec /usr/local/biotools/b/blast:2.14.0--h7d5a4b4_1 blastp \
  -outfmt 6 \
  -evalue 1e-4 \
  -db blst/$2.pep.fa \
  -query seq/${seq_id}.pep.fa \
  -out blst/${seq_id}_to_ref.blastp \
  -max_target_seqs 1
