#$ -S /bin/bash
#$ -t 1-132:1
#$ -cwd
#$ -V
#$ -o /dev/null
#$ -e ./job
#$ -tc 100

seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

[ ! -d ncbi_dataset ] && mkdir ./ncbi_dataset
[ ! -d seq ] && mkdir ./seq

datasets download genome accession ${seq_id} --filename ${seq_id}.zip --include gtf,cds,protein
unzip ${seq_id}.zip -d ncbi_dataset/${seq_id}
rm -r ${seq_id}.zip

bithon gls --indir ncbi_dataset/${seq_id}/ncbi_dataset/data/${seq_id} --outdir seq/${seq_id}
rm -r ncbi_dataset
