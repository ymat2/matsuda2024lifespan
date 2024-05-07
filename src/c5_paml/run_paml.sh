#$ -S /bin/bash
#$ -t 1-18874:1
#$ -cwd
#$ -V
#$ -l s_vmem=24G
#$ -l mem_req=24G
#$ -o /dev/null
#$ -e ./job/
#$ -tc 270

seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

echo ${seq_id} 1>&2

mkdir -p cds/${seq_id}/Aves
mkdir -p cds/${seq_id}/Chiroptera

python3 src/c5_paml/keep_clade.py \
  -a cds/${seq_id}/${seq_id}.cds.aln.nogap.fa \
  -t out/species.tsv \
  --use_species Aves Mammalia Reptilia \
  -o cds/${seq_id}/Aves/${seq_id}.aves.fa
python3 src/c5_paml/keep_clade.py \
  -a cds/${seq_id}/${seq_id}.cds.aln.nogap.fa \
  -t out/species.tsv \
  --use_species Chiroptera Mammalia Reptilia \
  -o cds/${seq_id}/Chiroptera/${seq_id}.chiroptera.fa

iqtree2 -s cds/${seq_id}/Aves/${seq_id}.aves.fa  --redo
iqtree2 -s cds/${seq_id}/Chiroptera/${seq_id}.chiroptera.fa --redo

python3 src/c5_paml/paml_branch.py \
  -a cds/${seq_id}/Aves/${seq_id}.aves.fa \
  -t cds/${seq_id}/Aves/${seq_id}.aves.fa.treefile \
  --tsv out/species.tsv \
  --foreground Aves \
  --codeml_path /home/yukimatsu/.guix-profile/bin/ > cds/${seq_id}/aves_result.txt
python3 src/c5_paml/paml_branch.py \
  -a cds/${seq_id}/Chiroptera/${seq_id}.chiroptera.fa \
  -t cds/${seq_id}/Chiroptera/${seq_id}.chiroptera.fa.treefile \
  --tsv out/species.tsv \
  --foreground Chiroptera \
  --codeml_path /home/yukimatsu/.guix-profile/bin/ > cds/${seq_id}/chiroptera_result.txt
