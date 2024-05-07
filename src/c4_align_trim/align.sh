#$ -S /bin/bash
#$ -t 1-19024:1
#$ -cwd
#$ -o /dev/null
#$ -e ./job/align.log.stderr
#$ -tc 250


seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

[ ! -d ./aln ] && mkdir ./aln

python3 src/c4_align_trim/concat.py \
	-p sco/pep/${seq_id}.pep.fa \
	-c sco/cds/${seq_id}.cds.fa \
	-o aln/${seq_id}.concat.fa

singularity exec /usr/local/biotools/m/mafft:7.520--h031d066_2 mafft \
	--auto --anysymbol --quiet \
	aln/${seq_id}.concat.fa > aln/${seq_id}.concat.aln.fa

python3 src/c4_align_trim/pal2nal.py \
	-a aln/${seq_id}.concat.aln.fa \
	-c sco/cds/${seq_id}.cds.fa \
	-o aln/${seq_id}

rm aln/${seq_id}.concat.aln.fa
