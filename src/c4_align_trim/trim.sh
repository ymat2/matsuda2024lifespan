#$ -S /bin/bash
#$ -t 1-19024:1
#$ -cwd
#$ -o /dev/null
#$ -e ./job/trim.log.stderr
#$ -tc 250


seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

[ ! -d ./trimal ] && mkdir ./trimal

if [ -s ./aln/${seq_id}.pep.aln.fa ]; then
	singularity exec /usr/local/biotools/t/trimal:1.4.1--0 trimal \
		-in aln/${seq_id}.pep.aln.fa \
		-out trimal/${seq_id}.pep.aln.tmp.fa \
		-resoverlap 0.05 \
		-seqoverlap 90

	singularity exec /usr/local/biotools/t/trimal:1.4.1--0 trimal \
		-in trimal/${seq_id}.pep.aln.tmp.fa \
		-out trimal/${seq_id}.pep.aln.nogap.fa \
		-nogaps

	singularity exec /usr/local/biotools/t/trimal:1.4.1--0 trimal \
		-in trimal/${seq_id}.pep.aln.tmp.fa \
		-out trimal/${seq_id}.pep.aln.autom.fa \
		-automated1

	rm trimal/${seq_id}.pep.aln.tmp.fa
fi

if [ -s ./aln/${seq_id}.cds.aln.fa ]; then
	singularity exec /usr/local/biotools/t/trimal:1.4.1--0 trimal \
		-in aln/${seq_id}.cds.aln.fa \
		-out trimal/${seq_id}.cds.aln.tmp.fa \
		-resoverlap 0.05 \
		-seqoverlap 90

	singularity exec /usr/local/biotools/t/trimal:1.4.1--0 trimal \
		-in trimal/${seq_id}.cds.aln.tmp.fa \
		-out trimal/${seq_id}.cds.aln.nogap.fa \
		-nogaps

	rm trimal/${seq_id}.cds.aln.tmp.fa
fi
