#$ -S /bin/bash
#$ -t 1-18874:1
#$ -cwd
#$ -V
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -o /dev/null
#$ -e ./job/
#$ -tc 45


seq_ids=()
while read -r x; do
  seq_ids+=($x)
done < $1
seq_id=${seq_ids[$SGE_TASK_ID-1]}

[ ! -d ./csubst/${seq_id} ] && mkdir -p ./csubst/${seq_id}

echo ${seq_id} 1>&2

python3 src/c6_csubst/prep_csubst_infile.py \
	-a cds/${seq_id}/${seq_id}.cds.aln.nogap.fa \
	-d out/species.tsv \
	-t out/RERconverge_master_tree.nogap.nwk \
	-o csubst/${seq_id}

cd csubst/${seq_id}

csubst analyze \
	--alignment_file ../../cds/${seq_id}/${seq_id}.cds.aln.nogap.fa \
	--rooted_tree_file pruned_tree.nwk \
	--foreground foreground.txt \
	--exhaustive_until 1 \
	--cutoff_stat 'OCNany2spe,0|omegaCany2spe,0' \
	--fg_exclude_wg yes \
	--fg_stem_only yes \
	--iqtree_exe iqtree2
