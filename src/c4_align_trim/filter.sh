#$ -S /bin/bash
#$ -cwd
#$ -o ./job/filt.log.stdout
#$ -e ./job/filt.log.stderr


ls trimal/ | while read x; do
  python3 src/c4_align_trim/filter_align_file.py \
  	-a trimal/$x \
  	--minsp 3 \
  	--minseqlen 10
done
