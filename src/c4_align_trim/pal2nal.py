import argparse
from bithon import fs

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-a", "--aln")
	parser.add_argument("-c", "--cds")
	parser.add_argument("-o", "--outpath")
	args = parser.parse_args()

	gn2aln = fs.fasta2dict(args.aln)
	gn2cds = fs.fasta2dict(args.cds)

	pep_seq = {k.removesuffix("_pep"): v for k, v in gn2aln.items() if "_pep" in k}
	cds_seq = {k.removesuffix("_cds"): v for k, v in gn2aln.items() if "_cds" in k}
	pep_seq = amb2gap(pep_seq)
	cds_seq = replace_dif_aa(pep_seq, cds_seq)

	nuc_seq = {gn: align(cds_seq[gn], gn2cds[gn]).lower() for gn in cds_seq}
	nuc_seq = {gn: del_last_stop_codon(nuc_seq[gn]) for gn in nuc_seq}

	fs.write_fasta(pep_seq, args.outpath+".pep.aln.fa")
	fs.write_fasta(cds_seq, args.outpath+".translated.aln.fa")
	fs.write_fasta(nuc_seq, args.outpath+".cds.aln.fa")


def amb2gap(dct):
  dct_amb2gap = dct
  amb_aa = ["B", "Z", "J", "U", "O"]
  for k in dct_amb2gap:
    for chr in amb_aa:
      dct_amb2gap[k] = dct_amb2gap[k].replace(chr, "-")
  return dct_amb2gap


def replace_dif_aa(pep_dct, cds_dct):
	cds_replaced = {}
	for gn in pep_dct.keys():
		pseq, cseq = pep_dct[gn], cds_dct[gn]
		cds_seq_replaced = ""
		for i in range(len(pseq)):
			if pseq[i] == cseq[i]:
				cds_seq_replaced += cseq[i]
			elif cseq[i] == "-":
				cds_seq_replaced += "-"
			else:
				cds_seq_replaced += "*"
		cds_replaced[gn] = cds_seq_replaced
	return cds_replaced

def align(aln, cds):
	diff_ = len(aln)*3 - len(cds)
	if diff_ > 0:
		cds = cds + "n"*diff_
	pos = 0
	aligned_cds = ""
	for aa in aln:
		if aa == "-":
			aligned_cds += "---"
		elif aa == "*":
			aligned_cds += "nnn"
			pos += 3
		else:
			codon = cds[pos:pos+3]
			if "N" in codon or "n" in codon:
				aligned_cds += "nnn"
			else:
				aligned_cds += codon
			pos += 3
	return(aligned_cds)


def del_last_stop_codon(_cds):
	if _cds[-3:] in ["TAA", "TGA", "TAG"]:
		return _cds[:-3]
	else:
		return _cds


if __name__ == "__main__":
	main()
