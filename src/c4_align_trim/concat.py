import argparse
from Bio.Seq import Seq
from bithon import fs


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--pep")
  parser.add_argument("-c", "--cds")
  parser.add_argument("-o", "--outfile")
  args = parser.parse_args()

  pep_seq = fs.fasta2dict(args.pep)
  pep_seq = {k+"_pep": v for k, v in pep_seq.items()}
  nuc_seq = fs.fasta2dict(args.cds)
  nuc_seq = {k+"_cds": nuc2pep(v) for k, v in nuc_seq.items()}

  out_seq = {}
  out_seq |= pep_seq
  out_seq |= nuc_seq

  fs.write_fasta(out_seq, args.outfile)


def nuc2pep(dna):
	if len(dna)%3 == 0:
		dna = dna
	elif (len(dna)+1)%3 == 0:
		dna = dna + "N"
	else:
		dna = dna + "NN"
	seq = Seq(dna)
	seq = str(seq.translate(stop_symbol="*"))
	return seq.rstrip("*")


if __name__ == "__main__":
  main()
