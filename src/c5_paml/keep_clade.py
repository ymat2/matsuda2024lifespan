import argparse
import pandas as pd

from bithon import fs


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--aln")
  parser.add_argument("-t", "--tsv")
  parser.add_argument("-u", "--use_species", nargs="*")
  parser.add_argument("-o", "--outfile")
  args = parser.parse_args()

  inaln = fs.fasta2dict(args.aln)
  sp2cls = tsv2dct(args.tsv)
  outaln = keep_clade(inaln, sp2cls, args.use_species)
  fs.write_fasta(outaln, args.outfile)


def tsv2dct(tsv):
  df = pd.read_csv(tsv, sep = "\t")
  sp2cls = dict(df[["Organism Name", "class"]].values)
  sp2cls = {k.replace(" ", "_"): v for k, v in sp2cls.items()}
  return sp2cls


def keep_clade(alignment, sp2cls, clade: list):
  new_aln = {}
  for sp in alignment:
    _class = sp2cls[sp]
    if _class in clade:
      new_aln[sp] = alignment[sp]
  return new_aln


if __name__ == "__main__":
  main()
