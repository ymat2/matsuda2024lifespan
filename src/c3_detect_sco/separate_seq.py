import argparse
import glob
from bithon import fs


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-f", "--fastadir")
  parser.add_argument("-l", "--sco")
  parser.add_argument("-st", "--seq_type", choices=["pep", "cds"])
  parser.add_argument("-o", "--outdir")
  args = parser.parse_args()

  seq_all = {}
  FASTA_FILES = glob.glob(args.fastadir+"/*."+args.seq_type+".fa")
  for f in FASTA_FILES:
    dct = fs.fasta2dict(f)
    seq_all.update(dct)

  with open(args.sco) as f:
    for line in f:
      genes = line.rstrip("\n").split("\t")
      hs_gene_sym = genes[0].removesuffix("_Homo_sapiens")
      dct_sco = {as_sp_name(gn): seq_all[gn] for gn in genes}
      outfile = args.outdir+"/"+hs_gene_sym+"."+args.seq_type+".fa"
      fs.write_fasta(dct_sco, outfile)


def as_sp_name(gn):
  gns = gn.split("_")
  sp = gns[-2]+"_"+gns[-1]
  return sp


if __name__ == "__main__":
  main()
