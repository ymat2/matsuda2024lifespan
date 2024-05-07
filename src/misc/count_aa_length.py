import argparse
import glob
from bithon import fs


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--indir")
  parser.add_argument("-o", "--output")
  args = parser.parse_args()

  fasta_files = glob.glob(args.indir + "/*.pep.aln.nogap.fa")

  with open(args.output, "w") as f:
    f.write("symbol\tlength\n")
    for fa in fasta_files:
      symbol = fa.split("/")[-1].split(".")[0]
      length = str(get_aa_length(fa))
      f.write(symbol + "\t" + length + "\n")


def get_aa_length(path):
  length = fs.fasta2dict(path)
  length = [len(aa) for aa in length.values()]
  if len(set(length)) != 1:
    print("Warning:", path, "contains different length AA sequence.")
  return length[0]


if __name__ == "__main__":
  main()
