import argparse
import os
from bithon import fs

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--aln")
  parser.add_argument("--minsp", type=int, default=3)
  parser.add_argument("--minseqlen", type=int, default=10)
  args = parser.parse_args()

  aln = fs.fasta2dict(args.aln)
  first_key = list(aln.keys())[0]
  seq = aln.get(first_key, "NA")

  if len(aln) < args.minsp:
    print(args.aln, ": Number of species is less than %d." % (args.minsp))
    os.remove(args.aln)

  if seq == "NA":
    print(args.aln, "has no sequences.")
  elif len(seq) < args.minseqlen:
    print(args.aln, ": Length of sequences is less than %d." % (args.minseqlen))
    os.remove(args.aln)


if __name__ == "__main__":
  main()
