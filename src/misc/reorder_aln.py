import argparse
import pandas as pd
from bithon import fs

def collapse(seq_dict):
  from collections import Counter
  new_seq_dict = {sp: "" for sp in seq_dict}
  for i in range(len(seq_dict["Homo_sapiens"])):
    seq_i = [seq_dict[sp][i] for sp in seq_dict]
    cons = Counter(seq_i)
    cons = cons.most_common(1)[0][0]
    for sp, seq in seq_dict.items():
      if seq[i] == cons:
        new_seq_dict[sp] += "."
      else:
        new_seq_dict[sp] += seq[i]
  return new_seq_dict

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment")
parser.add_argument("-t", "--tsv")
parser.add_argument("--collapse", default=1.0)
args = parser.parse_args()

df = pd.read_csv(args.tsv, sep = '\t')
sp2class = dict(zip(df["Organism Name"], df["class"]))
sp2class = {sp.replace(" ", "_"): sp2class[sp] for sp in sp2class}

fasta = fs.fasta2dict(args.alignment)
fasta = collapse(fasta)

with open("reordered.fa", "w") as f:
  for i in ["Mammalia", "Chiroptera", "Reptilia", "Aves"]:
    f.write(">"+i+"\n")
    f.write("-"*len(fasta["Homo_sapiens"])+"\n")
    for sp, seq in fasta.items():
      if sp2class[sp] == i:
        f.write(">"+sp+"\n")
        f.write(seq+"\n")
