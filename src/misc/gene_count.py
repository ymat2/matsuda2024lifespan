import argparse
import glob
import pandas as pd

from bithon import fs

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-d", "--dir")
  parser.add_argument("-t", "--tsv")
  parser.add_argument("-o", "--outdir")
  args = parser.parse_args()

  sco_files = glob.glob(args.dir+"/*.fa")
  print("Recruit", len(sco_files) ,"files.")
  df = pd.read_csv(args.tsv, sep = "\t")
  species = df["Organism Name"].to_list()
  species = [sp.replace(" ", "_") for sp in species]
  print("Number of species:", len(species))
  print(species)

  print("Counting genes per species...")
  gc_per_species = count_per_species(sco_files, species)

  print("Counting genes per clades...")
  gc_per_clade = count_per_clade(sco_files, df)

  print("Writing output files...")
  with open(args.outdir+"/gene_count_per_species.tsv", "w") as f:
    f.write("symbol\t")
    f.write("\t".join(species)+"\n")
    for gn, gc in gc_per_species.items():
      f.write(gn+"\t"+"\t".join(list(map(str, gc)))+"\n")

  with open(args.outdir+"/gene_count_per_clade.tsv", "w") as f:
    f.write("symbol\tAves\tChiroptera\tMammalia\tReptilia\n")
    for gn, gc in gc_per_clade.items():
      line = list(map(str, [gc["Aves"], gc["Chiroptera"], gc["Mammalia"], gc["Reptilia"]]))
      f.write(gn+"\t"+"\t".join(line)+"\n")


def count_per_species(files, species):
  gene_count_per_species = dict()
  for f in files:
    symbol = ".".join(f.split("/")[-1].split(".")[:-1])
    seqs = fs.fasta2dict(f)
    gene_count_per_species[symbol] = [1 if sp in seqs else 0 for sp in species]
  return gene_count_per_species


def count_per_clade(files, df):
  gene_count_per_clade = dict()
  sp2class = dict(zip(df["Organism Name"], df["class"]))
  for f in files:
    symbol = ".".join(f.split("/")[-1].split(".")[:-1])
    seqs = fs.fasta2dict(f)
    gene_count_per_clade[symbol] = {"Aves": 0, "Chiroptera": 0, "Mammalia": 0, "Reptilia": 0}
    for sp in seqs.keys():
      sp = sp.replace("_", " ")
      if sp in sp2class:
        gene_count_per_clade[symbol][sp2class[sp]] += 1
  return gene_count_per_clade


if __name__ == "__main__":
  main()
