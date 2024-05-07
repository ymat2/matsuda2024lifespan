from ete3 import Tree
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-t", "--tree")
  parser.add_argument("-d", "--tsv")
  parser.add_argument("-a", "--aln")
  parser.add_argument("-o", "--outdir")
  args = parser.parse_args()

  master_tree = Tree(args.tree)
  anc = master_tree.get_common_ancestor("Anolis_carolinensis", "Taeniopygia_guttata")
  master_tree.set_outgroup(anc)
  tips = read_tips(args.aln)

  master_tree.prune(tips)
  master_tree.write(outfile = args.outdir+"/pruned_tree.nwk", format=5)

  species = read_sp_class(args.tsv)
  write_foreground(species, args.outdir+"/foreground.txt")


def read_tips(aln):
  tips = []
  with open(aln) as f:
    for line in f:
      if line[0] == ">":
        tips.append(line.rstrip("\n")[1:])
  return tips


def read_sp_class(tsv):
  sp2class = {}
  with open(tsv) as f:
    next(f)
    for line in f:
      sp_name = line.split("\t")[1].replace(" ", "_")
      sp_class = line.split("\t")[7]
      sp2class[sp_name] = sp_class
  return sp2class


def write_foreground(dct, outpath):
  with open(outpath, "w") as f:
    for sp in dct:
      if dct[sp] == "Aves":
        f.write("1\t"+sp+"\n")
      elif dct[sp] == "Chiroptera":
        f.write("2\t"+sp+"\n")
      else:
        continue


if __name__ == "__main__":
  main()
