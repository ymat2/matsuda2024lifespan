import argparse
from ete3 import PhyloTree, TreeStyle


def main():

  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--alignment")
  parser.add_argument("-t", "--tree")
  args = parser.parse_args()

  tree, t_style = get_tree_style(args.alignment, args.tree)
  tree.show(tree_style = t_style)


def get_tree_style(aln_path, tree_path):

  with open(aln_path) as f:
    aln = f.read()
  with open(tree_path) as f:
    tree = f.read()

  tree = PhyloTree(tree)
  tree.link_to_alignment(aln)

  return tree, TreeStyle()


if __name__ =="__main__":
  main()
