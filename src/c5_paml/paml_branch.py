import argparse
import os
import warnings
import pandas as pd
from bithon import fs
from ete3 import EvolTree

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-a", "--aln")
  parser.add_argument("-t", "--tree")
  parser.add_argument("--tsv")
  parser.add_argument("-f", "--foreground", nargs="*")
  parser.add_argument("--codeml_path")
  parser.add_argument("-o", "--outdir")
  args = parser.parse_args()

  with open(args.aln) as f:
    aln = f.read()

  with open(args.tree) as f:
    tree = f.read()

  tree = EvolTree(tree, binpath=args.codeml_path)
  tree.link_to_alignment(aln)
  if args.outdir:
    workdir = os.path.dirname(args.aln) + "/" + args.outdir
  else:
    workdir = os.path.dirname(args.aln)
  tree.workdir = workdir

  marks = get_foreground_species(args.aln, args.tsv, args.foreground)

  if tree.check_monophyly(marks, target_attr="name")[0]:
    mark_foreground_m(marks, tree)
  else:
    txt = "Foreground species are not monophylic."
    warnings.warn(txt)
    mark_foreground_p(marks, tree)

  # print(tree.write())

  tree.run_model("b_free")
  tree.run_model("M0")

  b_free = get_omega(tree, "b_free")
  p_value = tree.get_most_likely("b_free", "M0")
  p_value = '{:.7g}'.format(p_value)

  summary = list(map(str, [args.aln, p_value, b_free[" #1"], b_free[" #0"]]))
  print("\t".join(summary))


def get_foreground_species(aln, tsv, clade):
  df = pd.read_csv(tsv, sep = "\t")
  aln = fs.fasta2dict(aln)
  sp2cls = dict(df[["Organism Name", "class"]].values)
  sp2cls = {k.replace(" ", "_"): v for k, v in sp2cls.items()}
  fg = [sp for sp, cls in sp2cls.items() if sp in aln and cls in clade]
  return fg


def get_node(tree, node):
  res = tree.search_nodes(name=node)
  if len(res) > 1:
    exit('ERROR: more than 1 node with name: %s' % node)
  elif len(res) < 1:
    try:
      res = tree.search_nodes(node_id=int(node))
    except ValueError:
      exit('ERROR: node %s not found' % node)
    if len(res) < 1:
      exit('ERROR: node %s not found' % node)
  return res[0]


def mark_foreground_m(marks, tree):
  nodes = [get_node(tree, sp) for sp in marks]
  anc = tree.get_common_ancestor(nodes)
  node_ids = [node.node_id for node in anc.get_descendants() + [anc]]
  tree.mark_tree(node_ids, ["#1"]*len(node_ids))


def mark_foreground_p(marks, tree):
  nodes = [get_node(tree, sp) for sp in marks]
  anc = tree.get_common_ancestor(nodes)
  mark_nodes = []
  for node in nodes:
    while node.up:
      if node in anc.get_descendants()+[anc]:
        mark_nodes.append(node.node_id)
        node = node.up
      else:
        break
  tree.mark_tree(mark_nodes, ["#1"]*len(mark_nodes))


def get_omega(tree, model):
  mark2omega = {}
  result = tree.get_evol_model(model)
  for attr in result.branches.values():
    mark = attr.get('mark')
    omega = attr.get('w')
    mark2omega[mark] = omega
  return mark2omega


if __name__ == "__main__":
  main()
