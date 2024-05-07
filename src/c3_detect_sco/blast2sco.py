import argparse
import glob


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-d", "--directory", type=str, help="directory where .blastp are")
  parser.add_argument("-r", "--reference", type=str, help="reference species")
  parser.add_argument("-m", "--min_sp", type=int, help="minimum number of species to have SCO", default=3)
  parser.add_argument("-o", "--output", help="output filename")
  args = parser.parse_args()

  BLAST_QUERY = glob.glob(args.directory+"/*_to_ref.blastp")
  ACCESSION = [i.split("/")[-1].removesuffix("_to_ref.blastp") for i in BLAST_QUERY]

  sco = {}
  for sp in ACCESSION:
    if sp != args.reference:
      qry_blast_file = args.directory+"/"+sp+"_to_ref.blastp"
      db_blast_file = args.directory+"/ref_to_"+sp+".blastp"
      qry2ref = blast2dict(qry_blast_file)
      ref2db = blast2dict(db_blast_file)
      for ref_gene in ref2db:
        if ref_gene not in sco:
          sco[ref_gene] = []
        forward_hit = ref2db[ref_gene]
        reverse_hit = qry2ref.get(forward_hit)
        if ref_gene == reverse_hit:
          sco[ref_gene].append(forward_hit)
    else:
      continue


  with open(args.output, "w") as f:
    for gn in sco:
      if len(sco[gn]) >= args.min_sp:
        f.write(gn+"\t"+"\t".join(sco[gn])+"\n")


def blast2dict(pth):
  dct = {}
  with open(pth) as f:
    for line in f:
      cols = line.rstrip("\n").split("\t")
      qry, db = cols[0], cols[1]
      if qry not in dct:
        dct[qry] = db
  return dct


def parse_tsv(pth):
  dct = {}
  with open(pth) as f:
    next(f)
    for line in f:
      tabs = line.split("\t")
      sp_accession, sp_class = tabs[0], tabs[7]
      dct[sp_accession] = sp_class
  return dct

if __name__ == "__main__":
  main()
