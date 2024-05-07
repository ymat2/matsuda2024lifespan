import os
import glob
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-d", "--dir")
  parser.add_argument("-o", "--outfile")
  args = parser.parse_args()

  dirs = glob.glob(args.dir+"/*")

  fg_lines = {}
  for d in dirs:
    if os.path.exists(d+"/csubst_cb_2.tsv"):
      res = d+"/csubst_cb_2.tsv"
      symbol = d.split("/")[-1]
      with open(res) as f:
        lns = f.readlines()
        if len(lns) == 1:
          print("No line found in",res)
          continue
        else:
          fg_lines[lns[0]] = "symbol"
          for l in lns[1:]:
            fg_lines[l] = symbol
      print(d)


  with open(args.outfile, "w") as f:
    for k,v in fg_lines.items():
      f.write(v+"\t"+k)


if __name__ == "__main__":
  main()
