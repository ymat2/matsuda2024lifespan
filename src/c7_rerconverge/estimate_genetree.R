library(RERconverge)

chr = commandArgs(trailingOnly=TRUE)[1]
ext = commandArgs(trailingOnly=TRUE)[2]
pattern = paste0("^", chr, ".*\\", ext)

# Estimate_phangorn_MLtrees ----------------------------------------------------
RERconverge::estimatePhangornTreeAll(
  alndir = "aln",
  pattern = pattern,
  treefile = "out/mastertree.nwk",
  output.file = paste0("aln/", chr, ".genetree.pep.trees"),
  submodel = "LG",
  type = "AA",
  format = "fasta"
)
