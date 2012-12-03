#!/usr/bin/python
"""
Merge GraphMining outputs into single GraphMining input file.
WARNING: This is a manual script written to work for a particular dataset.
in the future, I may extend it to be more general

EXAMPLE RUN:
/usr/bin/time python $HOME/pymod/biclique_workflow/merge_graphmine_graphmined.py outdir=$HOME/density_merge_bicliques/test dcor_filter=True dcor_thresh=0.08

SAMPLE GRAPHMERGE OUTPUT
56 190 277 312 334 ; 58 59 316 514 557 712 875 961 1109 1124 1134 1190 1266 1313 1426 1441 1494 1615 1763 1794 1828 2060 ; density: 0.762727

SAMPLE GRAPHMERGE INPUT
====================
20332 15237 23511 2439 ; 3587 3332 3466 4747 2070

"""
import os, subprocess, sys
import cPickle as pickle
import numpy as np

MERGEFILES = {
1: "/nfs/01/osu6683/net/gse15745_nov2/bicliques/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.740+_adj2.000_sup.biclique.rowmapped.0.650_mergetype1_graphmined.output",
2: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass2/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.750+_adj1.000_sup.biclique.rowmapped.0.600_mergetype1_graphmined.output",
3: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass3i11/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.720+_adj1.000_sup.biclique.rowmapped.0.550_mergetype1_graphmined.output.528.258.remapout",
4: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass4i5/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.690+_adj0.500_sup.biclique.rowmapped.0.500_mergetype1_graphmined.output.681.288.remapout",
5: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass5i3/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.670+_adj0.500_sup.biclique.rowmapped.0.480_mergetype1_graphmined.output.1081.378.remapout",
6: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass6i3/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.630+_adj0.500_sup.biclique.rowmapped.0.380_mergetype1_graphmined.output.1380.532.remapout",
7: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass7i5/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.580+_adj0.500_sup.biclique.rowmapped.0.400_mergetype1_graphmined.output.1981.803.remapout",
8: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass8i3/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.530+_adj0.500_sup.biclique.rowmapped.0.400_mergetype1_graphmined.output.2542.1100.remapout",
9: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass9i4/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.460+_adj0.500_sup.biclique.rowmapped.0.380_mergetype1_graphmined.output.3140.1481.remapout",
10: "/nfs/01/osu6683/net/gse15745_nov2/bicliques_pass10i2/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.400+_adj0.500_sup.biclique.rowmapped.0.280_mergetype1_graphmined.output.4131.2014.remapout",
}

DCOR_M = "/fs/lustre/osu6683/gse15745_nov2/dependency_dispatch/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.pkl"

def filter_bad(D, rows, cols, t):
  """Filter rows and columns that don't significantly contribute.

  Args:
    D: np.matrix of Dcor dependency scores
    rows: [int] of indices, from 1
    cols: [int] of indices, from 1
    t: dcor mean threshold
  """
  r = np.array(rows)-1
  c = np.array(cols)-1
  try:
    DD = D[r,:][:,c]
  except Exception:
    print rows, r
    print cols, c
    raise
  ru = np.mean(DD,1)
  cu = np.mean(DD,0)
  return (r[ru>t]+1, c[cu>t]+1)
  
def main(outdir="", m=24334, n=10277, dcor_filter=False, dcor_thresh=0.08):
  """Return file name of file written."""
  dcor_thresh = float(dcor_thresh)
  if dcor_filter:
    print "Loading DCOR matrix..."
    DCOR = pickle.load(open(DCOR_M))
  m, n = int(m), int(n)
  if dcor_filter:
    flag = "dcor_u%.3f" % dcor_thresh
  else:
    flag = ""
  in_fname = os.path.join(outdir, "Graphmined.%srowmapped" % (flag))
  
  # combine .R files
  lines = []
  for k, fname in MERGEFILES.items():
    print "Combining bicliques from iteration %d..." % k
    for i, line in enumerate(open(fname)):
      rows, cols = map(lambda s: map(int, s.split(' ')), line.split(' ; ')[:2])
      rows = filter(lambda x: x<=m, rows)
      cols = filter(lambda x: x<=n, cols)
      if dcor_filter:
        q,r=len(rows),len(cols)
        rows, cols = filter_bad(DCOR, rows, cols, dcor_thresh)
        print "Filtered %d rows, %d columns." % (q-len(rows), r-len(cols))
        if not len(rows) or not len(cols):
          print "All rows or columns filtered! Skip this biclique."
          continue
      q = "%s ; %s" % (" ".join(map(str, rows)), " ".join(map(str, cols)))
      lines.append(q)
    print "  Added %d lines." % (i+1)
        
  # write combined .R files
  open(in_fname, "w").write("\n".join(lines))
  print "Wrote combined GraphMining outputs to input format at %s" % (in_fname)
  return in_fname

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
