#!/usr/bin/python
"""Given a list of biclique indices, merge highly overlapping bicliques, filter out-of-bounds rows and columns.

EXAMPLES
python $HOME/pymod/biclique_workflow/merge_biclique_indices.py fname=$HOME/net/gse15745_nov2/bicliques_pass6i3/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.630+_adj0.500_sup.biclique.rowmapped.0.380_mergetype1_graphmined.output.1380.532.remapout k=6 n_rows=24334 n_cols=10277 > $HOME/gse15745_nov2012_experiments/bicliques/biclique6.R

python $HOME/pymod/biclique_workflow/merge_biclique_indices.py fname=$HOME/net/gse15745_nov2/bicliques_pass7i4/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.560+_adj1.000_sup.biclique.rowmapped.0.380_mergetype1_graphmined.output k=7 n_rows=24334 n_cols=10277 > $HOME/gse15745_nov2012_experiments/bicliques/biclique7.R
"""
from __future__ import division
import sys

def main(fname, overlap=.50, k=1, n_rows=24334, n_cols=10277):
  rows, cols = [], []
  n_rows, n_cols = int(n_rows), int(n_cols)
  all_rows, all_cols = set(), set()
  for line in open(fname):
    r,c,d = [s.split(' ') for s in line.split(' ; ')]
    r = set(filter(lambda x: x<=n_rows, map(int, r)))
    c = set(filter(lambda x: x<=n_cols, map(int, c)))
    rows.append(r); cols.append(c)
    all_rows.update(r); all_cols.update(c)

  n = len(rows)
  overlaps, merged = set(), set()
  for i in xrange(n):
    if i in overlaps: continue
    merged.add(i)
    for j in xrange(i,n):
      if len(rows[i] & rows[j]) / len(rows[i]) >= overlap and \
            len(cols[i] & cols[j]) / len(cols[i]) >= overlap:
        rows[i] = rows[i] | rows[j]
        cols[i] = cols[i] | cols[j]
        overlaps.add(j)

  print "# Biclique Set %s Generated from %s" % (k, fname)
  print "#rows: ", ",".join(map(str, sorted(all_rows)))
  print "#cols: ", ",".join(map(str, sorted(all_cols)))
  print "# %d to %d bicliques, overlap by %f" % (n, len(merged), overlap)
  for i, x in enumerate(merged):
    print "bc%s.%d.row<-c(%s)" % (k, i, ",".join(map(str, sorted(rows[x]))))
    print "bc%s.%d.col<-c(%s)" % (k, i, ",".join(map(str, sorted(cols[x]))))
    
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
  
