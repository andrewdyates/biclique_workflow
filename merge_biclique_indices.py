#!/usr/bin/python
"""Given a list of biclique indices, merge highly overlapping bicliques, filter out-of-bounds rows and columns.

Merge percentage is:
For each pair:
- calcuate number of boxes in each biclique A, B = |rows| * |columns|
- calculate number of boxes in overlap between bicliques I = |rows| * |columns|
- A*pct >= I or B*pct >= I, merge
- repeat in loop until convergence

If given a dependency matrix, filter non-contributing rows and columns.

Requires HeapDict 1.0.0:
  http://pypi.python.org/pypi/HeapDict


EXAMPLES:

python $HOME/pymod/biclique_workflow/merge_biclique_indices.py fname=$HOME/net/gse15745_nov2/bicliques_pass10i2/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.400+_adj0.500_sup.biclique.rowmapped.0.280_mergetype1_graphmined.output.4131.2014.remapout k=10 n_rows=24334 n_cols=10277 > $HOME/gse15745_nov2012_experiments/bicliques/biclique10.R

python $HOME/pymod/biclique_workflow/merge_biclique_indices.py fname=~/density_merge_bicliques/Graphmined.rowmapped.manual_mmm_0.5000_1.output k=best n_rows=24334 n_cols=10277 overlap=0.65 > bicliques_t.5_o.65_areamerge.R

TODO: run for most recent bicliques
"""
from __future__ import division
# Install from http://pypi.python.org/pypi/HeapDict
from heapdict import heapdict
import sys
import cPickle as pickle
import merge_graphmine_graphmined

def main(fname, overlap=.60, k=1, n_rows=24334, n_cols=10277, d_fname=None, dcor_thresh=0.08):
  if d_fname:
    D = pickle.load(open(d_fname))
  dcor_thresh = float(dcor_thresh)
  rows, cols = [], []
  overlap = float(overlap); assert overlap >= 0 and overlap <= 1
  n_rows, n_cols = int(n_rows), int(n_cols)
  all_rows, all_cols = set(), set()
  for i, line in enumerate(open(fname)):
    r,c,d = (s.split(' ') for s in line.split(' ; '))
    r = filter(lambda x:  x<=n_rows, map(int, r)); r=set(r)
    c = filter(lambda x: x<=n_cols, map(int, c)); c=set(c)
    # Filter non-contributing rows and columns
    if d_fname:
      s,t = len(r),len(c)
      rf, cf = merge_graphmine_graphmined.filter_bad(D, r, c, dcor_thresh)
      r, c = set(rf), set(cf)
      sys.stderr.write("Filtered %d rows, %d columns from unmerged biclique %d.\n" % \
          (s-len(r), t-len(c), i))
    rows.append(r); cols.append(c)
    all_rows.update(r); all_cols.update(c)
  assert len(rows) == len(cols)
  bicliques = zip(rows, cols)
  merged = merge_bicliques(bicliques, overlap)

  print "# Biclique Set %s Generated from %s" % (k, fname)
  print "# %d to %d bicliques, overlap by %f" % (len(bicliques), len(merged), overlap)
  print "# %d unique rows, %d unique columns" % (len(all_rows), len(all_cols))
  print "#rows: ", ",".join(map(str, sorted(all_rows)))
  print "#cols: ", ",".join(map(str, sorted(all_cols)))
  for i, bc in enumerate(merged):
    print "bc%s.%d.row<-c(%s)" % (k, i, ",".join(map(str, sorted(bc[0]))))
    print "bc%s.%d.col<-c(%s)" % (k, i, ",".join(map(str, sorted(bc[1]))))
  return merged


def merge_bicliques(bicliques, max_overlap=.55):
  """Greedily merge biclique set B by maximum item overlap.
  TODO: add minimum merged density constraint.

  Args:
    bicliques: [(row, col)], row=[int], col=[int] of bicliques to merge
    max_overlap: float 0 <= x <= 1 of maximum allowable overlap of items without merge
  Returns:
    {}
  """
  # Calculate all pairs overlap
  # a symmetric distance matrix that I can update
  # a max heap to keep track of most overlapping biclique pair

  B = dict(enumerate(bicliques)) # uniquely index bicliques, but allow constant access
  D = heapdict()
  for i in B:
    for j in filter(lambda j:j>i, B):
      # higher percentage overlap is lower "priority"
      D[tuple(sorted((i,j)))] = -get_overlap(B[i], B[j])

  # Merge top pair while it is above threshold
  pair, overlap = D.popitem()
  while -overlap > max_overlap:
    i, j = sorted(pair)
    B[i] = (B[i][0]|B[j][0], B[i][1]|B[j][1])
    del B[j]
    # Update all distances including i, delete those with j
    for k in filter(lambda k:k!=i and k!=j, B):
      D[tuple(sorted((i,k)))] = -get_overlap(B[i], B[k])  # distance to merged biclique i
      del D[tuple(sorted((j,k)))]  # biclique j no longer exists
    if len(D):
      pair, overlap = D.popitem()  # get next pair
    else:
      break
  return B.values()

def get_overlap(a, b):
  """Compute overlap between bicliques a and b; biclique = ([rows],[cols])."""
  area_i = len(a[0])*len(a[1])
  area_j = len(b[0])*len(b[1])
  area_ij = len(a[0]&b[0])*len(a[1]&b[1])
  return max(area_ij/area_i, area_ij/area_j)
    
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
  
