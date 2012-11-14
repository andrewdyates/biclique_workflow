#!/usr/bin/python
from __future__ import division
import sys
"""Given a list of biclique indices, merge highly overlapping bicliques"""

def main(fname, overlap=.50):
  rows, cols = [], []
  for line in open(fname):
    r,c,d = [set(s.split(' ')) for s in line.split(' ; ')]
    rows.append(r); cols.append(c)

  n= len(rows)
  overlaps, merged = set(), set()
  for i in xrange(n):
    if i in overlaps: continue
    merged.add(i)
    for j in xrange(i,n):
      x = len(rows[i] & rows[j]) / len(rows[i])
      if len(rows[i] & rows[j]) / len(rows[i]) >= overlap and \
            len(cols[i] & cols[j]) / len(cols[i]) >= overlap:
        rows[i] = rows[i] | rows[j]
        cols[i] = cols[i] | cols[j]
        overlaps.add(j)

  print "# Generated from %s" % fname
  print "# %d to %d bicliques, overlap by %f" % (n, len(merged), overlap)
  for i, x in enumerate(merged):
    print "bc%d.row<-c(%s)" % (i, ",".join(sorted(rows[x])))
    print "bc%d.col<-c(%s)" % (i, ",".join(sorted(cols[x])))
    
if __name__ == "__main__":
  main(sys.argv[1])
  
