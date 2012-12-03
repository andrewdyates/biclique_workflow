#!/usr/bin/python
"""
Merge GraphMining outputs, run GraphMining on merged outputs.
WARNING: This is a manual script written to work for a particular dataset.
in the future, I may extend it to be more general



EXAMPLE DENSITY MERGE RUN:
/usr/bin/time /nfs/01/osu6683/pythonbase/modules/biclique_workflow/MergeNetworkPatterns/SumNetwork/bin/GraphMining -g 2 -d 0.4000 -m 1 -p /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.800+_adj1.000_sup.biclique.rowmapped.0.400_mergetype1_graphmined.output /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.rfilt0.cfilt0.tab /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.800+_adj1.000_sup.biclique.rowmapped

SAMPLE GRAPHMERGE OUTPUT
56 190 277 312 334 ; 58 59 316 514 557 712 875 961 1109 1124 1134 1190 1266 1313 1426 1441 1494 1615 1763 1794 1828 2060 ; density: 0.762727

SAMPLE GRAPHMERGE INPUT
====================
20332 15237 23511 2439 ; 3587 3332 3466 4747 2070

"""
import os, re, subprocess, sys
import cPickle as pickle

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

CMD = "/usr/bin/time %(exe)s -g 2 -d %(density)f -m %(merge)s -p %(out_fname)s %(tab_fname)s %(in_fname)s"

DCOR_M = "/fs/lustre/osu6683/gse15745_nov2/dependency_dispatch/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.pkl"

def filter_bad(D, rows, cols, t):
  """Filter rows and columns that don't significantly contribute.

  Args:
    D: np.matrix of Dcor dependency scores
    rows: [int] of indices, from 1
    cols: [int] of indices, from 1
    t: dcor mean threshold
  """
  r = np.array(row)-1
  c = np.array(col)-1
  DD = D[r,c]
  ru = np.mean(DD,0)
  cu = np.mean(DD,1)
  return (r[ru>t]+1, c[cu>t]+1)
  
def main(density=0.2, outdir="", merge='1', m=24334, n=10277, dcor_filter=False, dcor_thresh=0.08):
  dcor_thresh = float(dcor_thresh)
  if dcor_filter:
    print "Loading DCOR matrix..."
    DCOR = pickle.load(open(DCOR_M))
  m, n = int(m), int(n)
  density = float(density)
  merge = str(merge)
  in_fname = os.path.join(outdir, "Graphmined.%.4f.rowmapped" % (density))
  out_fname = os.path.join(outdir, "Graphmined.rowmapped.manual_mmm_%.4f_%s.output" % (density, merge))
  params = {
      'density': density,
      'exe': "/nfs/01/osu6683/pythonbase/modules/biclique_workflow/MergeNetworkPatterns/SumNetwork/bin/GraphMining",
      'tab_fname': '/fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.rfilt0.cfilt0.tab',
      'out_fname': out_fname,
      'in_fname': in_fname,
      'merge': merge
  }
  
  # combine .R files
  lines = []
  for k, fname in MERGEFILES.items():
    print "Combining bicliques from iteration %d..." % k
    for i, line in enumerate(open(fname)):
      rows, cols = map(lambda s: map(int, s.split(' ')), line.split(' ; ')[:2])
      rows = filter(lambda x: x<=m, rows)
      cols = filter(lambda x: x<=n, cols)
      if dcor_filter:
        m,n=len(rows),len(cols)
        rows, cols = filter_bad(DCOR, rows, cols, dcor_thresh)
        print "Filtered %d rows, %d columns." % (m-len(rows), n-len(cols))
      q = "%s ; %s" % (" ".join(map(str, rows)), " ".join(map(str, cols)))
      lines.append(q)
    print "  Added %d lines." % (i+1)
        
  # write combined .R files
  open(in_fname, "w").write("\n".join(lines))
  print "Wrote combined .R files in GraphMining input format to %s" % (in_fname)
  
  # Run density merge
  cmd = CMD % params
  print cmd
  subprocess.call(cmd, shell=True)
  print "Complete!"

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
