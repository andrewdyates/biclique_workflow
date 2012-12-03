#!/usr/bin/python
"""
Run GraphMining on merged outputs.
WARNING: This is a manual script written to work for a particular dataset.
in the future, I may extend it to be more general

EXAMPLE RUN:
/usr/bin/time python $HOME/pymod/biclique_workflow/batch_graphmine_graphmined.py in_fname=/nfs/01/osu6683/density_merge_bicliques/test/Graphmined.dcor_u0.080rowmapped density=0.45 outdir=$HOME/density_merge_bicliques/test merge=1 dry=True

EXAMPLE DENSITY MERGE RUN:
/usr/bin/time /nfs/01/osu6683/pythonbase/modules/biclique_workflow/MergeNetworkPatterns/SumNetwork/bin/GraphMining -g 2 -d 0.4000 -m 1 -p /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.800+_adj1.000_sup.biclique.rowmapped.0.400_mergetype1_graphmined.output /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.rfilt0.cfilt0.tab /fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values_0.800+_adj1.000_sup.biclique.rowmapped

SAMPLE GRAPHMERGE OUTPUT
56 190 277 312 334 ; 58 59 316 514 557 712 875 961 1109 1124 1134 1190 1266 1313 1426 1441 1494 1615 1763 1794 1828 2060 ; density: 0.762727

SAMPLE GRAPHMERGE INPUT
====================
20332 15237 23511 2439 ; 3587 3332 3466 4747 2070

"""
import os, subprocess, sys
import cPickle as pickle
import numpy as np

CMD = "/usr/bin/time %(exe)s -g 2 -d %(density)f -m %(merge)s -p %(out_fname)s %(tab_fname)s %(in_fname)s"

def main(in_fname=None, density=0.2, outdir="", merge='1', dry=False):
  assert in_fname
  density = float(density)
  merge = str(merge)
  out_fname = os.path.join(outdir, "%s.mmm_%.4f_%s.output" % \
                             (os.path.basename(in_fname), density, merge))
  params = {
      'density': density,
      'exe': "/nfs/01/osu6683/pythonbase/modules/biclique_workflow/MergeNetworkPatterns/SumNetwork/bin/GraphMining",
      'tab_fname': '/fs/lustre/osu6683/test_recompiled_mafia_mbc_mmm_.8/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.rfilt0.cfilt0.tab',
      'out_fname': out_fname,
      'in_fname': in_fname,
      'merge': merge
  }
  
  # Run density merge
  if not dry:
    cmd = CMD % params
    print cmd
    subprocess.call(cmd, shell=True)
  print "Complete!"

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)
