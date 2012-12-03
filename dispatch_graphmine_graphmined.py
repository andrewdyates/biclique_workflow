#!/usr/bin/python
"""
python $HOME/pymod/biclique_workflow/dispatch_graphmine_graphmined.py outdir=$HOME/density_merge_bicliques dry=True
"""
from __future__ import division
import qsub
import sys
import numpy as np
import merge_graphmine_graphmined

CMD = "python /nfs/01/osu6683/pymod/biclique_workflow/batch_graphmine_graphmined.py density=%(density)s merge=%(merge)s outdir=%(outdir)s in_fname=%(in_fname)s"

def runloop(outdir="", dry=False, dcor_filter=True, dcor_thresh=0.08):
  if isinstance(dcor_filter,basestring) and dcor_filter.lower() in ('f', 'none', 'false'):
    dcor_filter = False
  dcor_thresh = float(dcor_thresh)
  # generate input
  in_fname = merge_graphmine_graphmined.main(outdir=outdir, dcor_filter=dcor_filter, dcor_thresh=dcor_thresh)
  print "Generated %s." % in_fname
  for merge in ('0', '1'):
    for density in map(lambda x: "%.3f"%x, np.arange(10,91,2.5)/100):
      jobname = "bcm_%s_%s" % (merge, density)
      Q = qsub.Qsub(n_nodes=1, n_ppn=12, hours=3, work_dir=outdir, email=False, jobname=jobname)
      cmd = CMD % {'density': density, 'merge': merge, 'outdir': outdir, 'in_fname': in_fname}
      Q.add(cmd)
      print Q.script()
      Q.submit(dry)

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  runloop(**args)
