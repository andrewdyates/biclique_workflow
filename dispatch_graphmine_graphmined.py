#!/usr/bin/python
"""
python loop_merge_bicliques_course.py outdir=$HOME/density_merge_bicliques
"""
import qsub
import sys

CMD = "/usr/bin/time python /nfs/01/osu6683/density_merge_bicliques/merge_bicliques_course.py density=%(density)s merge=%(merge)s outdir=%(outdir)s"

def runloop(outdir="", dry=False):
  for merge in ('0', '1'):
    for density in map(str, (0.1, 0.2, 0.3, 0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9)):
      jobname = "bcm_%s_%s" % (merge, density)
      Q = qsub.Qsub(n_nodes=1, n_ppn=12, hours=20, work_dir=outdir, email=True, jobname=jobname)
      cmd = CMD % {'density': density, 'merge': merge, 'outdir': outdir}
      Q.add(cmd)
      print Q.script()
      Q.submit(dry)

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  runloop(**args)
