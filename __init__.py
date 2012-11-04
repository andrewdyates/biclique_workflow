"""Execute

WARNING: TOO LOW OF SUPPORT WILL CRASH MBC!
"""
import numpy as np
import os, errno
from matrix_to_adjacency import *
#from mafia_bicliques import *
import subprocess


# add time and mpiexec prefixes in CMD call
# e.g., mpiexec -npernode %d time %s
MAFIA_CMD = "%(exe)s -fci %(support).4f -ascii %(in_fname)s %(out_fname)s"
MBC_CMD = "%(exe)s fci %(adj_fname)s %(support).2f %(out_fname)s 0"
GRAPH_MINING_CMD = "%(exe)s -g 2 -d %(density).4f -m %(merge_type)d -p %(out_fname)s %(graph_fname)s %(biclique_fname)s"

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PATHS = {
  "GraphMining": os.path.join(THIS_DIR, "MergeNetPackage/SumNetwork/bin/GraphMining"),
  "GraphVis": os.path.join(THIS_DIR, "MergeNetPackage/SumNetwork/bin/GraphVis"),
  "mafia": os.path.join(THIS_DIR, "MergeNetPackage/MAFIA-MBC/bin/mafia"),
  "MBC": os.path.join(THIS_DIR, "MergeNetPackage/MAFIA-MBC/bin/MBC"),
  "MAFIA-MBC": os.path.join(THIS_DIR, "MergeNetPackage/MAFIA-MBC"),
}

def adj_list_fnames(npy_fname, work_dir, absvalue, thresh_cmp, threshold):
  if absvalue:
    absvalue = "abs"
  else:
    absvalue=""
  if thresh_cmp == "greater":
    thresh_cmp = "+"
  elif thresh_cmp == "less":
    thresh_cmp = "-"
  npy_basename = os.path.basename(npy_fname)
  prefix = "%s_%.3f%s%s_adj" % (npy_basename.rpartition('.')[0], threshold, absvalue, thresh_cmp)
  adj_fname = os.path.join(work_dir, prefix+".dat")
  missing_fname = os.path.join(work_dir, prefix+".missing")
  return adj_fname, missing_fname

def call_mbc(support_percent=5, adj_fname=None, add_time=True, add_mpiexec=True, npernode=12, verbose=True):
  """Note: requires .mafia file in ./freq_itemsets/ dir. Support in %, not decimal.
  Calling this without an expected .mafia file will generate it automatically.
  
  Args:
    support_percent: float support percent, e.g., 5 = 5%
    adj_fname: str fpath to adjacency file name of input
  Returns:
    str of fname to MAFIA-MBC output
  """
  support_percent = float(support_percent)
  assert support_percent > 0 and support_percent <= 100 
  # change to MAFIA path
  os.chdir(PATHS['MAFIA-MBC'])
  out_fname = mbc_output_name(adj_fname, support_percent)
  cmd = MBC_CMD % {
    'exe': PATHS['MBC'], 'support': support_percent, 'out_fname': out_fname, 'adj_fname': adj_fname}
  if add_time:
    cmd = "time %s" % cmd
  if add_mpiexec:
    cmd = "mpiexec -npernode %d %s" % (npernode, cmd)
  if verbose:
    print "MBC command to call via subprocess: %s" % cmd
  p = subprocess.call(cmd, shell=True)
  if verbose:
    print "---SUBPROCESS CALL OUTPUT---"
    print p
  return out_fname

def mbc_output_name(adj_fname, support):
  return adj_fname.rpartition('.')[0] + "%.3f_sup.biclique" % (support)

def make_dir(outdir):
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise
  return outdir
