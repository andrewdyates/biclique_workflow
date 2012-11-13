"""
WARNING: TOO LOW OF SUPPORT WILL CRASH MBC!
"""
import numpy as np
import os, errno
from matrix_to_adjacency import *
from mafia_bicliques import *
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

def call_mbc(support_percent=5, adj_fname=None, add_time=True, add_mpiexec=False, verbose=True):
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
  curr_dir = os.getcwd()
  os.chdir(PATHS['MAFIA-MBC'])
  if verbose:
    print "Changed cwd from %s to %s" % (curr_dir, PATHS['MAFIA-MBC'])
  out_fname = mbc_output_name(adj_fname, support_percent)
  cmd = MBC_CMD % {
    'exe': PATHS['MBC'], 'support': support_percent, 'out_fname': out_fname, 'adj_fname': adj_fname}
  if add_time:
    cmd = "/usr/bin/time %s" % cmd
  if add_mpiexec:
    cmd = "mpiexec -npernode 1 %s" % (cmd)
  if verbose:
    print "MBC command to call via subprocess: %s" % cmd
  p = subprocess.call(cmd, shell=True)
  os.chdir(curr_dir)
  if verbose:
    print "---SUBPROCESS CALL OUTPUT---"
    print p
    print "Changed back to original directory %s" % (curr_dir)
  return out_fname

def mbc_output_name(adj_fname, support):
  return adj_fname.rpartition('.')[0] + "%.3f_sup.biclique" % (support)

def preprocess_biclique_name(biclique_fname):
  return biclique_fname + ".rowmapped"

def make_dir(outdir):
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise
  return outdir

def write_preprocess_biclique(n_rows=None, biclique_fname=None, missing_fname=None):
  n_rows = int(n_rows)
  assert n_rows > 0
  assert os.path.exists(biclique_fname)
  out_fname = preprocess_biclique_name(biclique_fname)
  if os.path.exists(out_fname):
    print "WARNING: target output file exists, will overwrite: %s" % (out_fname)

  missing_set = get_missing_set(open(missing_fname))
  print "Loaded %d missing row IDs from %s." % (len(missing_set), missing_fname)

  print "Loading %s..." % biclique_fname
  g = MafiaBiclique(fp=open(biclique_fname), missing_rows=missing_set, n=n_rows)
  print "Loaded raw MAFIA biclique from %s." % (biclique_fname)
  print "%d rows, %d cols covered in %d bicliques" % (len(g.row_set), len(g.col_set), len(g.bicliques))
  print "Row check: (%d + %d) == %d == %d?" % (len(g.row_set), len(g.missing_rows), \
                                    len(g.row_set)+len(g.missing_rows), n_rows)
  print "Writing remapped biclique for 'GraphMining' consumption at %s." % (out_fname)
  fp_out = open(out_fname, "w")
  n_lines_written = 0
  for line in g.yield_lines_to_density_merger():
    fp_out.write(line)
    n_lines_written += 1
  fp_out.close()
  print "Wrote %d lines to %s." % (n_lines_written, out_fname)
  return out_fname

def graphmining_name(biclique_fname, density, merge_type):
  return biclique_fname + ".%.3f_mergetype%d_graphmined.output" % (density, merge_type)

def remapped_name(graphmining_fname, nrows=0, ncols=0):
  return graphmining_fname + ".%d.%d.remapout" % (nrows, ncols)

def call_graphmining(density, graph_fname, biclique_fname, merge_type=1, add_time=True, add_mpiexec=False):
  density = float(density)
  merge_type = int(merge_type)
  assert os.path.exists(PATHS["GraphMining"])
  assert os.path.exists(biclique_fname)
  assert os.path.exists(graph_fname)
  assert merge_type in (0,1)

  out_fname = graphmining_name(biclique_fname, density, merge_type)
  if os.path.exists(out_fname):
    print "WARNING: target output file exists, will overwrite: %s" % out_fname

  assert density > 0 and density <= 1
  cmd = GRAPH_MINING_CMD % {
    'exe': PATHS["GraphMining"],
    'density': density,
    'out_fname': out_fname,
    'graph_fname': graph_fname,
    'biclique_fname': biclique_fname,
    'merge_type': merge_type,
    }
  if add_time:
    cmd = "/usr/bin/time %s" % cmd
  if add_mpiexec:
    cmd = "mpiexec -npernode 1 %s" % (cmd)
  print cmd
  p = subprocess.call(cmd, shell=True)
  print p
  return out_fname


def get_remap(n_all, missing_list, start=1):
  """Return dict of int to int index remapping."""
  if not missing_list:
    return None
  if not isinstance(missing_list[0], int):
    missing_list = [int(s) for s in missing_list]
  missing_list = sorted(missing_list)
  listmap = {}
  ii = 0
  for i in xrange(start, n_all-len(missing_list)+start):
    while ii < len(missing_list) and missing_list[ii] == i+ii:
      ii += 1
    listmap[i] = i + ii
  return listmap

def remap_graphmining_out(merged_fp, remap_fp, n_all_row, n_all_col, extract_rows, extract_cols, start=1):
  rowmap = get_remap(n_all_row, extract_rows)
  colmap = get_remap(n_all_col, extract_cols)
  for line in merged_fp:
    row_s, col_s, density_s = line.split(' ; ')
    row_idx = []
    for s in row_s.split(' '):
      i = int(s)
      if i-start >= n_all_row:
        break
      if rowmap:
        row_idx.append(str(rowmap[i]))
      else:
        row_idx.append(str(i))
    col_idx = []
    for s in col_s.split(' '):
      i = int(s)
      if i-start >= n_all_col:
        break
      if colmap:
        col_idx.append(str(colmap[i]))
      else:
        col_idx.append(str(i))
    remap_fp.write(' '.join(row_idx))
    remap_fp.write(' ; ')
    remap_fp.write(' '.join(col_idx))
    remap_fp.write(' ; ')
    remap_fp.write(density_s)
    
    

  
