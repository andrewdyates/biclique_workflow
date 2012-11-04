#!/usr/bin/python
"""Mine bicliques from a rectangular dependency matrix.

EXAMPLE USE:
  python $HOME/biclique_workflow/test.py depM_fname=$HOME/gse15745/gse15745_gpl6104_gpl8490_spearman.values.npy support=20.0 threshold=0.6 work_dir=/fs/lustre/osu6683/gse15745_bicliques thresh_cmp=greater 
"""
from __future__ import division
from __init__ import *
import sys
import numpy as np
from matrix_io import *
import matrix_to_adjacency


def workflow(depM_fname=None, work_dir=None, support=10.0, threshold=0.6, thresh_cmp="greater", absvalue=False, overwrite=True):
  threshold, support = float(threshold), float(support)

  # 1) Create work directory
  if os.path.exists(work_dir):
    print "Work directory %s exists..." % (work_dir)
  else:
    make_dir(work_dir)
    print "Created Work directory %s." % (work_dir)

  # 2) Generate adjancency list and missing list
  print "Loading dependency matrix %s..." % (npy_fname)
  M = load(npy_fname)['M']
  empty_lines = []
  adj_iter = matrix_to_adjacency.npy_to_mafia(empty_lines, M=M, threshold=threshold, thresh_cmp=thresh_cmp, absvalue=absvalue)
  adj_fname, missing_fname = adj_list_fnames(depM_fname, work_dir, absvalue, thresh_cmp, threshold)
  
  if os.path.exists(adj_fname) and os.path.exists(missing_fname):
    print "WARNING: %s and %s exist. Overwriting..." % (adj_fname, missing_fname)
  fp_adj, fp_missing = open(adj_fname, "w"), open(missing_fname, "w")
  print "Opened files for writing:\nadj_matrix: %s\nmissing_row_list: %s" % (adj_fname, missing_fname)

  n_lines = 0
  n_edges = 0
  col_set = set() 
  for line in adj_iter:
    fp_adj.write(line); fp_adj.write("\n")
    n_lines += 1
    row = line.split(' ')
    n_edges += len(row)
    for x in row: col_set.add(x)
    fp_adj.close()
    fp_missing.write("\n".join(map(str, empty_lines)))
    fp_missing.close()
    n_missing = len(empty_lines)
    n_all_edges = np.size(M.ravel())
    edge_density = n_edges / n_all_edges
    n_max_edges = len(col_set) * n_lines

  print "Wrote %d adjacency lines, %d missing lines. (row coverage %.4f%%)" % (n_lines, n_missing, n_lines/np.size(M,1)*100)
  print "Total: %d lines. Expected total %d lines." % (n_lines+n_missing, np.size(M,0))
  print "Wrote %d edges of %d total edges. (%.6f%%)" % (n_edges, n_all_edges, edge_density*100)
  print "%d of %d columns included. (col coverage %.4f%%)" % (len(col_set), np.size(M,1), len(col_set)/np.size(M,1)*100)

  # average node degree statistics...
  avg_inc_row_degree = n_edges / n_lines
  max_avg_inc_row_degree = n_max_edges / n_lines
  avg_inc_col_degree = n_edges / len(col_set)
  max_avg_inc_col_degree = n_max_edges / len(col_set)
  avg_inc_degree = n_edges*2 / (len(col_set) + n_lines)
  max_avg_inc_degree = n_max_edges*2 / (len(col_set) + n_lines)
  
  print "Average (included) row degree: %.2f. Max: %.2f (%.4f%%)" % \
      (avg_inc_row_degree, max_avg_inc_row_degree, avg_inc_row_degree/max_avg_inc_row_degree*100)
  print "Average (included) col degree: %.2f. Max: %.2f (%.4f%%)" % \
      (avg_inc_col_degree, max_avg_inc_col_degree, avg_inc_col_degree/max_avg_inc_col_degree*100)
  print "Average (included) all node degree: %.2f. Max %.2f. (%.4f%%)" % (avg_inc_degree, max_avg_inc_degree, avg_inc_degree/max_avg_inc_degree*100)
  print "===================="
  print "Adjacency graph generation complete!"
  print "===================="
  
  # 3) Change working directory to MAFIA-MBC
  print "cd to %s" % PATHS["MAFIA-MBC"]
  os.chdir(PATHS["MAFIA-MBC"])
  
  # 4) Call MBC on adjancency list
  print "Calling MAFIA-MBC/MBC on %s with support %f%%..." % (adj_fname, support)
  #  (this will automatically call mafia, create freq_itemsets dir, and create fmi file.)
  biclique_fname = call_mbc(adj_fname=adj_fname, support_percent=support)
  print "MBC biclique fname:", biclique_fname
  print "to be continued..."
  # 6) Density merge maximal bicliques
  # 6.2) If a text-representation of npy_fname does not yet exist, create it
  #matrix_text_path = os.path.join(work_dir, os.path.basename(npy_fname)+".txt")
  #if not os.path.exists(matrix_text_path):
  #    npy_to_text(matrix_text_path, npy_to_text)
  
  # 7) Remap IDs to row and column numbers

  #s = call_mafia(in_fname=adj_file, support=support)
  #print s
  # TODO: DELETE freq_itemsets dir
  # verify that number of lines in adj + missing equals number of rows

  # print "==="
  # s = call_graph(density=density, graph_fname=matrix_file, clique_fname=biclique)
  # print s
  # print "==="

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  workflow(**args)
  
