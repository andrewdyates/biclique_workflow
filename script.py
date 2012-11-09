#!/usr/bin/python
"""Mine bicliques from a rectangular dependency matrix.

EXAMPLE USE (in qsub script with 12 ppns and 1 node):
  
  python $HOME/biclique_workflow/script.py depM_fname=/fs/lustre/osu6683/gse15745_nov2/dependency_dispatch/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.pkl support=0.1 threshold=0.5 work_dir=/fs/lustre/osu6683/gse15745_nov2/bicliques thresh_cmp=greater
"""
from __future__ import division
from __init__ import *
import sys
import numpy as np
from matrix_io import *
import matrix_to_adjacency


def workflow(depM_fname=None, work_dir=None, support=10.0, threshold=0.6, thresh_cmp="greater", absvalue=False, add_mpiexec=False, overwrite=False, density=0.7, merge_type=1, extract_rows=None, extract_cols=None):
  threshold, support = float(threshold), float(support)
  if add_mpiexec in ('F', 'f', 'false', 'False', 'FALSE', 'None'):
    add_mpiexec = False
  density = float(density)
  merge_type = int(merge_type)
  assert density > 0 and density <= 1
  assert merge_type in (0,1)

  # 1) Create work directory
  if os.path.exists(work_dir):
    print "Work directory %s exists..." % (work_dir)
  else:
    make_dir(work_dir)
    print "Created Work directory %s." % (work_dir)

  # 2) Generate adjancency list and missing list
  print "Loading dependency matrix %s..." % (depM_fname)
  M = load(depM_fname)['M']
  print "Loaded matrix %s. Rows=%d, Cols=%d." % (os.path.basename(depM_fname), M.shape[0], M.shape[1])

  # 2.2) Handle row and column extractions
  row_filt_n, col_filt_n = 0, 0
  if extract_rows or extract_cols:
    print "Preshape:", M.shape
    if isinstance(extract_rows, basestring):
      extract_rows = extract_rows.split(',')
    if isinstance(extract_cols, basestring):
      extract_cols = extract_cols.split(',')
    if extract_rows:
      print "extract_rows(%d):" % len(extract_rows), extract_rows
      M = np.delete(M, extract_rows, 0)
      row_filt_n = len(extract_rows)
    if extract_cols:
      print "extract_cols(%d):" % len(extract_cols), extract_cols
      M = np.delete(M, extract_cols, 1)
      col_filt_n = len(extract_rows)
    print "Postshape:", M.shape
      
  empty_lines = []
  adj_iter = matrix_to_adjacency.npy_to_mafia(empty_lines, M=M, threshold=threshold, thresh_cmp=thresh_cmp, absvalue=absvalue)
  adj_fname, missing_fname = adj_list_fnames(depM_fname, work_dir, absvalue, thresh_cmp, threshold)

  if not (os.path.exists(adj_fname) and os.path.exists(missing_fname)) or overwrite:
    if os.path.exists(adj_fname) and os.path.exists(missing_fname):
      print "WARNING: %s and %s exist. Overwriting..." % (adj_fname, missing_fname)
    fp_adj = open(adj_fname, "w")
    fp_missing = open(missing_fname, "w")
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
  
    print "Wrote %d adjacency lines, %d missing lines. (row coverage %.4f%%)" % (n_lines, n_missing, n_lines/np.size(M,0)*100)
    print "%d of %d columns included. (col coverage %.4f%%)" % (len(col_set), np.size(M,1), len(col_set)/np.size(M,1)*100)
    print "SUM CHECK: Total: %d lines. Expected total %d lines. OK=%s" % (n_lines+n_missing, np.size(M,0), n_lines+n_missing == np.size(M,0))
    print "Wrote %d edges of %d all adj matrix entries. (%.6f%%)" % (n_edges, n_all_edges, edge_density*100)

    # average node degree statistics...
    avg_inc_row_degree = n_edges / n_lines
    max_avg_inc_row_degree = n_max_edges / n_lines
    avg_inc_col_degree = n_edges / len(col_set)
    max_avg_inc_col_degree = n_max_edges / len(col_set)
    avg_inc_degree = n_edges*2 / (len(col_set) + n_lines)
    max_avg_inc_degree = n_max_edges*2 / (len(col_set) + n_lines)
    
    print "Average (included) row degree: %.2f. Max: %.2f (%.4f%% of max)" % \
        (avg_inc_row_degree, max_avg_inc_row_degree, avg_inc_row_degree/max_avg_inc_row_degree*100)
    print "Average (included) col degree: %.2f. Max: %.2f (%.4f%% of max)" % \
        (avg_inc_col_degree, max_avg_inc_col_degree, avg_inc_col_degree/max_avg_inc_col_degree*100)
    print "Average (included) all node degree: %.2f. Max %.2f. (%.4f%% of max)" % (avg_inc_degree, max_avg_inc_degree, avg_inc_degree/max_avg_inc_degree*100)
    print "===================="
    print "Adjacency graph generation complete!"
    print "===================="
  else:
    print "adj_fname=%s and missing_fname=%s already exist and overwrite flag is false. Do not recompute adjacency matrix." \
        % (adj_fname, missing_fname)
  
  # 3) Call MBC on adjancency list
  biclique_fname = mbc_output_name(adj_fname, support)
  print "MBC biclique fname:", biclique_fname
  if not os.path.exists(biclique_fname) or overwrite:
    if os.path.exists(biclique_fname):
      print "WARNING: biclique_fname %s exists." % biclique_fname
    print "Calling MAFIA-MBC/MBC on %s with support %f%%..." % (adj_fname, support)
    #  (this will automatically call mafia, create freq_itemsets dir, and create fmi file.)
    call_mbc(adj_fname=adj_fname, support_percent=support, add_mpiexec=add_mpiexec)
  
  # 4) Density merge maximal bicliques
  # If a text-representation of npy_fname does not yet exist, create it

  bfname = os.path.basename(depM_fname).rpartition('.')[0]
  matrix_text_fname = os.path.join(work_dir, "%s.rfilt%d.cfilt%d.tab" % (bfname, len(row_filt_n), len(col_filt_n)))
  # SOMETHING HERE IS WRONG BECAUSE MATRIX IS BLANK
  if not os.path.exists(matrix_text_fname) or overwrite:
    if os.path.exists(matrix_text_fname):
      print "WARNING: %s exists. Overwriting..." % matrix_text_fname
    print "Saving dependency matrix as text to %s..." % matrix_text_fname
    save(M, open(matrix_text_fname, "w"), ftype="txt")
    print "Saved %s." % matrix_text_fname
  else:
    print "%s already exists. Do not resave." % (matrix_text_fname)
  # check for blank text matrix bug
  if open(matrix_text_fname).next().strip() == "":
    print "ERROR: %s is blank! Try to resave..." % matrix_text_fname 
    save(M, open(matrix_text_fname, "w"), ftype="txt")
    if open(matrix_text_fname).next().strip() == "":
      print "CRITICAL: %s is still blank after rewrite. Debug this problem before proceding." % matrix_text_fname
      sys.exit(1)

  # 4.1) Remap row numbers back to original matrix size (column numbers are never remapped)
  remapped_fname = preprocess_biclique_name(biclique_fname)
  if not os.path.exists(remapped_fname) or overwrite:
    s = write_preprocess_biclique(np.size(M,0), biclique_fname, missing_fname)
    assert s == remapped_fname
    print "Remapped MAFIA-MBC output row numbers to original matrix size in file %s == %s" % \
        (remapped_fname, preprocess_biclique_name(biclique_fname))
  else:
    print "%s exists, do not overwrite." % preprocess_biclique_name(biclique_fname)

  # 4.2) Call graphmining program
  merged_fname = graphmining_name(remapped_fname, density, merge_type)
  if not os.path.exists(merged_fname) or overwrite:
    print "Calling GraphMining on %s." % (remapped_fname)
    s = call_graphmining(graph_fname=matrix_text_fname, biclique_fname=remapped_fname, density=density, merge_type=merge_type)
    assert s == merged_fname
    print "Merged file saved to %s." % (merged_fname)
  else:
    print "GraphMining output file %s already exists. Do not overwrite." % (merged_fname)

  

  

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  workflow(**args)
  
