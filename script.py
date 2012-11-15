#!/usr/bin/python
"""Mine bicliques from a rectangular dependency matrix.

EXAMPLE USE (in qsub script with 12 ppns and 1 node):
  
  python $HOME/biclique_workflow/script.py depM_fname=/fs/lustre/osu6683/gse15745_nov2/dependency_dispatch/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.pkl support=0.1 threshold=0.5 work_dir=/fs/lustre/osu6683/gse15745_nov2/bicliques thresh_cmp=greater

  mkdir /fs/lustre/osu6683/gse15745_nov2/bicliques_pass2
  python $HOME/pymod/biclique_workflow/script.py depM_fname=/fs/lustre/osu6683/gse15745_nov2/dependency_dispatch/Methyl-correct-aligned.pkl_mRNA-correct-aligned.pkl.DCOR.values.pkl support=1.5 threshold=0.72 work_dir=/fs/lustre/osu6683/gse15745_nov2/bicliques_pass2 thresh_cmp=greater density=0.60 extract_rows=34,56,73,190,277,304,312,334,393,526,627,783,1000,1045,1082,1144,1155,1162,1242,1263,1464,1467,1526,1562,1647,1659,1663,1734,1736,1892,2007,2063,2206,2258,2265,2271,2283,2354,2439,2482,2840,2857,2890,2944,3180,3199,3200,3207,3233,3244,3249,3461,3559,3595,3735,3808,3859,3926,3928,4012,4034,4243,4281,4402,4491,4497,4539,4567,4805,4941,5059,5202,5246,5281,5320,5388,5643,5702,5737,5817,5824,5849,5932,5935,5981,5993,6015,6021,6051,6063,6114,6124,6240,6269,6471,6482,6615,6815,6832,6889,6923,6933,7058,7103,7205,7289,7290,7552,7650,7750,7834,7887,7966,8008,8133,8190,8262,8327,8341,8351,8538,8612,8736,8749,8918,8968,9008,9020,9177,9211,9332,9442,9516,9699,9871,9968,10310,10321,10446,10528,10595,10600,10742,10815,10946,11018,11136,11170,11174,11201,11262,11266,11345,11349,11387,11465,11518,11597,11685,11847,11848,11862,11942,11962,12027,12060,12062,12067,12112,12214,12260,12343,12379,12398,12433,12505,12597,12622,12650,12762,12785,12818,13080,13111,13150,13254,13371,13483,13588,13607,13672,13761,13776,13800,13846,13854,13872,13877,13998,14080,14155,14201,14245,14263,14268,14276,14401,14477,14481,14487,14555,14585,14587,14592,14620,14675,14683,14708,14730,14804,14908,14932,14964,15010,15088,15200,15237,15247,15251,15500,15670,15692,15850,15853,15889,15947,15995,16000,16003,16014,16016,16023,16028,16032,16109,16112,16166,16187,16353,16361,16423,16456,16518,16527,16574,16585,16664,16763,16974,17225,17252,17292,17591,17644,17652,17750,17859,17913,17970,17989,18014,18078,18084,18103,18108,18125,18153,18163,18164,18199,18324,18341,18504,18519,18619,18661,18676,18722,18727,18734,18843,18894,18901,18988,19106,19292,19300,19325,19347,19391,19447,19476,19636,19643,19656,19752,19779,19813,19965,20019,20038,20063,20106,20146,20169,20189,20217,20228,20315,20332,20382,20407,20600,20721,20831,20840,20901,20940,20945,21149,21163,21182,21244,21327,21328,21437,21444,21474,21524,21540,21543,21560,21711,21793,21807,21829,21879,21906,21910,21930,21955,22113,22143,22271,22283,22402,22511,22512,22518,22572,22579,22607,22647,22711,22775,22795,22818,22866,22876,22957,23032,23100,23287,23308,23361,23431,23479,23511,23600,23668,23719,23757,23859,23989,24022,24028,24089,24126,24151,24183,24201 extract_cols=58,59,316,514,557,712,875,961,1109,1124,1134,1190,1266,1313,1426,1441,1494,1615,1711,1763,1794,1828,2060,2070,2093,2103,2216,2312,2431,2480,2700,2978,2997,3001,3064,3092,3118,3243,3263,3332,3436,3466,3568,3587,3684,3717,3915,4173,4333,4365,4366,4480,4487,4524,4601,4632,4690,4747,4927,4932,4961,5018,5157,5242,5268,5291,5320,5331,5388,5471,5661,5714,5880,6068,6096,6116,6147,6222,6228,6255,6293,6331,6379,6445,6461,6510,6646,6670,6780,6789,6797,7022,7024,7062,7210,7314,7617,7865,7879,8093,8094,8095,8106,8127,8140,8217,8218,8483,8557,8730,8784,8802,8841,8842,8998,9160,9326,9388,9399,9810,9876,10035,10037,10152,10193,10226,10261 overwrite=True >> /fs/lustre/osu6683/gse15745_nov2/bicliques_pass2
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
  log_fname = os.path.join(work_dir, 'log.txt')
  fp_log = open(log_fname, 'a')

  # 2) Generate adjancency list and missing list
  print "Loading dependency matrix %s..." % (depM_fname)
  M = load(depM_fname)['M']
  n_all_row, n_all_col = M.shape
  print "Loaded matrix %s. Rows=%d, Cols=%d." % (os.path.basename(depM_fname), M.shape[0], M.shape[1])

  # 2.2) Handle row and column extractions
  # WARNING: numpy indexes from 0, but MAFIA-MBC and R index from 1.
  # Only convert indices to 0-based for numpy matrix deletion!
  row_filt_n, col_filt_n = 0, 0
  if extract_rows or extract_cols:
    print "Preshape:", M.shape
    if isinstance(extract_rows, basestring):
      extract_rows = [int(s) for s in extract_rows.split(',')]
    if isinstance(extract_cols, basestring):
      extract_cols = [int(s) for s in extract_cols.split(',')]
    if extract_rows:
      print "extract_rows(%d):" % len(extract_rows), extract_rows
      M = np.delete(M, [x-1 for x in extract_rows], 0)
      row_filt_n = len(extract_rows)
    if extract_cols:
      print "extract_cols(%d):" % len(extract_cols), extract_cols
      M = np.delete(M, [x-1 for x in extract_cols], 1)
      col_filt_n = len(extract_cols)
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
  matrix_text_fname = os.path.join(work_dir, "%s.rfilt%d.cfilt%d.tab" % (bfname, row_filt_n, col_filt_n))
  print "Density merge matrix input: matrix_text_fname: %s" % matrix_text_fname
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

  # 5) Map back any extracted rows or columns
  if extract_rows or extract_cols:
    remap_fname = remapped_name(merged_fname, row_filt_n, col_filt_n)
    print "Remapping %s to %s by reinserting extracted rows and columns" % (merged_fname, remap_fname)
    remap_graphmining_out(open(merged_fname), open(remap_fname,'w'), n_all_row, n_all_col, extract_rows, extract_cols)
  print "Biclique mining complete."
  fp_log.close()
  

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  workflow(**args)
  
