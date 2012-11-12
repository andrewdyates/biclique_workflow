from __init__ import *
from StringIO import *

BICLIQUE = """1 2 3 ; 1 2 3 21163 ; density: 0.762677
"""
extract_rows = [1,3,4] # 1 2 3 => 2 5 6
extract_cols = [2]


def main():
  in_fp = StringIO(BICLIQUE)
  out_fp = StringIO()
  remap_graphmining_out(in_fp, out_fp, 6, 6, extract_rows, extract_cols)
  print BICLIQUE
  print out_fp.getvalue()


if __name__ == "__main__":
  main()
