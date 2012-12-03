import matplotlib
matplotlib.use('Agg')

def save_vis(D, fname, width=20, height=20):
  plt.clf(); plt.cla()
  fig = matplotlib.pyplot.figure(figsize=(width, height))
  q = fig.add_subplot(111)
  a = plt.imshow(D, axes=q)
  fig.savefig(fname, bbox_inches=0)

