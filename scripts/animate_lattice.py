#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os
import sys


def get_lattice(filename):
    with open(filename, "r") as f:
        w, h, E, T = f.readline().split(" ")
        w = int(w)
        h = int(h)
        E = float(E)
        T = float(T)
        lattice = []
        for i in range(w*h):
            x, y, z = f.readline().split(" ")
            x = float(x)
            y = float(y)
            z = float(z)
            lattice.append([x, y, z])
    return np.reshape(np.array(lattice), (w, h, 3)), E, T

def plot_animation(lattice_fnames, video_fname):
  lattice, E, T = get_lattice(lattice_fnames[0])
  fig, ax1 = plt.subplots(1, 1)
  im1 = ax1.imshow(lattice[:,:,2], vmax=1, vmin=-1)
  ax1.set_title("{}x{} lattice, E = {:.3f}/spin, T = {:.5f}".format(lattice.shape[0], lattice.shape[1], E/lattice.shape[0]/lattice.shape[1], T))
  fig.colorbar(im1)

  # animation function.  This is called sequentially
  def animate(i):
      lattice, E, T = get_lattice(lattice_fnames[i])
      im1.set_array(lattice[:,:,2])
      ax1.set_title("{}x{} lattice, E = {:.3f}/spin, T = {:.5f}".format(lattice.shape[0], lattice.shape[1], E/lattice.shape[0]/lattice.shape[1], T))

  anim = animation.FuncAnimation(fig, animate,
                                frames=len(lattice_fnames), interval=200)

  mywriter = animation.FFMpegWriter()
  anim.save(video_fname,writer=mywriter)

if __name__ == "__main__":
  sim_dir = sys.argv[1]
  video_fname = sys.argv[2]

  filenames = os.listdir(sim_dir)
  lattice_fnames = []
  for fname in filenames:
    if fname.startswith("lattice.dump"):
      i = int(fname[len("lattice.dump"):])
      lattice_fnames.append((i, sim_dir + "/" + fname))
  lattice_fnames_sorted = ['']*len(lattice_fnames)
  for i, fname in lattice_fnames:
    lattice_fnames_sorted[i] = fname
  plot_animation(lattice_fnames_sorted, video_fname)
