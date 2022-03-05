#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
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

def plot_image(lattice, E, T, image_name):
  plt.imshow(lattice[:,:,2], vmin=-1, vmax=1)
  plt.colorbar()
  plt.title("{}x{} lattice, E = {:.3f}/spin, T = {:.5f}".format(lattice.shape[0], lattice.shape[1], E/lattice.shape[0]/lattice.shape[1], T))
  plt.savefig(image_name)

def plot_lattice(lattice_filename, image_filename):
  lattice, E, T = get_lattice(lattice_filename)
  plot_image(lattice, E, T, image_filename)

if __name__ == "__main__":
  plot_lattice(sys.argv[1], sys.argv[2])