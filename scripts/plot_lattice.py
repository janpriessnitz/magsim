#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import os
import sys


def get_lattice_dump(filename):
    lattices = []
    Es = []
    Ts = []
    steps_list = []
    with open(filename, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            w, h, E, T, steps = header.split(" ")
            w = int(w)
            h = int(h)
            E = float(E)
            T = float(T)
            Es.append(E)
            Ts.append(T)
            steps_list.append(int(steps))
            for i in range(w*h):
                x, y, z = f.readline().split(" ")
                x = float(x)
                y = float(y)
                z = float(z)
                lattices.append([x, y, z])
    return np.reshape(np.array(lattices), (-1, w, h, 3)), Es, Ts, steps_list

def plot_image(lattice, E, T, image_name):
  plt.imshow(lattice[:,:,2], vmin=-1, vmax=1)
  plt.colorbar()
  N = lattice.shape[0]*lattice.shape[1]
  plt.title("E/cell={}".format(E/N*2))

  plt.title("{}x{} lattice, E = {:.3f}/spin, T = {:.5f}".format(lattice.shape[0], lattice.shape[1], E/lattice.shape[0]/lattice.shape[1], T))
  plt.savefig(image_name)

def plot_lattice(lattice_filename, image_filename):
  lattices, Es, Ts, steps = get_lattice(lattice_filename)
  plot_image(lattice, E, T, image_filename)

if __name__ == "__main__":
  plot_lattice(sys.argv[1], sys.argv[2])