#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.stats
import os
import sys

def tanh_model(xs, x0, delta, me):
    # return np.abs(me*np.tanh(np.pi*(xs-x0)/delta))
    return me*np.tanh(np.pi*(xs-x0)/delta)

def fit_domain_wall(xs, magzs):
    (fit_x0, fit_delta, fit_me), res = scipy.optimize.curve_fit(tanh_model, xs, magzs, p0=[240, 30, 1])
    return (fit_x0, fit_delta, fit_me)


interface_width_fixup = 1.7011094153574284 - (413.963767 - 402.532814) # actual interface width vs interface width in ASD calculation
def fix_interface_width(profile):
    zmin = np.min(profile[:,0])
    zmax = np.max(profile[:,0])
    zmiddle = (zmax + zmin)/2
    profile[profile[:,0] > zmiddle,0] += interface_width_fixup
    return profile

from matplotlib.animation import FuncAnimation

# output_dir = 'scratch/CoRuCo/temp_sweep/0K_nowall/output'
# output_dir = 'scratch/CoRuCo/temp_sweep/Ru1_0K/output/'
# output_dir = 'scratch/CoRuCo/temp_sweep/no3_nowall/1000/output'
output_dir = sys.argv[1]

prof = np.genfromtxt(os.path.join(output_dir, "profile.out0"))
prof = fix_interface_width(prof)

fig, ax = plt.subplots()
ln, = ax.plot(prof[:,0], prof[:,1], linestyle='', marker='.')

fit_xs = np.array(sorted(prof[:,0]))
ln_fit, = ax.plot(fit_xs, tanh_model(fit_xs, 240, 30, 1))
ax.set_ylim(-1, 1)

scale = 2.47e-1


widths = []

def init():
    return ln,

def update(i):
    prof = np.genfromtxt(os.path.join(output_dir, "profile.out{}".format(i)))
    prof = fix_interface_width(prof)
    # spins = np.genfromtxt(os.path.join(expdir, "lattice.out{}".format(i)))
    # xs = pos[:,0]
    # ys = spins[:,2]
    # bin_means, bin_edges, binnumber = scipy.stats.binned_statistic(xs, ys, bins=200)
    # scale = 2.47e-1
    x0, delta, ms = fit_domain_wall(prof[:,0], prof[:,1])
    # print("domain wall width [nm]: ", scale*delta)
    widths.append(scale*delta)
    ln.set_data(prof[:,0], prof[:,1])

    ln_fit.set_data(fit_xs, tanh_model(fit_xs, x0, delta, ms))

    return ln,

n = 100
ani = FuncAnimation(fig, update, frames=n,
                    init_func=init, blit=True)
plt.show()
ani.save(os.path.join(output_dir, "profile.gif"))
np.savetxt(os.path.join(output_dir, "widths.txt"), widths)

# spins = np.genfromtxt(os.path.join(expdir, "lattice.out{}".format(n-1)))
# xs = pos[:,0]
# ys = spins[:,2]
# x0, delta, ms = fit_domain_wall(xs, ys)
# print("domain wall width [nm]: ", scale*delta)
