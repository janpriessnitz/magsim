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

def fit_domain_wall(xs, magzs, p0=[300, 30, 1]):
    (fit_x0, fit_delta, fit_me), res = scipy.optimize.curve_fit(tanh_model, xs, magzs, p0=p0)
    return (fit_x0, fit_delta, fit_me)


# interface_width_fixup = 1.7011094153574284 - (413.963767 - 402.532814) # actual interface width vs interface width in ASD calculation
# def fix_interface_width(profile):
#     zmin = np.min(profile[:,0])
#     zmax = np.max(profile[:,0])
#     zmiddle = (zmax + zmin)/2
#     profile[profile[:,0] > zmiddle,0] += interface_width_fixup
#     return profile

from matplotlib.animation import FuncAnimation

output_dir = sys.argv[1]

prof = np.genfromtxt(os.path.join(output_dir, "profile.out0"))

fig, ax = plt.subplots()
ln, = ax.plot(prof[:,0], prof[:,3], linestyle='', marker='.')

fit_xs = prof[:,0]
ln_fit, = ax.plot(fit_xs, tanh_model(fit_xs, 240, 30, 1))
ax.set_ylim(-1, 1)

scale = 2.47e-1


widths = []

def init():
    return ln,

x0 = 300
delta = 30
ms = 1

def update(i):
    global x0, delta, ms
    prof = np.genfromtxt(os.path.join(output_dir, "profile.out{}".format(i)))
    # spins = np.genfromtxt(os.path.join(expdir, "lattice.out{}".format(i)))
    # xs = pos[:,0]
    # ys = spins[:,2]
    # bin_means, bin_edges, binnumber = scipy.stats.binned_statistic(xs, ys, bins=200)
    scale = 2.47e-1
    x0, delta, ms = fit_domain_wall(prof[:,0], prof[:,3], [x0, delta, ms])
    print("domain wall width [nm]: ", scale*delta)
    widths.append(scale*delta)
    ln.set_data(prof[:,0], prof[:,3])
    ln_fit.set_data(fit_xs, tanh_model(fit_xs, x0, delta, ms))

    return ln,

n = 200
if len(sys.argv) > 2:
    n = int(sys.argv[2])
ani = FuncAnimation(fig, update, frames=n,
                    init_func=init, blit=True)
plt.show()
ani.save(os.path.join(output_dir, "anim_profile.gif"))
np.savetxt(os.path.join(output_dir, "widths.txt"), widths)

# spins = np.genfromtxt(os.path.join(expdir, "lattice.out{}".format(n-1)))
# xs = pos[:,0]
# ys = spins[:,2]
# x0, delta, ms = fit_domain_wall(xs, ys)
# print("domain wall width [nm]: ", scale*delta)
