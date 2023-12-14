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
    (fit_x0, fit_delta, fit_me), res = scipy.optimize.curve_fit(tanh_model, xs, magzs, p0=[80, 30, 1])
    return (fit_x0, fit_delta, fit_me)


rundir = sys.argv[1]

Ts = []
x0s = []
deltas = []
mes = []

measuredTs = [0.01, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200]
measuredTs = [str(T) for T in measuredTs]

# for tempdir in os.listdir(rundir):
for tempdir in measuredTs:
  try:
    T = float(tempdir)
    prof = np.genfromtxt(os.path.join(rundir, tempdir, "averege_profile.out"))
    x0, delta, me = fit_domain_wall(prof[:,0], prof[:,1])
  except:
     continue
  Ts.append(T)
  x0s.append(x0)
  deltas.append(delta)
  mes.append(me)

fig, (ax1, ax2) = plt.subplots(1, 2, dpi=250, figsize=(9, 4))
ax1.scatter(Ts, deltas)
ax2.scatter(Ts, mes)

ax1.set_xlabel("T [K]")
ax1.set_ylabel(r"$\delta_{dw}$ [nm]")
ax2.set_xlabel("T [K]")
ax2.set_ylabel(r"$M_s/|M|$")
plt.savefig(os.path.join(rundir, "dw_fit.png"))

np.savetxt(os.path.join(rundir, "dw_fit.out"), np.array([Ts, x0s, deltas, mes]).T, header="T, x_0, delta, m_e")