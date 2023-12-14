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


interface_width_fixup = 1.7011094153574284 - (413.963767 - 402.532814) # actual interface width vs interface width in ASD calculation for 1Ru
def fix_interface_width(profile):
    zmin = np.min(profile[:,0])
    zmax = np.max(profile[:,0])
    zmiddle = (zmax + zmin)/2
    profile[profile[:,0] > zmiddle,0] += interface_width_fixup
    return profile

output_dir = sys.argv[1]

print("averaging", output_dir)

n = 200
startN = 50
scale = 2.47e-1
prof = np.genfromtxt(os.path.join(output_dir, "profile.out{}".format(startN)))
prof = fix_interface_width(prof)
prof[:,0] *= scale

avgprof = prof
nprofs = 1

for i in range(startN+1, n):
  try:
    prof = np.genfromtxt(os.path.join(output_dir, "profile.out{}".format(i)))
    prof = fix_interface_width(prof)
    avgprof[:,1] += prof[:,1]
    nprofs += 1
  except:
    print("warning: profile not found:", os.path.join(output_dir, "profile.out{}".format(i)))
    break

avgprof[:,1] /= nprofs

np.savetxt(os.path.join(output_dir, "../", "average_profile.out"), avgprof)

fig, ax = plt.subplots()
ln, = ax.plot(avgprof[:,0], avgprof[:,1], linestyle='', marker='.')

ax.set_ylim(-1, 1)
plt.savefig(os.path.join(output_dir, "..", "average_profile.png"))
