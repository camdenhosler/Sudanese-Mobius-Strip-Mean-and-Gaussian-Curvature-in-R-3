#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project: Sudanese Mobius Strip Visualization
Author: Camden Hosler
Date: December 2025

Description:
This script calculates the mean curvature of a sudanese mobius strip projected
from S^3 into R^3 and visualizes it using Matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

res=200
u_range = np.linspace(-1/2 * np.pi, 1/2 * np.pi, res)
v_range = np.linspace(0, 2 * np.pi, res)

u, v = np.meshgrid(u_range, v_range)

du = u_range[1] - u_range[0]
dv = v_range[1] - v_range[0]

#Parameterization of 4D mobius strip
t = 0.5
x_1 = np.cos(u) * np.cos(v)
x_2 = np.cos(u) * np.sin(v)
x_3 = np.sin(u) * np.cos(t*v)
x_4 = np.sin(u) * np.sin(t*v)

#Rotated mobius strip off the singularity at the pole by a rotation matrix
x_n1 = 1/2*(x_1 - x_2 - x_3 - x_4)
x_n2 = 1/2*(x_1 + x_2 - x_3 + x_4)
x_n3 = 1/2*(x_1 + x_2 + x_3 - x_4)
x_n4 = 1/2*(x_1 - x_2 + x_3 + x_4)

#Stereographically projected strip
c_1 = x_n1/(1 - x_n4)
c_2 = x_n2/(1 - x_n4)
c_3 = x_n3/(1 - x_n4)

#Projection rotated back
c_n1 = c_1
c_n2 = np.sqrt(2)/2 * c_2 - np.sqrt(2)/2 * c_3
c_n3 = np.sqrt(2)/2 * c_2 + np.sqrt(2)/2 * c_3

#Mean Curvature Calculations in R^3
#axis 1 = u, axis 0 = v
xu = np.gradient(c_n1, du, axis=1)
xv = np.gradient(c_n1, dv, axis=0)

yu = np.gradient(c_n2, du, axis=1)
yv = np.gradient(c_n2, dv, axis=0)

zu = np.gradient(c_n3, du, axis=1)
zv = np.gradient(c_n3, dv, axis=0)

xuu = np.gradient (xu, du, axis=1)
xvv = np.gradient (xv, dv, axis=0)
xuv = 0.5*(
    np.gradient(xu, dv, axis=0) +
    np.gradient(xv, du, axis=1)
)


yuu = np.gradient (yu, du, axis=1)
yvv = np.gradient (yv, dv, axis=0)
yuv = 0.5*(
    np.gradient(yu, dv, axis=0) +
    np.gradient(yv, du, axis=1)
)

zuu = np.gradient (zu, du, axis=1)
zvv = np.gradient (zv, dv, axis=0)
zuv = 0.5*(
    np.gradient(zu, dv, axis=0) +
    np.gradient(zv, du, axis=1)
)

#Normal Vector
nx_org = yu*zv - zu*yv
ny_org = zu*xv - xu*zv
nz_org = xu*yv - yu*xv

norm_n = np.sqrt(nx_org**2 + ny_org**2 + nz_org**2)
norm_n[norm_n == 0] = 1e-12

nx = nx_org / norm_n
ny = ny_org / norm_n
nz = nz_org / norm_n

#First Fundamental Form Coefs
E = xu**2 + yu**2 + zu**2
F = xu*xv + yu*yv + zu*zv
G = xv**2 + yv**2 + zv**2

#Second Fundamental Form Coefs
L = xuu*nx + yuu*ny + zuu*nz
M = xuv*nx + yuv*ny + zuv*nz
N = xvv*nx + yvv*ny + zvv*nz

#Mean Curvature
denom = 2*(E*G - F**2)
denom[denom == 0] = 1e-12

H = (L*G - 2*M*F + N*E) / denom

K = (L*N - M**2) / (E*G - F**2)

#Color for Plot
Hmin = np.min(H)
Hmax = np.max(H)

fin_norm = colors.TwoSlopeNorm(vmin=Hmin, vcenter=0, vmax=Hmax)

raw_cmap = plt.colormaps['inferno']

fin_cmap = raw_cmap(fin_norm(H))

fin_mappable = plt.cm.ScalarMappable(norm=fin_norm, cmap=raw_cmap)

fin_mappable.set_array([])


#Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

plt.colorbar(fin_mappable, ax=ax, label="Mean Curvature")


surf = ax.plot_surface(c_n1, c_n2, c_n3, 
                       facecolors=fin_cmap, 
                       edgecolor='none', 
                       shade=False)

ax.set_box_aspect([np.ptp(c_n1),np.ptp(c_n2),np.ptp(c_n3)])
ax.set_title("Sudanese Mobius Strip Visualization")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

