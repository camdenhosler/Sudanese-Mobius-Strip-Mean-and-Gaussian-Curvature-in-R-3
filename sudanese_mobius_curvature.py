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
    
def param_uv(u,v, TWIST_PARAMETER=0.5):
    #Parameterization of 4D mobius strip
    x_1 = np.cos(u) * np.cos(v)
    x_2 = np.cos(u) * np.sin(v)
    x_3 = np.sin(u) * np.cos(TWIST_PARAMETER*v)
    x_4 = np.sin(u) * np.sin(TWIST_PARAMETER*v)

    #Rotated mobius strip off the singularity at the pole by a rotation matrix
    x_n1 = 1/2*(x_1 - x_2 - x_3 - x_4)
    x_n2 = 1/2*(x_1 + x_2 - x_3 + x_4)
    x_n3 = 1/2*(x_1 + x_2 + x_3 - x_4)
    x_n4 = 1/2*(x_1 - x_2 + x_3 + x_4)
    
    return x_n1, x_n2, x_n3, x_n4

def stereoproj(x_n1, x_n2, x_n3, x_n4):
    #Stereographically projected strip
    c_1 = x_n1/(1 - x_n4)
    c_2 = x_n2/(1 - x_n4)
    c_3 = x_n3/(1 - x_n4)

    #Projection rotated back
    c_n1 = c_1
    c_n2 = np.sqrt(2)/2 * c_2 - np.sqrt(2)/2 * c_3
    c_n3 = np.sqrt(2)/2 * c_2 + np.sqrt(2)/2 * c_3
    
    return c_n1, c_n2, c_n3

def curvature_comp(c_n1, c_n2, c_n3, du, dv):
    #Mean Curvature Calculations in R^3
    c_n1_u, c_n1_v = np.gradient(c_n1, du, dv, edge_order=2)
    c_n2_u, c_n2_v = np.gradient(c_n2, du, dv, edge_order=2)
    c_n3_u, c_n3_v = np.gradient(c_n3, du, dv, edge_order=2)
    
    #c_n Notation is Shortened to just x y and z of Parameterization
    xu, xv = c_n1_u, c_n1_v 
    yu, yv = c_n2_u, c_n2_v
    zu, zv = c_n3_u, c_n3_v
    
    xuu, raw_xuv = np.gradient(xu, du, dv, edge_order=2)
    raw_xvu, xvv = np.gradient(xv, du, dv, edge_order=2)
    xuv = 0.5*(raw_xuv + raw_xvu)
    
    yuu, raw_yuv = np.gradient(yu, du, dv, edge_order=2)
    raw_yvu, yvv = np.gradient(yv, du, dv, edge_order=2)
    yuv = 0.5*(raw_yuv + raw_yvu)

    zuu, raw_zuv = np.gradient(zu, du, dv, edge_order=2)
    raw_zvu, zvv = np.gradient(zv, du, dv, edge_order=2)
    zuv = 0.5*(raw_zuv + raw_zvu)
    
    #Normal Vector
    nx_org = yu*zv - zu*yv
    ny_org = zu*xv - xu*zv
    nz_org = xu*yv - yu*xv
    
    norm_n = np.sqrt(nx_org**2 + ny_org**2 + nz_org**2)
    
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
    
    H = (L*G - 2*M*F + N*E) / denom
    
    K = (L*N - M**2) / (E*G - F**2)
    
    return norm_n, E, F, G, H, K

def singularity_check(norm_n, E, F, G, eps=1e-8):
    #Singularity Error (No Singularities Expected)
    if np.any(norm_n < eps) or np.any(E*G - F**2 < eps):
        raise ValueError("Unexpected near-zero value in curvature computation")

def plot_surface(norm_n, c_n1, c_n2, c_n3, H, K, eps=1e-8):
    H_masked = np.ma.masked_where(norm_n < eps, H)
    
    #Color for Plot
    Hmin = np.min(H)
    Hmax = np.max(H)
    Hmid = (Hmax + Hmin) / 2
    
    
    fin_norm = colors.TwoSlopeNorm(vmin=Hmin, vcenter=Hmid, vmax=Hmax)
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
    plt.show()

if __name__ == "__main__":
    
    RESOLUTION=200
    u_range = np.linspace(-1/2 * np.pi, 1/2 * np.pi, RESOLUTION)
    v_range = np.linspace(0, 2 * np.pi, RESOLUTION)

    u, v = np.meshgrid(u_range, v_range, indexing='ij')

    du = u_range[1] - u_range[0]
    dv = v_range[1] - v_range[0]
    
    x_n1, x_n2, x_n3, x_n4 = param_uv(u, v)
    c_n1, c_n2, c_n3 = stereoproj(x_n1, x_n2, x_n3, x_n4)
    norm_n, E, F, G, H, K = curvature_comp(c_n1, c_n2, c_n3, du, dv)
    singularity_check(norm_n, E, F, G)
    plot_surface(norm_n, c_n1, c_n2, c_n3, H, K)