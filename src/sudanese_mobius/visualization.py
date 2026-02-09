#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors 
from matplotlib.widgets import RadioButtons
from typing import Optional
from .curvature import CurvatureData

def plot_surface(
    c: np.ndarray,
    curv: CurvatureData,
    K_S3: Optional[np.ndarray] = None,
    scale: Optional[np.ndarray] = None,
    stride: int = 4,
    initial_key: str = "|H|"
):
    """
    Plots the surface with various color maps to demonstrate differences in
    curvature and metric between S3 and R3.

    The parameter `stride` decides the resolution of the surface. Setting 
    `stride=4` was the most stable for my computer but it can be lowered as 
    long as the resolution%stride = 1
    Scale is provided in both the curvature data and directly into the function
    so it can be set as optional.  Meaning if Scale is not provided the
    simulation will still run

    Parameters
    ----------
    c : ndarray
        Array of points (shape: (..., 3) representing the projection of the 
        surface in R3.
    curv : CurvatureData
        Data class containg
        -mean_curvature: ndarray
            Array of points (shape: (...,) representing the Mean Curvature of 
            the surface in R3.
        -gaussian_curvature : ndarray
            Array of points (shape: (...,) representing the Gaussian Curvature 
            of the surface in R3.
        -scale_factor: ndarray
            Array of scalar values at each point (shape: (...,) representing 
            the surface's scale factor in R3.
    K_S3 : ndarray, optional
        Array of points (shape: (...,) representing the Gaussian Curvature of 
        the surface in S3.
    scale : ndarray, optional
        Array of scalar values at each point (shape: (...,) representing 
        the surface's scale factor in R3.    
    initial_key : string
        Used to set the initial radio button value

    Returns
    -------

    fig : matplotlib.figure.Figure
        The figure which the surface is plotted on
    ax : mpl_toolkits.mplot3d.Axes3D
        The 3D axes
    radio : RadioButtons
        The radio button used to control which colormap is shown
    
    """
    data_lib_full = {
        "|H|": np.abs(curv.mean_curvature),
        "$K_{R^3}$": curv.gaussian_curvature,
        "λ": curv.scale_factor,
    }
    
    colorbar_titles = {
    "|H|": "Absolute Mean Curvature (|H|)",
    "$K_{R^3}$": "Gaussian Curvature ($K_{R^3}$)",
    "λ": "Scale Factor (λ)",
    "$K_{S^3} - K_{R^3}$": "Difference in Gaussian Curvature ($K_{S^3} - K_{R^3}$)"
    }
    
    if K_S3 is not None:
        data_lib_full["$K_{S^3} - K_{R^3}$"] = (K_S3 - curv.gaussian_curvature)
        
    if scale is not None:
        data_lib_full["λ"] = scale
        
    data_lib = {k: v[::stride, ::stride] for k, v in data_lib_full.items()}
    c_s = c[::stride, ::stride, 0], c[::stride, ::stride, 1], c[::stride, ::stride, 2]
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    plt.subplots_adjust(left=0.25)
    
    current_data = data_lib[initial_key]
    
    vmin, vmax = np.nanmin(current_data), np.nanmax(current_data)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.inferno
        
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(current_data)
    cbar = plt.colorbar(mappable, ax=ax, shrink=0.7)
    cbar.set_label(colorbar_titles[initial_key])

    def draw_surface(data2d, color_norm, cmap_obj):
        # facecolors must be (n,m,4)
        face_colors = cmap_obj(color_norm(data2d))
        X, Y, Z = c_s
        surf = ax.plot_surface(
            X, Y, Z,
            facecolors=face_colors,
            linewidth=0,
            antialiased=False,
            rstride=1, cstride=1, shade=False
        )
        return surf
    
    surf = draw_surface(current_data, norm, cmap)
    
    
    
    def update_plot(label):
        nonlocal surf, mappable
        new_data = data_lib[label]
        vmin, vmax = np.nanmin(new_data), np.nanmax(new_data)
        
        # Choose normalization & cmap based on data
        if vmin < 0 < vmax:
            new_norm = colors.TwoSlopeNorm(vmin=-max(abs(vmin), abs(vmax)), vcenter=0.0, vmax=max(abs(vmin), abs(vmax)))
            new_cmap = plt.cm.coolwarm
        else:
            new_norm = colors.Normalize(vmin=vmin, vmax=vmax)
            new_cmap = plt.cm.inferno
                
                
        # Remove old surface and redraw
        if surf is not None:
            try:
                surf.remove()
            except Exception:
                pass

        surf = draw_surface(new_data, new_norm, new_cmap)

        # Update colorbar
        mappable.set_norm(new_norm)
        mappable.set_cmap(new_cmap)
        mappable.set_array(new_data)
        cbar.update_normal(mappable)
        cbar.set_label(colorbar_titles[label])

        fig.canvas.draw_idle()

    
    rax = plt.axes([0.02, 0.4, 0.18, 0.25], facecolor="#EAEAF2")
    radio = RadioButtons(rax, list(data_lib.keys()), active=list(data_lib.keys()).index(initial_key))
    radio.on_clicked(update_plot)
    
    fig.radio = radio
    
    # Cosmetic: set aspect
    Xflat, Yflat, Zflat = c[:, :, 0], c[:, :, 1], c[:, :, 2]
    ax.set_box_aspect([np.ptp(Xflat), np.ptp(Yflat), np.ptp(Zflat)])
    ax.set_title("Visualization of Sudanese Mobius Strip")
    
    plt.style.use("seaborn-v0_8")
    plt.show()

    return fig, ax, radio

