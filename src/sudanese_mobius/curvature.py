#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from dataclasses import dataclass

@dataclass
class CurvatureData:
    mean_curvature: np.ndarray
    gaussian_curvature: np.ndarray
    scale_factor: np.ndarray


def gaussian_curvature_S3(x, du, dv):
    #Gaussian Curvature in S^3
    xu = np.gradient(x, du, axis=0, edge_order=2)
    xv = np.gradient(x, dv, axis=1, edge_order=2)

    xuu = np.gradient(xu, du, axis=0, edge_order=2)
    xuv = np.gradient(xu, dv, axis=1, edge_order=2)
    xvv = np.gradient(xv, dv, axis=1, edge_order=2)

    def det3(a, b, c):
        return (
            a[...,0]*b[...,1]*c[...,2] +
            a[...,1]*b[...,2]*c[...,0] +
            a[...,2]*b[...,0]*c[...,1] -
            a[...,2]*b[...,1]*c[...,0] -
            a[...,1]*b[...,0]*c[...,2] -
            a[...,0]*b[...,2]*c[...,1]
        )

    #Normal Vector
    n = np.stack([
        det3(x[...,1:], xu[...,1:], xv[...,1:]),
       -det3(x[..., [0,2,3]], xu[..., [0,2,3]], xv[..., [0,2,3]]),
        det3(x[..., [0,1,3]], xu[..., [0,1,3]], xv[..., [0,1,3]]),
       -det3(x[...,0:3], xu[...,0:3], xv[...,0:3])
    ], axis=-1)

    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    
    #First Fundamental Form Coefs
    E = np.sum(xu * xu, axis=-1)
    F = np.sum(xu * xv, axis=-1)
    G = np.sum(xv * xv, axis=-1)
    denom = E * G - F**2
    
    #Second Fundamental Form Coefs
    L = np.sum(xuu * n, axis=-1)
    M = np.sum(xuv * n, axis=-1)
    N = np.sum(xvv * n, axis=-1)

    return (L*N - M**2) / denom + 1


def comp_curve_data(c, scale, K_S3, du, dv):
    #Mean Curvature Calculations in R^3
    cu = np.gradient(c, du, axis=0, edge_order=2)
    cv = np.gradient(c, dv, axis=1, edge_order=2)

    cuu = np.gradient(cu, du, axis=0, edge_order=2)
    cuv = np.gradient(cu, dv, axis=1, edge_order=2)
    cvu = np.gradient(cv, du, axis=1, edge_order=2)
    cvv = np.gradient(cv, dv, axis=1, edge_order=2)

    
    #c_n Notation is Shortened to just x y and z of Parameterization
    xu, xv = cu[..., 0], cv[..., 0] 
    yu, yv = cu[..., 1], cv[..., 1]
    zu, zv = cu[..., 2], cv[..., 2]
    
    
    xuu, raw_xuv = cuu[..., 0], cuv[..., 0]
    raw_xvu, xvv = cvu[..., 0], cvv[..., 0]
    xuv = 0.5*(raw_xuv + raw_xvu)
    
    yuu, raw_yuv = cuu[..., 1], cuv[..., 1]
    raw_yvu, yvv = cvu[..., 1], cvv[..., 1]
    yuv = 0.5*(raw_yuv + raw_yvu)

    zuu, raw_zuv = cuu[..., 2], cuv[..., 2]
    raw_zvu, zvv = cvu[..., 2], cvv[..., 2]
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
    denom = (E*G - F**2)
    
    H = (L*G - 2*M*F + N*E) / (2 * denom)
    #Gaussian
    K_R3 = (L*N - M**2) / denom
    
    Scale_Factor = scale
    
    return CurvatureData(
    mean_curvature=H,
    gaussian_curvature=K_R3,
    scale_factor=Scale_Factor,
    )