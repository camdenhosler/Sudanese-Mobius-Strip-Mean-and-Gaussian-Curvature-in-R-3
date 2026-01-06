#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def mobius_strip_s3(u, v, twist=0.5):
    """
    Parameterization of the Mobius strip in S3 as devised by Lawson.

    The parameter `twist` controls how the surface is embedded in S3.
    Setting `twist=0.5` produces a Mobius strip; setting `twist=1.0`
    yields a torus instead. The resulting surface is rotated slightly
    to avoid contact with the north pole of the 3-sphere.

    Parameters
    ----------
    u : ndarray
        First grid parameter (e.g., `u_min ≤ u ≤ u_max`).
    v : ndarray
        Second grid parameter (e.g., `v_min ≤ v ≤ v_max`).
    twist : float
        Embedding twist factor set to 0.5.

    Returns
    -------
    x : ndarray
        Array of points (shape: (..., 4) representing the surface in S3.
    """ 
    x1 = np.cos(u) * np.cos(v)
    x2 = np.cos(u) * np.sin(v)
    x3 = np.sin(u) * np.cos(twist * v)
    x4 = np.sin(u) * np.sin(twist * v)

    x = np.stack([x1, x2, x3, x4], axis=-1)

    # Rotation to avoid pole singularity
    R = 0.5 * np.array([
        [ 1, -1, -1, -1],
        [ 1,  1, -1,  1],
        [ 1,  1,  1, -1],
        [ 1, -1,  1,  1]
    ])

    return x @ R.T
