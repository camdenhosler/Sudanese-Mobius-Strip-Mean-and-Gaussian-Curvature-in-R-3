#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def stereographic_projection(x):
    """
    Stereographically Projects the orginial parameterization in S3 to R3 and
    calculates the scale factor at each point due to the projection

    Parameters
    ----------
    x : ndarray
        Array of points (shape: (..., 4)

    Returns
    -------
    c : ndarray
        Array of points (shape: (..., 3) representing the projection of the 
        surface in R3.
    scale: ndarray
        Array of scalar values at each point (shape: (...,) representing 
        the surface's scale factor in R3.
    """ 
    denom = 1.0 - x[..., 3]
    c = x[..., :3] / denom[..., None]

    # Rotate coordinates for visualization
    rot = np.array([
        [1, 0, 0],
        [0, np.sqrt(2)/2, -np.sqrt(2)/2],
        [0, np.sqrt(2)/2,  np.sqrt(2)/2]
    ])

    c = c @ rot.T
    scale = 1.0 / denom

    return c, scale