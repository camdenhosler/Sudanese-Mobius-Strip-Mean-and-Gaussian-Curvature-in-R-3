#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def mobius_strip_s3(u, v, twist=0.5):

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
