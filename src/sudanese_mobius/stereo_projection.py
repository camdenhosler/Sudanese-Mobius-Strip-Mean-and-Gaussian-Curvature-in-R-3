#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def stereographic_projection(x):
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