# %%
#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import numpy as np

from sudanese_mobius.parameters import GridParameters
from sudanese_mobius.S3_parameterization import mobius_strip_s3
from sudanese_mobius.stereo_projection import stereographic_projection
from sudanese_mobius.curvature import gaussian_curvature_S3
from sudanese_mobius.curvature import comp_curve_data
from sudanese_mobius.visualization import plot_surface

def main():
    params = GridParameters()
    u, v = np.meshgrid(params.u, params.v, indexing="ij")
    du = params.u[1] - params.u[0]
    dv = params.v[1] - params.v[0]
    
    x = mobius_strip_s3(u, v)
    c, scale = stereographic_projection(x)
    K_S3 = gaussian_curvature_S3(x, du, dv)
    curvatureData = comp_curve_data(c, scale, K_S3, du, dv)
    
    plot_surface(c, curvatureData, K_S3=K_S3, scale=scale, stride=4)

if __name__ == "__main__":
    main()