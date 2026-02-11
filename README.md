
# Sudanese Möbius Strip Visualization

## Overview
This project numerically generates and visualizes a Möbius strip embedded in 4D and projected into 3D, showing how curvature changes under stereographic projection.

 This project was built as a learning exercise in differential geometry, numerical curvature estimation, and interactive scientific visualization. 

The code computes local scale factors, the absolute mean curvature in $\mathbb{R}^3$, the Gaussian curvature in $\mathbb{R}^3$, and the difference in Gaussian curvature induced by the projection.  An interactive visualization is included which allows for switching between the scalar fields.

No new theoretical results are claimed.

## Key Features
* Parameterization of the Sudanese Möbius strip in $S^3$ (Lawson 1970)

* Stereographic projection of the Sudanese Möbius strip into $\mathbb{R}^3$

* Numerical Estimation of:
    *   Local scale factors in $\mathbb{R}^3$
    * Absolute mean curvature in $\mathbb{R}^3$
    * Gaussian curvature in $\mathbb{R}^3$
    * Difference in Gaussian curvature between $S^3$ and $\mathbb{R}^3$

* Interactive visualization of scalar fields with radio buttons

<p align="center">
  <img src="images/figure1.png" width="550" alt="Sudanese Möbius Strip Mean Curvature">
  <br>
  <i>Stereographic projection of the Sudanese Möbius strip into $\mathbb{R}^3$. The surface is color-coded by Absolute Mean Curvature $|H|$.</i>
</p>

## Quick Start

**Dependencies**

- Python 3.10+
- NumPy 1.26+
- Matplotlib 3.7+

**Install**

```bash
git clone https://github.com/camdenhosler/Sudanese-Mobius-Strip-Mean-and-Gaussian-Curvature-in-R-3.git
cd Sudanese-Mobius-Strip-Mean-and-Gaussian-Curvature-in-R-3

# Create a virtual environment
python -m venv venv

# On macOS/Linux:
source venv/bin/activate
# On Windows: 
# venv\Scripts\activate

pip install --upgrade pip
pip install -e .
```

**Run**

```
python scripts/generate_figure.py
```

Interactive plotting works automatically if a GUI backend is active. If running in a Jupyter notebook or IPython a Qt backend could be used (e.g. ```%matplotlib qt5```) 

## Mathematical Context

Since the Stereographic Projection is conformal but not isometric (angles are preserved but distances are not) curvature will be distorted.  The purpose of this project is to visualize this distortion.  The Sudanese Möbius strip was chosen for this project since it is both non-orientable and is a minimal surface in $S^3$.  The surface being minimal implies its mean curvature in $S^3$ is zero everywhere.  Due to the extrinsic nature of mean curvature (dependent on the surface's orientation in $S^3$ before projection) this property allows for easier visualization of distortion brought about by the projection.  Due to the non-orientability of the surface the absolute value of the mean curvature was taken.

## Limitations

* Curvatures are computed numerically using finite differences on a parameter grid

* Mean curvature values depend on the orientation in $S^3$ (since it's extrinsic)

* Although an orthogonal rotation matrix (unitary in $\mathbb{R}^4$) was used to rotate the surface of the singularity at the north pole, some artifacts may still remain

## References

1. Lawson, H. B. (1970). *Complete minimal surfaces in* $S^3$. Annals of Mathematics, 92(2), 335–374.

## License

[MIT License](LICENSE)
