
# Sudanese-Mobius-Strip-Mean-and-Gaussian-Curvature-in-R-3

## Overview
This project constructs the Sudanese Möbius strip in $S^3$ using the Lawson parameterization, stereographically projects it into $\mathbb{R}^3$, and visualizes the induced mean and Gaussian curvature. Although stereographic projection is conformal, it is not isometric; distances are therefore not preserved, and changes in both mean curvature and Gaussian curvature are expected. 

The mean curvature measures how the surface bends in space and is sensitive to both stereographic projection and the surface's orientation in $S^3$. In contrast, Gaussian curvature is intrinsic to the surface and depends only on distances along it; its distortion arises solely from the non-isometric nature of the projection and is unaffected by rotation in $S^3$. 

The mean curvature is still included to illustrate the surface's bending in space, though its precise values depend on the surface's orientation in $S^3$ and therefore are somewhat rotation-dependent. In contrast, Gaussian curvature provides a rotation-independent measure of intrinsic distortion.

This project focuses on numerical visualization and geometric interpretation rather than closed-form curvature analysis.

<p align="center">
  <img src="images/figure1.png" width="850" alt="Sudanese Mobius Strip Mean Curvature">
  <br>
  <i><b>Figure 1:</b> Stereographic projection of the Sudanese Möbius strip into $\mathbb{R}^3$. The surface is color-coded by Absolute Mean Curvature $|H|$.</i>
</p>

## Mathematical Background

**Sudanese Mobius Strip**

The Sudanese Mobius Strip is a minimal surface mobius strip in $S^3$ orginally parameterized by Lawson.  Its stereographic projection is unique in the fact that its boundary is coplanar and forms a circle.

**Stereographic Projection**

The Stereographic Projection from $S^3$ to $\mathbb{R}^3$ is calculated as below

$$ 
X_s = \frac{X_1}{1 - X_4} ,  
Y_s = \frac{X_2}{1 - X_4} ,
Z_s = \frac{X_3}{1 - X_4} $$

**Rotation in $S^3$**

The below matrix was used to rotate the surface in $S^3$ before projection

$$
\frac{1}{2}
\begin{bmatrix}
1 & -1 & -1 & -1 \\
1 & 1 & -1 & 1 \\
1 & 1 & 1 & -1 \\
1 & -1 & 1 & 1
\end{bmatrix}
$$

This matrix performs a rotation in $\mathbb{R}^4$ that reorients the Sudanese Möbius strip in $S^3$. 
It helps move the surface away from singularities (near the north pole) and prepares it for stereographic projection into $\mathbb{R}^3$. The factor of $\frac{1}{2}$ ensures that the transformation preserves scale appropriately.

**$\mathbb{R}^3$ Curvature**

For both of the following calculations the coefficients are equal since they are both in $\mathbb{R}^3$, each subscript represents a partial derrivative with respect to the subscripted variable, and n represents the unit normal vector of the surface.

The First Fundamental Form is given by

$$I =
\begin{bmatrix}
E & F \\
F & G
\end{bmatrix} =
\begin{bmatrix}
\mathbf{x}_u \cdot \mathbf{x}_u & \mathbf{x}_u \cdot \mathbf{x}_v \\
\mathbf{x}_u \cdot \mathbf{x}_v & \mathbf{x}_v \cdot \mathbf{x}_v
\end{bmatrix}
$$

The Second Fundamental Form is given by

$$
II =
\begin{bmatrix}
L & M \\
M & N
\end{bmatrix} =
\begin{bmatrix}
\mathbf{x}_{uu} \cdot \mathbf{n} & \mathbf{x}_{uv} \cdot \mathbf{n} \\
\mathbf{x}_{uv} \cdot \mathbf{n} & \mathbf{x}_{vv} \cdot \mathbf{n}
\end{bmatrix}
$$


The surface's mean curvature in $\mathbb{R}^3$ is written as $H$ and is calculated as shown below

$H = \frac{LG - 2MF + NE}{2(EG - F^2)}$.

The surface's Gaussian curvature in $\mathbb{R}^3$ is written as $K$ and is calculated as shown below

$K = \frac{LN - M**2}{(EG - F^2)}$.

$\mathbb{S}^3$ **curvature**

Since the strip is a minimal surface in $S^3$ it has zero mean curvature

The surface's Gaussian curvature is now calculated as shown below

$K = \frac{LN - M**2}{(EG - F^2)}$ + $K_0$

The first and second fundamental forms are the same as described above but are calculated using the original 4D parameterization in $S^3$ provided by Lawson. $K_0$ describes the curvature inherint to the space $S^3$ and is known to be 1.0.

**Notes**

Since Mobius Strips are non-orientable a plot of absolute mean curvature is shown, to remove the flipping due to the normal vector.

## Project Structure

**Dependencies**

-Numpy
-Matplotlib

**Usage**

-src/sudanese_mobius/ contains the parameterization, projection, and curvature calculation code

-scripts/ contains the script used to plot the surface

To generate the visualization run the plotting script directly after cloning.  The script imports code from src.

## Upcoming
Working on a notebooks folder with a jupyter lab simulation

## References

1. Lawson, H. B. (1970). *Complete minimal surfaces in* $S^3$. Annals of Mathematics, 92(2), 335–374.
