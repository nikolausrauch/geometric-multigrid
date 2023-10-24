# Geometric Multigrid Summary

The idea of multigrid methods is motivated by the "nature" of the error of the approximation an iterative solver produces.
While an iterative solver tends to "smooth out" the high frequency components of the error early on in the iteration process, the low frequencies still dominate.
However, a low frequency error can be represented in a coarser grid where it is computationally less expensive to solve for and then "transfered" back to the finer grid to correct the approximation.  
The aforementioned idea is not just limited to two grids but can be further extended to a number of intermediate grids, up to the point where we have to solve for a single grid point.
In the following I'll outline the inner workings of a **geometric multigrid** solver applied to the *Poisson's equation* in 1D and 2D.

## 1D Poisson's equation 
Given is the one dimensional Poisson's equation with zero Dirichlet boundary conditions
$$
    \partial_{x}^2u=f(x), \quad \Omega=(0,1)\\
    u(0) = u(1) = 0
$$
We can compute numerically the steady-state solution, using the central difference to discretize the equation
$$
    \frac{-u_{i-1} + 2u_i - u_{i+1}}{h^2} = f_i
$$
where $h$ represents the grid spacing.   
In case $u(0)$ and $u(1)$ are given directly by the Dirichlet boundary conditions $u(0) = u_0,~ u(1) = u_{N+1}$, we can move them to the right hand side for the respective discretization points, which leaves us with the following system of linear equations for the unknowns $u_i, i = 1...N$:
$$
    \begin{bmatrix}
        2 & -1                          \\
        -1 & 2 & -1                     \\
           & -1 & \ddots & \ddots       \\
           &    & \ddots & \ddots & -1  \\
           &    &        & -1     & 2   \\
    \end{bmatrix}
    \begin{bmatrix}
        u_1     \\
        u_2     \\
        \vdots  \\
        \vdots  \\
        u_{N}
    \end{bmatrix}
    =
    \begin{bmatrix}
        h^2 f_1 + u_0       \\
        h^2 f_2             \\
        \vdots              \\
        \vdots              \\
        h^2 f_N + u_{N+1}
    \end{bmatrix}
$$

## 2D Poisson's equation
In two dimension the Poisson's equation is given as
$$
    -\partial_{x}^2 u - \partial_{y}^2 u = f(x, y) \qquad \Omega = (0, 1)^2\\
    u = 0 \text{ on } \partial\Omega
$$
and analogously to the one dimensional case, the central differences yields the following system of linear
equations
$$
    \frac{-u_{i-1,j} + 2u_{ij} - u_{i+1, j}}{h_x^2} + \frac{-u_{i,j-1} + 2u_{ij} - u_{i, j+1}}{h_y^2} = f_{ij}.
$$
Written in matrix form (for $h_x = h_y$) we get
$$
    \begin{bmatrix}
        \mathbf{D} & -\mathbf{I}                        \\
        -\mathbf{I} & \mathbf{D} & -\mathbf{I}          \\
           & -\mathbf{I} & \ddots & \ddots              \\
           &    & \ddots & \ddots & -\mathbf{I}         \\
           &    &        & -\mathbf{I}     & \mathbf{D} \\
    \end{bmatrix}
    \begin{bmatrix}
        u_{1,1}     \\
        u_{1,2}     \\
        \vdots  \\
        \vdots  \\
        u_{N,N}
    \end{bmatrix}
    =
    \begin{bmatrix}
        h^2 f_{1,1} + u_{0,1} + u_{1, 0}    \\
        h^2 f_{1,2} + u_{2,0}               \\
        \vdots                              \\
        \vdots                              \\
        h^2 f_{N,N} + u_{N,N+1} + u_{N+1,N}
    \end{bmatrix}
$$
where $\mathbf{I}$ and $\mathbf{D}$ are $N \times N$ matrices given by
$$
    \mathbf{I}
    =
    \begin{bmatrix}
        1  &                         \\
           & 1 &                     \\
           &   & \ddots &            \\
           &   &        & \ddots & 0 \\
           &   &        &  0     & 1 \\
    \end{bmatrix},
    \qquad
    \mathbf{D}
    =
    \begin{bmatrix}
        4 & -1                          \\
        -1 & 4 & -1                     \\
           & -1 & \ddots & \ddots       \\
           &    & \ddots & \ddots & -1  \\
           &    &        & -1     & 4   \\
    \end{bmatrix}
$$

## Geometric Multigrid Method

Iterative solvers (such as *Jacobi, Gauß-Seidel, SOR*), are capable of eliminating the high frequency errors in the solution rather quickly, i.e. after a few iterations.
However, they perform poorly in reducing low frequency errors - especially on high resolution grids, and as such their convergence rate stalls once the error becomes "smooth".
The fundamental concept behind multigrid methods is to maintain the fast convergence rate, by employing a hierarchy of physical grids at decreasing resolutions at which the **low frequency error appears to be of high frequency**.   
Consider a system of linear equations written in matrix form
$$
    \mathbf{Au} = \mathbf{f}
$$
Let $\mathbf{\tilde{u}}$ be the approximation to the solution $\mathbf{u}$, then the error $\mathbf{e}$ and the residual $\mathbf{r}$ are given by
$$
    \mathbf{e} = \mathbf{u} - \mathbf{\tilde{u}}\\
    \mathbf{r} = \mathbf{f} - \mathbf{A\tilde{u}}
$$
and by plugging them into the equation above we get the following relation between the error $\mathbf{e}$ and the residual $\mathbf{r}$:
$$
    \mathbf{A(u - e)} = \mathbf{r} - \mathbf{A\tilde{u}}\\
    \mathbf{Ae} = \mathbf{r}
$$
Considering this relation, the core idea is to first approximate the solution $\mathbf{\tilde{u}}$ on a fine grid with a conventional iterative solver in few iterations (**pre-smoothing**), and then transfer the residual to a coarse grid, by a process typically referred to as **restriction**.   
After solving $\mathbf{Ae} = \mathbf{r}$ with the down scaled residual on the coarse grid, the solution i.e. the error, is transferred back to the fine grid (**prolongation or interpolation**), and the approximation $\mathbf{\tilde{u}}$ is corrected by
$$
    \mathbf{u} = \mathbf{\tilde{u} + e}
$$
Often, the corrected solution is again smoothed with a couple of iterations (**post-smoothing**).
Applying these steps recursively on a whole hierarchy of grids forms the basis of a multigrid-cycle.
Depending on the order the solver visits the different grids i.e. the call-graph, we refer to the cycles as either a V-cyle ($\mu = 1$) or W-cycle ($\mu = 2$).

![algo](https://github.com/nikolausrauch/geometric_multigrid/assets/13553309/9cc0be10-f7f1-4c6b-a6fe-5a5d499dfad5)

![multigrid_cycles](https://github.com/nikolausrauch/geometric_multigrid/assets/13553309/71623710-865a-4ef7-906f-80800f49a971)

## Prolongation and Restriction
At its core the multigrid method relies on the two operations **prolongation and restriction**, to transfer the residual from fine to coarse, and the error from coarse to fine grids, respectively.
In the following I refer to a fine grid with grid spacing $h$ as $\Omega^h$, and a coarse grid with spacing $2h$ as $\Omega^{2h}$.
In **one** dimension the number of grid points in $\Omega^h$ is $2^{n+1} - 1$, and $2^n - 1$ for the corresponding $\Omega^{2h}$.
In two dimension we consider square grids, with $(2^{n+1} - 1)^2$ points in $\Omega^h$ and $(2^n - 1)^2$ points in $\Omega^{2h}$.
Therefore, the grid size of the original problem (the highest level) needs to be chosen with care.
For the remainder of this section, the vector containing values on the fine grid is denoted as $\mathbf{u}$ and the coarse grid with $\mathbf{v}$.
Then the interpolation operator $\mathbf{I}^h_{2h}$ and the restriction operator $\mathbf{R}^{2h}_{h}$ are defined as matrices, satisfying the following relations
$$
    \mathbf{I}^h_{2h} \mathbf{v} = \mathbf{u},\\
    \mathbf{R}^{2h}_{h} \mathbf{u} = \mathbf{v}\\
$$
In addition to the the projection of the residual and the error between grid levels, we need the discretization matrices $\mathbf{A}_{2ih}$ for each coarse grid $\Omega^{2ih$}.
They can either be constructed from discretising the original problem, or be computed with the interpolation and restriction matrices from each grid.
The matrix $\mathbf{A}_{2h}$ for $\Omega^{2h}$ is given by
$$
    \mathbf{A}_{2h} = \mathbf{R}^{2h}_{h} \mathbf{A}_{h} \mathbf{I}^h_{2h}
$$

### Interpolation and restriction in 1D

A simple choice for the prolongation step is to use *linear interpolation*. 
Fine grid points that are directly aligned with points from the coarse grid, adopt the corresponding values, while the remaining points take the average of their neighbours:
$$
    u_{2i} = v_i,\\
    u_{2i+1} = \frac{1}{2}(v_i + v_{i+1}), \qquad 0 \leq i \leq 2^n - 1
$$
Note, that the values on the borders $v0$, $v_{2^n}$ are assumed to be zero.
The **restriction and interpolation** only concern the interior values, but are included in the image below due to their role in the interpolation.
The **interpolation** operator $\mathbf{I}^h_{2h}$ can now be written as the $(2^{n+1} - 1) \times (2^n - 1)$ matrix:
$$
    \mathbf{I}^h_{2h}
    =
    \begin{bmatrix}
        1  &                           \\
        2  & 1 &                       \\
        1  & 2 &                       \\
           & 1 & \ddots &              \\
           &   &        & \ddots  & 1  \\
           &   &        &         & 2  \\
           &   &        &         & 1  \\
    \end{bmatrix}.
$$
The **restriction** operator can be implemented by taking weighted averages of neighbouring fine grid points (in this context referred to as *full weighting*):
$$
    u_i = \frac{1}{4}(v_{2i-1} + 2v_{2i} + v_{2i+1}).
$$
Incidentally, the restriction matrix $\mathbf{R}^{2h}_{h}$ can be retrieved by *transposition* of the previous interpolation operator:
$$
    \mathbf{R}^{2h}_{h}
    =
    \frac{1}{2}(\mathbf{I}^h_{2h})^T
    =
    \begin{bmatrix}
        1  & 2 & 1                            \\
           &   & 1 & 2 & 1                    \\
           \\
           &   &   &   & \ddots               \\
           &   &   &   &         & \ddots     \\
           &   &   &    &        &  1 & 2 & 1 \\
    \end{bmatrix}.
$$

![1D_transfer](https://github.com/nikolausrauch/geometric_multigrid/assets/13553309/f5a8a8c4-6f52-47cc-9576-bd9be550e63b)

### Interpolation and restriction in 2D
The prolongation in two dimension is realized by a *bilinear interpolation*:
$$
    u_{2i, 2j} = v_{i,j} \\
    u_{2i+1, 2j} = \frac{1}{2}(v_{i, j} + v_{i+1, j}), \qquad u_{2i, 2j+1} = \frac{1}{2}(v_{i, j} + v_{i, j+}) \\
    u_{2i+1, 2j+1} = \frac{1}{4}(v_{i, j} + v_{i+1, j} + v_{i, j+1} + v_{i+1, j+1})
$$
The resulting interpolation matrix $\mathbf{I}^h_{2h}$ has dimensions $(2^{n+1} - 1)^2 \times (2^n - 1)^2$.
In a similar manner to the one dimensional case, we can compute the restriction $\mathbf{R}^{2h}_{h}$ by multiplying the *transpose* of $\mathbf{I}^h_{2h}$ with $\frac{1}{4}$, yielding:
$$
    v_{ij} = \frac{1}{16} ( 4u_{2i, 2j} + 2(u_{2i-1, 2j} + u_{2i+1, 2j} + u_{2i, 2j-1} + u_{2i, 2j+1}) \\ 
                                        +  (u_{2i-1, 2j-1} + u_{2i+1, 2j+1} + u_{2i+1, 2j-1} + u_{2i-1, 2j+1}) )
$$
Therefore, each value in $\Omega^{2h}$ is a weighted average of nine grid values of $\Omega^h$:

![2D_transfer](https://github.com/nikolausrauch/geometric_multigrid/assets/13553309/4c2b08da-ed1d-47b2-85e9-49a1086c0272)
