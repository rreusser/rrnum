# rrnum

> Gobs of numerical computing code with maybe a couple salvageable pieces

## Summary

For my own sake, just posting this here so I don't lose it. This is from ca. 2005-2007 and reflects an attempt to re-use code. Because I thought numerical experiments in C were a good idea at the time. Yikes. There are probably like three or four useful things in here that can be salvaged.

## What's here

- `calc`
  - Adams Bashforth second order integration
  - Euler integration
  - Midpoint integration
  - RK4 integration
  - (Pretty naive) numerical differentiation
  - Finite difference stencil - I can't even tell exactly what this does...
- `data`
  - seems to be routines for calculating something like a linspace
- `fileio`
  - Tecplot output
  - Gnuplot output
- `general`
  - Array operations (add, subtract, fill)
  - HSV function
  - Functions (just some basics. Nothing interesting here)
- `grid`
  - Structured grid operations
    - Transfinite interpolation (filling in a 2D grid given border contours)
    - Elliptic smoothing (iteratively smooth a 2D structured grid) 
- `linalg`
  - Sparse matrices
    - row- and column-indexed sparse matrices. Basically queues up insertions, then sorts them and bakes the structure. Not ideal, but once baked, works reasonably well. This should be redone with the proper... interval tree? segment tree? Just RLE compression? So that the matrix doesn't have to be baked and reindexed.
    - Incomplete LU preconditioner (Requires converting row-indexed to column-indexed, which can probably be avoided at a log(N) cost by using a tree.)
    - Biconjugate Gradient solver (does this one work?)
    - Biconjugate Gradient Stabilized (BiCG-Stab) solver. Uses ILU preconditioning and sparse matrix multiplication. I recall this working pretty well, albeit slowly, on a 1000 &times; 1000 structured-curvilinear-grid PDE problem.)
    - Tridiagonal solver
    - Periodic tridiagonal solver
    - Pentadiagonal solver
  - Multigrid - Implements basic prolongation and restriction with Gauss-Seidel iteration on a rectilinear grid. Oof.
- `opengl` - Just direct draw mode for plotting. Nothing interesting here anymore.
  - Contour plots
  - Vector field plots
- `postscript`
  - An pluggable replacement for the OpenGL plotting that outputs to postscript. What was I thinking?

## License

&copy; 2016 Ricky Reusser. WTFPL.
