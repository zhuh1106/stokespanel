# stokespanel
stokescloseeval.m - Stokes layer potential close evaluation via panel quadrature. (both interior and exterior case at this point)

test_stokesSD.m - routine for self convergence test. Used a default density tau=(sin t,cos t). In practice, tau need to be first solved on the panel nodes; can use "quadr.m" to set up the (Legendre) panel quadrature.

stokespanelcor.m - Panel correction scheme for Stokes, gives special quadrature value.

quadr.m - set up Legendre panel quadrature

gauss.m - find Gauss nodes and weights and diff matrix on [-1, 1]
