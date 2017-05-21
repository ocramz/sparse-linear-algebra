# sparse-linear-algebra

Numerical computation in native Haskell

[![Hackage](https://img.shields.io/hackage/v/sparse-linear-algebra.svg)](https://hackage.haskell.org/package/sparse-linear-algebra)  [![Build Status](https://travis-ci.org/ocramz/sparse-linear-algebra.png)](https://travis-ci.org/ocramz/sparse-linear-algebra)

This library provides common numerical analysis functionality, without requiring any external bindings. It is not optimized for performance (yet), but it serves as an experimental platform for scientific computation in a purely functional setting.

Contents :

* Iterative linear solvers (`<\>`)

    * Generalized Minimal Residual (GMRES) (non-Hermitian systems) 

    * BiConjugate Gradient (BCG)

    * Conjugate Gradient Squared (CGS)

    * BiConjugate Gradient Stabilized (BiCGSTAB) (non-Hermitian systems)

    * Moore-Penrose pseudoinverse (`pinv`) (rectangular systems)

* Direct linear solvers

    * LU-based (`luSolve`); forward and backward substitution (`triLowerSolve`, `triUpperSolve`)
    
* Matrix factorization algorithms

    * QR (`qr`)

    * LU (`lu`)

    * Cholesky (`chol`)

    * Arnoldi iteration (`arnoldi`)

* Eigenvalue algorithms

    * QR (`eigsQR`)

    * QR-Arnoldi (`eigsArnoldi`) 



* Utilities : Vector and matrix norms, matrix condition number, Givens rotation, Householder reflection

* Predicates : Matrix orthogonality test (A^T A ~= I)



### Under development

* Eigenvalue algorithms

    * Rayleigh quotient iteration (`eigRayleigh`)

* Matrix factorization algorithms

    * Golub-Kahan-Lanczos bidiagonalization (`gklBidiag`)
   
    * Singular value decomposition (SVD)

* Iterative linear solvers

    * Transpose-Free Quasi-Minimal Residual (TFQMR)

---------

## Examples

The module `Numeric.LinearAlgebra.Sparse` contains the user interface.

### Creation of sparse data

The `fromListSM` function creates a sparse matrix from a collection of its entries in (row, column, value) format. This is its type signature:

    fromListSM :: Foldable t => (Int, Int) -> t (IxRow, IxCol, a) -> SpMatrix a

and, in case you have a running GHCi session (the terminal is denoted from now on by `λ>`), you can try something like this:

    λ> amat = fromListSM (3,3) [(0,0,2),(1,0,4),(1,1,3),(1,2,2),(2,2,5)] :: SpMatrix Double

Similarly, `fromListSV` is used to create sparse vectors: 

    fromListSV :: Int -> [(Int, a)] -> SpVector a
    

Alternatively, the user can copy the contents of a list to a (dense) SpVector using

    fromListDenseSV :: Int -> [a] -> SpVector a



### Displaying sparse data

Both sparse vectors and matrices can be pretty-printed using `prd`:

    λ> prd amat

    ( 3 rows, 3 columns ) , 5 NZ ( density 55.556 % )

    2.00   , _      , _      
    4.00   , 3.00   , 2.00   
    _      , _      , 5.00       

*Note (sparse storage)*: sparse data should only contain non-zero entries not to waste memory and computation.

*Note (approximate output)*: `prd` rounds the results to two significant digits, and switches to scientific notation for large or small values. Moreover, values which are indistinguishable from 0 (see the `Numeric.Eps` module) are printed as `_`. 


### Matrix factorizations, matrix product

There are a few common matrix factorizations available; in the following example we compute the LU factorization of matrix `amat` and verify it with the matrix-matrix product `##` of its factors :

    λ> (l, u) <- lu amat
    λ> prd $ l ## u
    
    ( 3 rows, 3 columns ) , 9 NZ ( density 100.000 % )

    2.00   , _      , _      
    4.00   , 3.00   , 2.00   
    _      , _      , 5.00       


Notice that the result is _dense_, i.e. certain entries are numerically zero but have been inserted into the result along with all the others (thus taking up memory!).
To preserve sparsity, we can use a sparsifying matrix-matrix product `#~#`, which filters out all the elements x for which `|x| <= eps`, where `eps` (defined in `Numeric.Eps`) depends on the numerical type used (e.g. it is 10^-6 for `Float`s and 10^-12 for `Double`s).

    λ> prd $ l #~# u
    
    ( 3 rows, 3 columns ) , 5 NZ ( density 55.556 % )

    2.00   , _      , _      
    4.00   , 3.00   , 2.00   
    _      , _      , 5.00 


A matrix is transposed using the `transpose` function.

Sometimes we need to compute matrix-matrix transpose products, which is why the library offers the infix operators `#^#` (i.e. matrix transpose * matrix) and `##^` (matrix * matrix transpose):

    λ> amat' = amat #^# amat
    λ> prd amat'
    
    ( 3 rows, 3 columns ) , 9 NZ ( density 100.000 % )

    20.00  , 12.00  , 8.00   
    12.00  , 9.00   , 6.00   
    8.00   , 6.00   , 29.00      

    
    λ> lc <- chol amat'
    λ> prd $ lc ##^ lc
    
    ( 3 rows, 3 columns ) , 9 NZ ( density 100.000 % )

    20.00  , 12.00  , 8.00   
    12.00  , 9.00   , 10.80  
    8.00   , 10.80  , 29.00      


In the last example we have also shown the Cholesky decomposition (M = L L^T where L is a lower-triangular matrix), which is only defined for symmetric positive-definite matrices.

### Linear systems

Large sparse linear systems are best solved with iterative methods. `sparse-linear-algebra` provides a selection of these via the `<\>` (inspired by Matlab's "backslash" function. Currently this method uses GMRES as default) :

    λ> b = fromListDenseSV 3 [3,2,5] :: SpVector Double
    λ> x <- amat <\> b
    λ> prd x

    ( 3 elements ) ,  3 NZ ( density 100.000 % )

    1.50   , -2.00  , 1.00      


The result can be verified by computing the matrix-vector action `amat #> x`, which should (ideally) be very close to the right-hand side `b` :

    λ> prd $ amat #> x

    ( 3 elements ) ,  3 NZ ( density 100.000 % )

    3.00   , 2.00   , 5.00       
    

The library also provides a forward-backward substitution solver (`luSolve`) based on a triangular factorization of the system matrix (usually LU). This should be the preferred for solving smaller, dense systems. Using the LU factors defined previously we can cross-verify the two solution methods:

    λ> x' <- luSolve l u b
    λ> prd x'

    ( 3 elements ) ,  3 NZ ( density 100.000 % )

    1.50   , -2.00  , 1.00     









## License

GPL3, see LICENSE

## Credits

Inspired by

* `linear` : https://hackage.haskell.org/package/linear
* `vector-space` : https://hackage.haskell.org/package/vector-space
* `sparse-lin-alg` : https://github.com/laughedelic/sparse-lin-alg

## References

[1] Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., 2000

[2] G.H. Golub and C.F. Van Loan, Matrix Computations, 3rd ed., 1996

[3] T.A. Davis, Direct Methods for Sparse Linear Systems, 2006

[4] L.N. Trefethen, D. Bau, Numerical Linear Algebra, SIAM, 1997

[5] W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery, Numerical Recipes in Fortran 77, 2nd ed., 1992 