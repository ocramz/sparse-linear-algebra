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

    * LU-based (`luSolve`)
    
* Matrix factorization algorithms

    * QR (`qr`)

    * LU (`lu`)

    * Cholesky (`chol`)

    * Arnoldi iteration (`arnoldi`)

* Eigenvalue algorithms

    * QR (`eigsQR`)

    * Rayleigh quotient iteration (`eigRayleigh`)

* Utilities : Vector and matrix norms, matrix condition number, Givens rotation, Householder reflection

* Predicates : Matrix orthogonality test (A^T A ~= I)



### Under development

* Matrix factorization algorithms

    * Golub-Kahan-Lanczos bidiagonalization
   
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
    ( 3 rows, 3 columns ) , 5 NZ ( sparsity 0.5555555555555556 )

    2.0  _   _ 
    4.0 3.0 2.0
     _   _  5.0

Note: sparse data should only contain non-zero entries not to waste memory and computation.

### Matrix operations

There are a few common matrix factorizations available; in the following example we compute the LU factorization of matrix `amat` and verify it with the matrix-matrix product `##` of its factors :

    λ> (l, u) <- lu amat
    λ> prd $ l ## u
    ( 3 rows, 3 columns ) , 9 NZ ( sparsity 1.0 )

    2.0  _   _ 
    4.0 3.0 2.0
     _   _  5.0

Notice that the result is _dense_, i.e. certain entries are numerically zero but have been inserted into the result along with all the others (thus taking up memory!).
To preserve sparsity, we can use a sparsifying matrix-matrix product `#~#`, which filters out all the elements x for which `|x| <= eps`, where `eps` (defined in `Numeric.Eps`) depends on the numerical type used (e.g. it is 10^-6 for `Float`s and 10^-12 for `Double`s).

    λ> prd $ l #~# u
    ( 3 rows, 3 columns ) , 5 NZ ( sparsity 0.5555555555555556 )

    2.0  _   _ 
    4.0 3.0 2.0
     _   _  5.0    

A matrix is transposed using the `transpose` function.

Sometimes we need to compute matrix-matrix transpose products, which is why the library offers the infix operators `#^#` (i.e. matrix transpose * matrix) and `##^` (matrix * matrix transpose):

    λ> amat' = amat #^# amat
    λ> prd amat'
    ( 3 rows, 3 columns ) , 9 NZ ( sparsity 1.0 )

    20.0 12.0 8.0
    12.0 9.0 6.0
    8.0 6.0 29.0
    
    λ> l <- chol amat'
    λ> prd $ l ##^ l
    ( 3 rows, 3 columns ) , 9 NZ ( sparsity 1.0 )

    20.000000000000004 12.0 8.0
    12.0 9.0 10.8
    8.0 10.8 29.0

In the last example we have also shown the Cholesky decomposition (M = L L^T where L is a lower-triangular matrix), which is only defined for symmetric positive-definite matrices.

### Linear systems

Large sparse linear systems are best solved with iterative methods. `sparse-linear-algebra` provides a selection of these via the `<\>` (inspired by Matlab's "backslash" function. Here we use GMRES as default solver method) :

    λ> b = fromListDenseSV 3 [3,2,5] :: SpVector Double
    λ> x <- amat <\> b
    λ> prd x
    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    1.4999999999999998 -1.9999999999999998 0.9999999999999998

The result can be verified by computing the matrix-vector action `amat #> x`, which should (ideally) be very close to the right-hand side `b` :

    λ> prd $ amat #> x
    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    2.9999999999999996 1.9999999999999996 4.999999999999999

The library also provides a forward-backward substitution solver (`luSolve`) based on a triangular factorization of the system matrix (usually LU). This should be the preferred for solving smaller, dense systems. Using the data defined above we can cross-verify the two solution methods:

    λ> x' <- luSolve l u b
    λ> prd x'

    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    1.5 -2.0 1.0








## License

GPL3, see LICENSE

## Credits

Inspired by

* `linear` : https://hackage.haskell.org/package/linear
* `vector-space` : https://hackage.haskell.org/package/vector-space
* `sparse-lin-alg` : https://github.com/laughedelic/sparse-lin-alg

## References

[1] : Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., 2000

[2] : L. N. Trefethen, D. Bau, Numerical Linear Algebra, SIAM, 1997