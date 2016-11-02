# sparse-linear-algebra

Numerical computation in native Haskell

TravisCI : [![Build Status](https://travis-ci.org/ocramz/sparse-linear-algebra.png)](https://travis-ci.org/ocramz/sparse-linear-algebra)

This library provides common numerical analysis functionality, without requiring any external bindings. It is not optimized for performance (yet), but it serves as an experimental platform for scientific computation in a purely functional setting.

Contents :

* Iterative linear solvers

    * BiConjugate Gradient (`bcg`)

    * Conjugate Gradient Squared (`cgs`)

    * BiConjugate Gradient Stabilized (`bicgstab`) (non-Hermitian systems)

    * Transpose-Free Quasi-Minimal Residual (`tfqmr`)

* Direct linear solvers

    * LU-based (`luSolve`)

* Matrix factorization algorithms

    * QR (`qr`)

    * LU (`lu`)

    * Cholesky (`chol`)

* Eigenvalue algorithms

    * QR (`eigsQR`)

    * Rayleigh quotient iteration (`eigRayleigh`)

* Utilities : Vector and matrix norms, matrix condition number, Givens rotation, Householder reflection

* Predicates : Matrix orthogonality test (A^T A ~= I)


---------

## Examples

The module `Numeric.LinearAlgebra.Sparse` contains the user interface.

### Creation of sparse data

The `fromListSM` function creates a sparse matrix from an array of its entries we use :

    fromListSM :: Foldable t => (Int, Int) -> t (IxRow, IxCol, a) -> SpMatrix a

e.g.

    > amat = fromListSM (3,3) [(0,0,2),(1,0,4),(1,1,3),(1,2,2),(2,2,5)]

and similarly

    fromListSV :: Int -> [(Int, a)] -> SpVector a

can be used to create sparse vectors.

### Displaying sparse data

Both sparse vectors and matrices can be pretty-printed using `prd`:

    > prd amat
    ( 3 rows, 3 columns ) , 5 NZ ( sparsity 0.5555555555555556 )

    [2,0,0]
    [4,3,2]
    [0,0,5]

The zeros are just added at printing time; sparse vectors and matrices should only contain non-zero entries.

### Matrix operations

Matrix factorizations are available as `lu` and `qr` respectively, and are straightforward to verify by using the matrix product `##`  :

    > (l, u) = lu amat
    > prd $ l ## u
    ( 3 rows, 3 columns ) , 9 NZ ( sparsity 1.0 )

    [2.0,0.0,0.0]
    [4.0,3.0,2.0]
    [0.0,0.0,5.0]

Notice that the result is _dense_, i.e. certain entries are numerically zero but have been inserted into the result along with all the others (thus taking up memory!).
To preserve sparsity, we can use a sparsifying matrix-matrix product `#~#`, which filters out all the elements x for which `|x| <= eps`, where `eps` (defined in `Numeric.Eps`) depends on the numerical type used.

    > prd $ l #~# u
    ( 3 rows, 3 columns ) , 5 NZ ( sparsity 0.5555555555555556 )

    [2.0,0.0,0.0]
    [4.0,3.0,2.0]
    [0.0,0.0,5.0]


### Linear systems

Large sparse linear systems are best solved with iterative methods. `sparse-linear-algebra` provides a selection of these via the `linSolve` function, or alternatively `<\>` (which uses BiCGSTAB as default) :

    > b = fromListSV 3 [(0,3),(1,2),(2,5)]
    > x = amat <\> b
    > prd x
    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    [1.4999999999999998,-1.9999999999999998,0.9999999999999998]

The result can be verified by computing the matrix-vector action `amat #> x`, which should (ideally) be very close to the right-hand side `b` :

    > prd $ amat #> x
    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    [2.9999999999999996,1.9999999999999996,4.999999999999999]

The library also provides a forward-backward substitution solver (`luSolve`) based on a triangular factorization of the system matrix (usually LU). This should be the preferred for solving smaller, dense systems. Using the data defined above:

    > xlu = luSolve l u b
    > prd xlu

    ( 3 elements ) ,  3 NZ ( sparsity 1.0 )

    [1.5,-2.0,1.0]





----------

This is also an experiment in principled scientific programming :

* set the stage by declaring typeclasses and some useful generic operations (normed linear vector spaces, i.e. finite-dimensional spaces equipped with an inner product that induces a distance function),

* define appropriate data structures, and how they relate to those properties (sparse vectors and matrices, defined internally via `Data.IntMap`, are made instances of the VectorSpace and Additive classes respectively). This allows to decouple the algorithms from the actual implementation of the backend,

* implement the algorithms, following 1:1 the textbook [1] 


## License

GPL3, see LICENSE

## Credits

Inspired by

* `linear` : https://hackage.haskell.org/package/linear
* `sparse-lin-alg` : https://github.com/laughedelic/sparse-lin-alg

## References

[1] : Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., 2000