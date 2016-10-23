# sparse-linear-algebra

Numerical computation in native Haskell

TravisCI : [![Build Status](https://travis-ci.org/ocramz/sparse-linear-algebra.png)](https://travis-ci.org/ocramz/sparse-linear-algebra)

This library provides common numerical analysis functionality, without requiring any external bindings. It is not optimized for performance (yet), but it serves as an experimental platform for scientific computation in a purely functional setting.

Algorithms :

* Iterative linear solvers

    * BiConjugate Gradient (BCG)

    * Conjugate Gradient Squared (CGS)

    * BiConjugate Gradient Stabilized (BiCGSTAB) (non-Hermitian systems)

    * Transpose-Free Quasi-Minimal Residual (TFQMR)

* Matrix decompositions

    * QR factorization

    * LU factorization

* Eigenvalue algorithms

    * QR algorithm

    * Rayleigh quotient iteration

* Utilities : Vector and matrix norms, matrix condition number, Givens rotation, Householder reflection

* Predicates : Matrix orthogonality test (A^T A ~= I)


---------

## Usage

The module `Numeric.LinearAlgebra.Sparse` contains the interface functions:

To create a sparse matrix from an array of its entries

    fromListSM :: Foldable t => (Int, Int) -> t (IxRow, IxCol, a) -> SpMatrix a


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