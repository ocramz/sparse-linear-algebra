# Cholesky Factorization Test Suite

This document describes the comprehensive property-based test suite for the Cholesky factorization implementation in `sparse-linear-algebra`.

## Overview

The Cholesky factorization decomposes a positive definite matrix **A** into a product **L L^H**, where **L** is a lower triangular matrix and **L^H** is the conjugate transpose of **L**. For real matrices, this becomes **L L^T**.

## Implementation Fix

The original Cholesky implementation had two bugs for complex-valued matrices:

### Fix 1: Diagonal Element Computation
Changed from using `**2` to proper complex conjugate multiplication:

```haskell
-- Before (incorrect for complex numbers):
let l = sum (fmap (**2) lrow)

-- After (correct for both real and complex):
let l = sum (fmap (\x -> x * conj x) lrow)
```

This ensures that for complex numbers, we compute `x * conj(x)` which gives the squared magnitude, rather than `x^2` which would be incorrect.

### Fix 2: Subdiagonal Element Computation  
Changed from using `contractSub` (which doesn't use conjugation) to using inner product:

```haskell
-- Before (incorrect for complex numbers):
inn = contractSub ll ll i j (j - 1)

-- After (correct for both real and complex):
lrow_i = ifilterSV (\k _ -> k <= j - 1) (extractRow ll i)
lrow_j = ifilterSV (\k _ -> k <= j - 1) (extractRow ll j)
inn = lrow_i <.> lrow_j
```

The inner product operator `<.>` correctly uses conjugation for complex numbers: `v <.> w = sum (v_i * conj(w_i))`.

These fixes enable the Cholesky factorization to work correctly for both real symmetric positive definite matrices and complex Hermitian positive definite matrices.

## Test Structure

The test suite consists of two main categories:

### 1. Real Number Tests (Symmetric Positive Definite Matrices)

#### Specific Test Cases:
- **3x3 SPD matrix**: A small dense SPD matrix for basic verification
- **5x5 tridiagonal SPD matrix**: A sparse tridiagonal matrix to test sparsity handling

#### Property-based Tests:
1. **Reconstruction Property**: `L L^T = A`
   - Verifies that the factorization correctly reconstructs the original matrix
   - Uses Frobenius norm to measure reconstruction error
   
2. **Lower Triangular Property**: `L is lower triangular`
   - Ensures the factor L has no non-zero elements above the diagonal
   
3. **Positive Diagonal Property**: All diagonal elements of `L` are positive
   - Validates that diagonal entries `L[i,i] > 0` for all i

### 2. Complex Number Tests (Hermitian Positive Definite Matrices)

#### Specific Test Cases:
- **2x2 HPD matrix**: A small Hermitian matrix for basic verification
- **3x3 HPD matrix**: A slightly larger Hermitian matrix

#### Property-based Tests:
1. **Reconstruction Property**: `L L^H = A`
   - Verifies that the factorization correctly reconstructs the original Hermitian matrix
   - Uses Frobenius norm to measure reconstruction error
   
2. **Lower Triangular Property**: `L is lower triangular`
   - Ensures the factor L has no non-zero elements above the diagonal
   
3. **Positive Real Diagonal Property**: All diagonal elements of `L` have positive real parts
   - Validates that `realPart(L[i,i]) > 0` for all i
   - For Hermitian matrices, the diagonal of L is typically real and positive

## QuickCheck Generators

### PropMatSPD (Symmetric Positive Definite)
Generates random SPD matrices using the construction:
```
SPD = M^T M + εI
```
where M is a random matrix and ε = 0.1 is a small constant ensuring positive definiteness.

### PropMatHPD (Hermitian Positive Definite)
Generates random HPD matrices using the construction:
```
HPD = M^H M + εI
```
where M is a random complex matrix and ε = 0.1 + 0i.

## Matrix Sizes

The property-based tests use QuickCheck's `sized` combinator to generate matrices of varying dimensions (typically 2x2 to ~10x10), ensuring the factorization works correctly across different scales.

## Running the Tests

```bash
stack test
```

Or with make:
```bash
make test
```

## Test Coverage

The test suite provides extensive coverage of:
- ✅ Real-valued symmetric positive definite matrices (sparse and dense)
- ✅ Complex-valued Hermitian positive definite matrices (dense)
- ✅ Various matrix sizes (2x2 through ~10x10)
- ✅ Mathematical properties (reconstruction, triangularity, positive diagonal)
- ✅ Edge cases via QuickCheck property-based testing

## Notes

1. The Cholesky factorization requires the input matrix to be positive definite. The algorithm will throw a `NeedsPivoting` exception if a zero diagonal element is encountered.

2. For numerical stability, the test generators add a small diagonal perturbation (ε = 0.1) to ensure matrices are strictly positive definite.

3. All tests use the `nearZero` function from the `Numeric.Eps` module to handle floating-point comparison tolerances appropriately.

## References

- Golub, G. H., & Van Loan, C. F. (2013). Matrix Computations (4th ed.). Johns Hopkins University Press.
- The Cholesky-Banachiewicz algorithm implemented proceeds row-by-row from the upper-left corner.
