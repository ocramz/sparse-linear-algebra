# Summary of Functions Commented Out for -Wall -Werror Compliance

This document lists all the functions that were commented out during the process of enabling strict `-Wall -Werror` compilation flags. These functions were unused in the codebase but may be useful for future development.

## File: src/Numeric/LinearAlgebra/Sparse.hs

### 1. BCG (Biconjugate Gradient) Solver Functions

#### `bcgInit`
**Lines**: 882-889 (commented)
```haskell
bcgInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> BCG a
bcgInit aa b x0 = BCG x0 r0 r0hat p0 p0hat where
  r0 = b ^-^ (aa #> x0)    
  r0hat = r0
  p0 = r0
  p0hat = r0
```
**Purpose**: Initializes the BCG solver state with initial guess and computes initial residuals

#### `bcgStep`
**Lines**: 891-908 (commented)
```haskell
bcgStep :: (MatrixType (SpVector a) ~ SpMatrix a,
      LinearVectorSpace (SpVector a), InnerSpace (SpVector a),
      MatrixRing (SpMatrix a), Fractional (Scalar (SpVector a))) =>
     SpMatrix a -> BCG a -> BCG a 
bcgStep aa (BCG x r rhat p phat) = BCG x1 r1 rhat1 p1 phat1 where
    aap = aa #> p
    alpha = (r `dot` rhat) / (aap `dot` phat)
    x1 = x ^+^ (alpha .* p)
    r1 = r ^-^ (alpha .* aap)
    rhat1 = rhat ^-^ (alpha .* (transpose aa #> phat))
    beta = (r1 `dot` rhat1) / (r `dot` rhat)
    p1 = r1 ^+^ (beta .* p)
    phat1 = rhat1 ^+^ (beta .* phat)
```
**Purpose**: Performs one iteration step of the BCG algorithm

---

### 2. CGS (Conjugate Gradient Squared) Solver Functions

#### `cgsInit`
**Lines**: 917-921 (commented)
```haskell
cgsInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> CGS a
cgsInit aa b x0 = CGS x0 r0 r0 r0 where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
```
**Purpose**: Initializes the CGS solver state

#### `cgsStep`
**Lines**: 923-945 (commented)
```haskell
cgsStep :: (V (SpVector a), Fractional (Scalar (SpVector a))) =>
     MatrixType (SpVector a) -> SpVector a -> CGS a -> CGS a
cgsStep aa rhat (CGS x r p u) = CGS xj1 rj1 pj1 uj1
    where
    aap = aa #> p
    alphaj = (r `dot` rhat) / (aap `dot` rhat)
    q = u ^-^ (alphaj .* aap)
    xj1 = x ^+^ (alphaj .* (u ^+^ q))         -- updated solution
    rj1 = r ^-^ (alphaj .* (aa #> (u ^+^ q))) -- updated residual
    betaj = (rj1 `dot` rhat) / (r `dot` rhat)
    uj1 = rj1 ^+^ (betaj .* q)
    pj1 = uj1 ^+^ (betaj .* (q ^+^ (betaj .* p)))
```
**Purpose**: Performs one iteration step of the CGS algorithm

---

### 3. BiCGSTAB (Biconjugate Gradient Stabilized) Solver Functions

#### `bicgsInit`
**Lines**: 950-954 (commented)
```haskell
bicgsInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> BICGSTAB a
bicgsInit aa b x0 = BICGSTAB x0 r0 r0 where
  r0 = b ^-^ (aa #> x0)   -- residual of initial guess solution
```
**Purpose**: Initializes the BiCGSTAB solver state

#### `bicgstabStep`
**Lines**: 956-971 (commented)
```haskell
bicgstabStep :: (V (SpVector a), Fractional (Scalar (SpVector a))) =>
     MatrixType (SpVector a) -> SpVector a -> BICGSTAB a -> BICGSTAB a
bicgstabStep aa r0hat (BICGSTAB x r p) = BICGSTAB xj1 rj1 pj1 where
     aap = aa #> p
     alphaj = (r <.> r0hat) / (aap <.> r0hat)
     sj = r ^-^ (alphaj .* aap)
     aasj = aa #> sj
     omegaj = (aasj <.> sj) / (aasj <.> aasj)
     xj1 = x ^+^ (alphaj .* p) ^+^ (omegaj .* sj)    -- updated solution
     rj1 = sj ^-^ (omegaj .* aasj)
     betaj = (rj1 <.> r0hat)/(r <.> r0hat) * alphaj / omegaj
     pj1 = rj1 ^+^ (betaj .* (p ^-^ (omegaj .* aap)))
```
**Purpose**: Performs one iteration step of the BiCGSTAB algorithm

---

### 4. Moore-Penrose Pseudoinverse

#### `pinv`
**Lines**: 980-985 (commented)
```haskell
pinv :: (LinearSystem v, MatrixRing (MatrixType v), MonadThrow m, MonadIO m) =>
     MatrixType v -> v -> m v
pinv aa b = (aa #^# aa) <\> atb where
  atb = transpose aa #> b
```
**Purpose**: Computes least-squares approximation of a rectangular system using Moore-Penrose pseudoinverse
**Note**: This function was also removed from the module export list (line 20)

---

### 5. IterativeSolver Typeclass

#### `IterativeSolver` class
**Lines**: 1027-1030 (commented)
```haskell
class IterativeSolver s where
  -- solver :: 
```
**Purpose**: Placeholder typeclass for iterative solver abstraction (was incomplete)

---

## Summary Statistics

- **Total functions commented out**: 10
  - BCG-related: 2 functions
  - CGS-related: 2 functions  
  - BiCGSTAB-related: 2 functions
  - Pseudoinverse: 1 function
  - Typeclass: 1 class definition
  - Show instances: 3 instances (for BCG, CGS, BICGSTAB)

## Reason for Commenting Out

All these functions were commented out because they triggered the `-Wunused-top-binds` warning under `-Wall -Werror` compilation. While these are sophisticated iterative solver implementations that could be valuable for solving sparse linear systems, they were not currently used anywhere in the codebase.

## Recommendations for Re-enabling

To re-enable these functions in the future:

1. **Export them from the module** - Add them to the export list in the module header
2. **Use them in examples or tests** - Create usage examples demonstrating their functionality
3. **Document them properly** - Add comprehensive Haddock documentation
4. **Add to the public API** - If they're meant to be part of the public API, ensure they're properly exposed

## Related Data Types

The following data types remain active and are used by the commented-out functions:

- `data BCG a = BCG { _xBcg, _rBcg, _rHatBcg, _pBcg, _pHatBcg :: SpVector a }`
- `data CGS a = CGS { _x, _r, _p, _u :: SpVector a}`
- `data BICGSTAB a = BICGSTAB { _xBicgstab, _rBicgstab, _pBicgstab :: SpVector a}`

These data types could be commented out as well if they're truly unused, but they've been retained in case they're referenced elsewhere or needed for the future implementation of these solvers.
