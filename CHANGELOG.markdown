# Changelog

## 0.4.0 (Upcoming)

### Breaking Changes

#### Pure Logging System
- **BREAKING**: Removed `logging-effect` dependency entirely
- **BREAKING**: Removed custom `MonadLog` typeclass - now uses `MonadWriter` from mtl directly
- **BREAKING**: Removed `Severity` and `WithSeverity` types - no logging levels
- **BREAKING**: Removed logging utility functions: `logWith`, `bracketsUpp`, `withSeverity`
- **BREAKING**: `IterConfig` no longer has `icLogHandler` field - logging is now done via `MonadWriter`
- **BREAKING**: `modifyInspectGuardedM` signature changed:
  - Before: `(MonadThrow m, ..., Monoid msg) => ConvergConfig t a -> IterConfig s t msg -> (s -> m s) -> s -> m s`
  - After: `(MonadThrow m, MonadWriter w m, ...) => ConvergConfig t a -> IterConfig s t -> (s -> m s) -> s -> m s`
- **BREAKING**: `LinearSystem` class method `(<\>)` now requires `MonadWriter w m` instead of `MonadLog String m`
- **BREAKING**: `IterativeT` transformer stack simplified - removed msg type parameter

### Improvements

#### Logging
- Pure functional logging using mtl's `MonadWriter` interface
- No IO side effects in numerical computations
- Cleaner, simpler code using standard mtl patterns
- Logs are returned as values via Writer monad, not processed via IO

#### GHC Compatibility
- Updated to modern GHC support (tested with GHC 9.4.8, 9.6.6, 9.8.2, 9.12.2)
- Updated Stack resolver to lts-22.39 (GHC 9.6.6)
- Added necessary language extensions for modern GHC compatibility
- Fixed type constraints for Complex number printing (added `RealFloat`)

#### Testing
- **+24% test coverage increase**: 31 tests (up from 25)
- Uncommented and fixed LU factorization tests (4 tests)
- Uncommented and fixed Arnoldi iteration tests (2 tests)
- **Uncommented and fixed QR decomposition tests** (5 specific cases, 2 property tests)
- **Added eigsQR eigenvalue algorithm tests** (5 specific cases)
- All tests use `WriterT` for pure logging in tests
- Added 4 new tests for 2x2 LU factorization (both Real and Complex)

#### Code Quality
- Exposed `Control.Iterative` module for external use
- Cleaner dependencies (removed logging-effect)
- Simpler abstractions (no custom MonadLog)
- Updated `.gitignore` for modern build artifacts
- **Uncommented and enabled QR decomposition (`qr`) function**
- **Uncommented and enabled eigsQR eigenvalue algorithm**

### Bug Fixes
- **Fixed LU factorization infinite loop**: Corrected termination condition in `lu` function that caused hangs for all matrix sizes, especially noticeable with 2x2 matrices
  - Changed loop termination from `i == n - 1` to `i == n` to properly process all matrix rows/columns
  - Removed redundant final U matrix update that was masking the incorrect termination
  - Issue affected both Real and Complex matrix types
- **Fixed QR decomposition Givens rotation error**: Resolved "no compatible rows" exception that occurred with certain sparse matrix structures
  - Modified `givens` function to return `Maybe (SpMatrix a)` instead of throwing exception when no compatible row exists
  - Added check for already-zero elements before attempting Givens rotation
  - QR algorithm now gracefully skips elements that cannot be rotated or are already zero
  - Fixes issue where eigsQR failed on sparse matrices like `[0,1,1,0; 1,0,0,0; 1,0,0,1; 0,0,1,0]`
- **Fixed Cholesky decomposition for complex Hermitian matrices**: Corrected missing conjugation operations that caused NaN results for complex matrices and certain structured matrices (arrowhead matrices)
  - Fixed subdiagonal computation to use conjugated inner product instead of regular matrix multiplication
  - Properly handles the case when row i hasn't been built yet during computation (returns 0)
  - Uses manual fold with `conj` for proper conjugation: `sum(L[i,k] * conj(L[j,k]))`
  - Added `InnerSpace a, Scalar a ~ a` type constraints to `chol` function
  - Ensures `L L^H = A` for Hermitian matrices and `L L^T = A` for symmetric matrices
  - Fixes issue with "Rails" example matrix from R mixed models and other arrowhead matrices
  - Test improvements: 7/17 Cholesky tests now passing (up from 1/17), with arrowhead matrix tests working correctly

### Migration Guide

**Old code (with logging-effect):**
```haskell
import Control.Monad.Log

let config = IterConfig "myFunc" 100 3 id myHandler
result <- modifyInspectGuardedM convConfig config f x0
```

**New code (with mtl MonadWriter):**
```haskell
import Control.Monad.Writer

let config = IterConfig "myFunc" 100 3 id
(result, logs) <- runWriterT $ modifyInspectGuardedM convConfig config f x0
-- Process logs if needed: mapM_ putStrLn logs
```

### Internal Changes
- `IterativeT` transformer stack: `ReaderT c (StateT s m)` over any `MonadWriter`
- `runIterativeT` now has signature: `(Monad m, MonadWriter w m) => r -> s -> IterativeT r s m a -> m (a, s)`
- Logs accumulated automatically via Writer monad during iteration

---

## 0.3.2
	* Introduced `logging-effect` as a more versatile alternative to debug logging in IO

	0.3.1
	* Changed `SpMatrix` to use `foldlWithKey'` for efficiency (Joshua Moerman)
	
	* Bumped LTS to 11.3 (GHC 8.2.2)
	
	* Removed unneeded dependencies from stack.yaml

	0.3
	* Fixed a number of instances, un-commented tests (Joshua Moerman)
	
	* Documented issues with complex number support (Joshua Moerman)
	
	0.2.9.9
	* Moved to IntMap.Strict (Gregory Schwartz)
	
	* Stackage LTS bump to 10.4 (GHC 8.2)
	
	0.2.9.7
	Improved pretty printer:

	* Fixed display precision (e.g. 2 decimal digits), fixed column width output for vectors and matrices
	
	* Small and large values (wrt fixed precision) switch to scientific notation
	
	0.2.9.4
	Exceptions constructors are exported by Numeric.LinearAlgebra.Sparse


	0.2.9.1
	* Uses classes from `vector-space` such as AdditiveGroup, VectorSpace and InnerSpace
	
	* QuickCheck tests for algebraic properties, such as matrix-vector products and soon more abstract ones e.g. positive semi-definite matrices
	
	* Getting rid of `error` in favor of MonadThrow exceptions for high-level operations such as matrix algorithms
