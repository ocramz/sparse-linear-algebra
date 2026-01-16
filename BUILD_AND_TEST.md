# Build and Test Instructions for CGS Re-enablement

## Overview
This document provides instructions for building and testing the re-enabled CGS (Conjugate Gradient Squared) solver functionality.

## Prerequisites
- Stack build tool installed
- Network connectivity for downloading dependencies

## Build Steps

### 1. Clean build (optional but recommended)
```bash
make clean
# or
rm -rf .stack-work/
```

### 2. Build the project
```bash
stack build
# or
make build
```

Expected outcome: Successful compilation with no errors.

### 3. Run the test suite
```bash
stack test
# or
make test
```

Expected outcome: All tests pass, including the new CGS tests:
- `cgsInit creates initial CGS state`
- `cgsStep performs one iteration`
- `CGS converges on 2x2 system`
- `CGS converges on 3x3 SPD system`
- `prop_cgs : CGS converges for SPD systems` (property test)

## Troubleshooting

### Compilation Errors
If you encounter compilation errors related to CGS:

1. **Missing imports**: Ensure `V` constraint is imported from `Numeric.LinearAlgebra.Class`
2. **Type errors**: Check that the type signatures match the expected constraints
3. **Module export errors**: Verify CGS items are properly exported

### Test Failures

#### Test: "cgsInit creates initial CGS state"
- **Issue**: Initial state not computed correctly
- **Debug**: Check that `r0 = b ^-^ (aa #> x0)` is correct
- **Fix**: Verify matrix-vector product and subtraction operations

#### Test: "CGS converges on systems"
- **Issue**: Algorithm doesn't converge within iteration limit
- **Debug**: 
  - Check if the test matrix is well-conditioned
  - Verify the algorithm implementation matches the standard CGS
  - Increase iteration count if needed (currently 50 for unit tests, 100 for property tests)
- **Possible fixes**:
  - Adjust tolerance in `nearZero` check
  - Increase maximum iterations
  - Add preconditioner (future enhancement)

#### Property Test: "prop_cgs"
- **Issue**: Random systems don't converge
- **Debug**: 
  - Property test uses SPD matrices: `m #^# m` for any matrix `m`
  - These should be well-conditioned
  - May fail for very ill-conditioned matrices
- **Fix**: May need to adjust the matrix generation or add condition number check

## Performance Notes

CGS typically converges faster than CG for many problems but can be less stable.
Expected behavior:
- Faster convergence than BiCG on well-conditioned SPD systems
- May exhibit irregular convergence behavior
- Residual may not decrease monotonically

## Verification Checklist

- [ ] Project builds without errors
- [ ] All existing tests still pass
- [ ] All new CGS unit tests pass
- [ ] Property tests pass consistently (run multiple times)
- [ ] No warnings related to CGS code
- [ ] Documentation is accurate

## Next Steps After Successful Build/Test

1. Run linter: `make lint`
2. Run security checks (CodeQL)
3. Request code review
4. Update README if needed
