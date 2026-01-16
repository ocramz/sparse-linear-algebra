# Build and Test Instructions

This document provides instructions for building and testing the BiCGSTAB and linSolve0 re-enablement.

## Prerequisites

- Stack build tool installed
- Internet connection for downloading dependencies

## Build

```bash
# Build the project
stack build

# Or using make
make build
```

Expected output: Successful compilation with no errors.

## Test

### Run All Tests

```bash
# Run full test suite
stack test

# Or using make
make test
```

### Run Specific Test Suites

```bash
# Test BiCGSTAB only
stack test --test-arguments "-m BiCGSTAB"

# Test CGS only
stack test --test-arguments "-m CGS"

# Test linSolve0 interface
stack test --test-arguments "-m linSolve"
```

## Expected Test Results

### BiCGSTAB Tests
- ✅ bicgsInit creates initial BiCGSTAB state
- ✅ bicgstabStep performs one iteration
- ✅ BiCGSTAB converges on 2x2 system
- ✅ BiCGSTAB converges on 3x3 SPD system
- ✅ prop_bicgstab : BiCGSTAB converges for SPD systems (100 QuickCheck cases)

### linSolve0 Tests
- ✅ BiCGSTAB (2 x 2 dense)
- ✅ BiCGSTAB (3 x 3 sparse, symmetric pos.def.)
- ✅ CGS (2 x 2 dense)
- ✅ CGS (3 x 3 sparse, SPD)
- ✅ CGNE (2 x 2 dense)
- ✅ CGNE (3 x 3 sparse, SPD)

### CGS Tests (Existing, should still pass)
- ✅ cgsInit creates initial CGS state
- ✅ cgsStep performs one iteration
- ✅ CGS converges on 2x2 system
- ✅ CGS converges on 3x3 SPD system
- ✅ prop_cgs : CGS converges for SPD systems (100 QuickCheck cases)

## Troubleshooting

### If Tests Fail

1. **Check tolerance settings**: The tests use lenient tolerances (1e-4 relative, 1e-6 absolute)
   - If tests fail due to convergence issues, consider increasing the number of iterations
   - Location: `test/LibSpec.hs`, functions `checkBiCGSTAB`, `checkCGS`

2. **Check system properties**: Property tests guard against degenerate cases:
   - Systems smaller than 3x3
   - Nearly zero RHS or solution vectors
   - Very sparse matrices with few non-zeros
   - Large sparse systems (n > 20 with density < 0.1)

3. **Debug with trace output**: Use `checkCGSDebug` function for detailed iteration output

4. **Verify preconditioning**: Consider adding preconditioning for ill-conditioned systems

### Network Issues

If `stack build` fails due to network timeouts:

1. Try again - sometimes Stackage servers are temporarily unavailable
2. Check your internet connection
3. Try using a different resolver in `stack.yaml`
4. Use cached dependencies if available: `stack build --no-install-ghc`

## Verification Checklist

Before considering the task complete:

- [ ] `stack build` completes successfully with no warnings
- [ ] All BiCGSTAB tests pass
- [ ] All linSolve0 tests pass
- [ ] All existing tests (CGS, etc.) still pass
- [ ] No new compiler warnings introduced
- [ ] Documentation is up to date

## Performance Notes

- BiCGSTAB typically converges faster than CGS on most systems
- BiCGSTAB is more stable for ill-conditioned matrices
- linSolve0 uses 200 fixed iterations - convergence usually happens much sooner
- Property tests may take a few seconds due to QuickCheck generating 100 test cases

## Next Steps

After successful build and test:

1. Review test output for any warnings
2. Check if any tests are flaky or timeout
3. Consider adding preconditioner support
4. Benchmark performance against other solvers
5. Add more comprehensive documentation/examples
