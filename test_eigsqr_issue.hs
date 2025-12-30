{-# LANGUAGE FlexibleContexts #-}
-- Test script to verify the eigsQR issue is fixed
import Numeric.LinearAlgebra.Sparse
import Data.Sparse.Common

main :: IO ()
main = do
  let mat = sparsifySM $ fromListDenseSM 4 [0,1,1,0, 1,0,0,0, 1,0,0,1, 0,0,1,0] :: SpMatrix Double
  
  putStrLn "Testing eigsQR with the issue matrix..."
  putStrLn "Matrix:"
  prd mat
  
  putStrLn "\nRunning eigsQR 100 False mat..."
  result <- eigsQR 100 False mat
  
  putStrLn "\nEigenvalues (approximate):"
  prd result
  
  putStrLn "\nTest passed! No 'Givens : no compatible rows' error."
  
  -- Expected eigenvalues from R: 1.618034, 0.618034, -0.618034, -1.618034
  putStrLn "\nExpected eigenvalues (from R):"
  putStrLn "[1.618034, 0.618034, -0.618034, -1.618034]"
