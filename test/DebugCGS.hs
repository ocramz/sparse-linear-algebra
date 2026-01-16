{-# LANGUAGE FlexibleContexts, TypeFamilies #-}
module Main where

import Numeric.LinearAlgebra.Sparse
import Data.Sparse.Common

aa0 :: SpMatrix Double
aa0 = fromListDenseSM 2 [1,3,2,4]

b0 :: SpVector Double
b0 = mkSpVR 2 [8,18]

x0true :: SpVector Double
x0true = mkSpVR 2 [2,3]

main :: IO ()
main = do
  let x0 = fromListSV (dim b0) []  -- initial guess (zero vector)
      rhat = b0 ^-^ (aa0 #> x0)
      initState = cgsInit aa0 b0 x0
      
  putStrLn "=== CGS Debug Test ==="
  putStrLn "\nMatrix aa0:"
  prd aa0
  putStrLn "\nRHS b0:"
  prd b0
  putStrLn "\nTrue solution x0true:"
  prd x0true
  putStrLn "\n--- Check: aa0 #> x0true should equal b0 ---"
  prd (aa0 #> x0true)
  
  putStrLn "\n=== Initial state ==="
  putStrLn $ "x0 = " ++ show (_x initState)
  putStrLn $ "r0 = " ++ show (_r initState)
  putStrLn $ "p0 = " ++ show (_p initState)
  putStrLn $ "u0 = " ++ show (_u initState)
  putStrLn $ "||r0|| = " ++ show (norm2 (_r initState))
  
  putStrLn "\n=== CGS Iterations ==="
  let showState i s = do
        let residual = (aa0 #> _x s) ^-^ b0
            rnorm = norm2 residual
            xnorm = norm2 (_x s)
            xerr = norm2 (_x s ^-^ x0true)
        putStrLn $ "Iter " ++ show i ++ ":"
        putStrLn $ "  ||r|| = " ++ show rnorm
        putStrLn $ "  ||x|| = " ++ show xnorm
        putStrLn $ "  ||x-x*|| = " ++ show xerr
        putStrLn $ "  x = " ++ show (_x s)
  
  showState 0 initState
  
  let s1 = cgsStep aa0 rhat initState
  showState 1 s1
  
  let s10 = iterate (cgsStep aa0 rhat) initState !! 10
  showState 10 s10
  
  let s100 = iterate (cgsStep aa0 rhat) initState !! 100
  showState 100 s100
