module Data.Array.Accelerate.Sparse.Static where

import Data.Array.Accelerate (Segments)

import Data.Array.Accelerate.Sparse.SVector


-- | statically-typed dimensions. `m` and `n` could be Nat
data SSMatrix i m n e = SSMatrix {
    ssmNrows :: m
  , ssmNcols :: n
  , ssmSegments :: Segments i
  , ssmEntries :: SVector e
                                 }
