{-# language FlexibleContexts, TypeFamilies #-}
module Numeric.LinearAlgebra.Core.Sparse where

import Data.List (unfoldr)
import Control.Monad.State.Strict

import Numeric.LinearAlgebra.Core.Class




-- | Givens' rotation matrix
rotMtx :: (a ~ ScElem m, Floating a, SparseMatrix m) => Int -> Int -> Int -> a -> m
rotMtx m ii jj angle = smFromList (m, m) arr
  where
    m' = m - 3 -- 2 on-diagonal values will be /= 1
    c = cos angle
    s = sin angle
    arr0 = [(i, i, 1) | i <- [0.. m']] -- identity matrix elements
    arr1 = [(jj, jj, c), (jj, ii, - s) -- 2D rotation matrix elements
           ,(ii, jj, s), (ii, ii, c) ]
    arr = arr1 ++ arr0



