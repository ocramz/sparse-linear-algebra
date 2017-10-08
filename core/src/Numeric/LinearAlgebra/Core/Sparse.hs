{-# language FlexibleContexts, TypeFamilies #-}
module Numeric.LinearAlgebra.Core.Sparse where

import Data.List (unfoldr)
import Control.Monad.State.Strict

import Numeric.LinearAlgebra.Core.Class




-- | Givens' rotation matrix
rotMtx :: (a ~ SpElem m, Floating a, SparseMatrix m) => Int -> Int -> Int -> a -> m
rotMtx m ii jj angle = smFromList (m, m) arr
  where
    m' = m - 3 -- 2 on-diagonal values will be /= 1
    c = cos angle
    s = sin angle
    arr0 = [(i, i, 1) | i <- [0.. m']] -- identity matrix elements
    arr1 = [(jj, jj, c), (jj, ii, - s) -- 2D rotation matrix elements
           ,(ii, jj, s), (ii, ii, c) ]
    arr = arr1 ++ arr0






-- | Test sparse matrix 

data SpM a = SpM { spmDims :: (Int, Int), spmNnz :: Int, spmData :: [(Int, Int, a)] } deriving (Eq, Show)

instance FiniteDim (SpM a) where
  type FDSize (SpM a) = (Int, Int)
  dim = spmDims

instance Sparse (SpM a) where
  -- type SpIx (SpM a) = (Int, Int)
  type SpElem (SpM a) = a
  nnz = spmNnz
  dens c = fromIntegral (nnz c) / fromIntegral (m * n) where (m, n) = dim c


instance SparseMatrix (SpM a) where
  smFromList dims ll = SpM dims (length ll) ll
  smToList (SpM _ _ ll) = ll



m0 = SpM (2,2) 2 [(0,0, exp 1), (1,1, pi)]
