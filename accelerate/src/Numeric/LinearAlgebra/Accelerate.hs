{-# language FlexibleContexts #-}
{-# language ScopedTypeVariables #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Accelerate
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  BSD3 (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- `accelerate` instances for linear algebra
--
-----------------------------------------------------------------------------
module Numeric.LinearAlgebra.Accelerate where


import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)))
import Data.Array.Accelerate.Interpreter (run)

-- | Vector as newtype
-- newtype Vector e = Vector (Array DIM1 e) deriving (Eq, Show)


-- | sparse vectors
newtype SVector e = SVector (Int, e) deriving (Eq, Show)

-- -- | segments : vector of segment lengths 
-- newtype Segments i = Segments (Vector i) deriving (Eq, Show)


data SMatrix i e = SMatrix {
    smNrows :: Int
  , smNcols :: Int
  , smSegments :: Segments i
  , smEntries :: SVector e } deriving (Eq, Show)



data SSMatrix i m n e = SSMatrix {
    ssmNrows :: m
  , ssmNcols :: n
  , ssmSegments :: Segments i
  , ssmEntries :: SVector e
                                 }





-- Sparse-matrix vector multiplication
-- -----------------------------------

type SparseVector e = Vector (A.Int32, e)
type SparseMatrix e = (Segments A.Int32, SparseVector e)


smvm :: A.Num a => Acc (SparseMatrix a) -> Acc (Vector a) -> Acc (Vector a)
smvm smat vec
  = let (segd, svec)    = A.unlift smat
        (inds, vals)    = A.unzip svec

        -- vecVals         = A.gather (A.map A.fromIntegral inds) vec
        vecVals         = A.backpermute
                             (A.shape inds)
                             (\i -> A.index1 $ A.fromIntegral $ inds A.! i)
                             vec
        products        = A.zipWith (*) vecVals vals
    in
    A.foldSeg (+) 0 products segd
