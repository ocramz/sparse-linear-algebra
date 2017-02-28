{-# language DeriveFunctor #-}
module Data.Sparse.Internal.CSB where

import Control.Applicative
-- import Control.Monad.Primitive
import Control.Monad.ST

import Data.Foldable (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import qualified Data.IntMap as IM

import Control.Monad
import Data.Maybe
import Data.Ord (comparing)

-- import Data.Complex
-- import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)

import Control.Concurrent
-- import qualified Control.Monad.Par as Par
-- import Control.DeepSeq


import Data.Sparse.Utils
import Data.Sparse.Types
-- import Data.Sparse.Internal.CSRVector
import Data.Sparse.Internal.Utils

import Numeric.LinearAlgebra.Class


{-
Specification of the "Compressed Sparse Block" matrix storage format, as published in [1]:

Let f(i, j) be the bijection from pairs of block indices to integers in the range 0,1,..., n^2/beta^2 − 1 that describes the ordering among blocks.
That is, f(i,j) < f(i',j') if and only if Aij appears before Ai'j' in val.

The row_ind and col_ind arrays store the row and column indices,
respectively, of the elements in the val array. These indices
are relative to the block containing the particular element, not the
entire matrix, and hence they range from 0 to β−1.

If val[k] stores the matrix element aiβ+r,jβ+c, which is located in the rth row
and cth column of the block Aij, then row_ind = r and col_ind = c.

The blk_ptr array stores the index of each block in the val array, which is analogous to the row_ptr array for CSR. If val[k] stores a matrix element falling in block Aij, then blk_ptr[ f(i,j) ] ≤ k < blk_ptr[ f(i,j)+1 ].

[1] A. Buluc, et al., Parallel Sparse Matrix-Vector and Matrix-Transpose-Vector Multiplication Using Compressed Sparse Blocks
-}


-- | CSB Block
-- 
-- Invariants :
-- 1) rowIx, colIx and val have same length
-- 2) " have at most (bRows x bCols) NZ
--
data CsbMatrix a = CSB {
  csbNrows, csbNcols :: {-# UNPACK #-} !Int,
  csbVal :: V.Vector a,
  csbBlkPtr, csbRowIx, csbColIx :: V.Vector Int
                       } deriving (Eq, Functor)


csbParams :: (Int, Int)  -- ^ Matrix size
          -> Int         -- ^ Block parameter (i.e. block edge length)
          -> (Int, Int)  -- ^ # of blocks per side
csbParams (m,n) beta = (ceiling m', ceiling n')
  where m' = fromIntegral m / fromIntegral beta
        n' = fromIntegral n / fromIntegral beta

-- | Which block (in row, col format) do matrix indices (i,j) fall in ?
bBin, bCoords :: Integral t => t -> t -> t -> (t, t)
bBin beta i j = (i `div` beta, j `div` beta)

-- | Block coordinates
bCoords beta i j = (i `mod` beta, j `mod` beta)

-- | Block index (row-major order)
blockIx :: (Int, Int) -> Int -> Int -> Int -> Int
blockIx dims beta i j = bx + by*nbx where
  (bx, by) = bBin beta i j
  (nbx, _) = csbParams dims beta



type Block i a = [(i,i,a)]

{-# inline consBlockElem #-}
consBlockElem ::
  IM.Key -> (i, i, a) -> IM.IntMap (Block i a) -> IM.IntMap (Block i a)
consBlockElem i x bb = IM.insert i (x : blocki) bb where
  blocki = fromMaybe [] (IM.lookup i bb)


consBlocks :: Foldable t =>
     (Int, Int) -> Int -> t (Int, Int, a) -> IM.IntMap (Block Int a)
consBlocks dims beta ijx = foldl' ins IM.empty ijx where
  ins accb (i,j,x) = consBlockElem ib (i,j,x) accb where
    ib = blockIx dims beta i j 
    



