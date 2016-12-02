module Data.Sparse.Internal.CSB where

import Control.Applicative
import Control.Monad.Primitive
import Control.Monad.ST

import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Control.Monad
import Data.Maybe

import Data.Complex
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)

import Control.Concurrent
-- import qualified Control.Monad.Par as Par
-- import Control.DeepSeq


import Data.Sparse.Utils
import Data.Sparse.Types
import Data.Sparse.Internal.CSRVector
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


-- * CSB Block

-- | Invariants :
-- 1) rowIx, colIx and val have same length
-- 2) " have at most (bRows x bCols) NZ

data Block a = Block {
  bRows, bCols :: {-# UNPACK #-} !Int,
  rowIx :: V.Vector Int,
  colIx :: V.Vector Int,
  val :: V.Vector a      -- invariant: 
  } deriving (Eq, Show)

instance Functor Block where
  fmap f (Block br bc ri ci x) = Block br bc ri ci (fmap f x)

instance Foldable Block where
  foldr f z bl = foldr f z (val bl)

instance Traversable Block where
  traverse f (Block br bc ri ci v) = Block br bc ri ci <$> traverse f v


sortWith :: Ord o => (t -> o) -> V.Vector t -> V.Vector t
sortWith f v = V.modify (VA.sortBy q) v where
  q x y = compare (f x) (f y)

sortWith3 :: Ord o =>
   ((a, b, c) -> o) ->
   V.Vector a ->
   V.Vector b ->
   V.Vector c ->
   V.Vector (a, b, c)
sortWith3 f x y z = sortWith f $ V.zip3 x y z


sortByFst3 :: Ord a => V.Vector a -> V.Vector b -> V.Vector c -> V.Vector (a, b, c)
sortByFst3 = sortWith3 fst3

sortBySnd3 :: Ord b => V.Vector a -> V.Vector b -> V.Vector c -> V.Vector (a, b, c)
sortBySnd3 = sortWith3 snd3







-- * CSB Matrix

-- | Invariants:
-- 1) csbNrows = sum bRows, csbNcols = sum bCols

data CsbMatrix a =
  CB { csbNrows :: {-# UNPACK #-} !Int,
       csbNcols :: {-# UNPACK #-} !Int,
       -- csbBlkPtr :: V.Vector Int,
       csbBlocks :: V.Vector (Block a) -- block ordering function f(i,j)
     } deriving Eq

instance Functor CsbMatrix where
  fmap f (CB m n bb) = CB m n $ fmap (fmap f) bb

-- instance Foldable CsbMatrix where
--   foldr f z cm = foldr f z (csbBlocks cm)

-- instance Traversable CsbMatrix where
--   traverse f (CB m n bp x) = CB m n bp <$> traverse f x





-- * fromList
-- m n [(i, j, x)] -> C






-- * Mat-vec (y = A x)
-- y[i] = < A[i,:] , x >


-- * MatT-vec (y = A^T x)
-- y[i] = < A[:, j] , x>


-- * Mat-vec (y = A x, parallel over CSB and block-partitioned x)
-- | Invariants:
-- 1) # cols block Aij == # entries j'th partition of x


