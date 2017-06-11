{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}
module Data.Sparse.Internal.CSC where

import Control.Monad (forM_, when)

import qualified Data.Vector as V
import qualified Data.Vector.Mutable as VM

import Data.Sparse.Types
-- import Data.Sparse.Internal.CSRVector
import Data.Sparse.Internal.Utils
import Numeric.LinearAlgebra.Class


data CscMatrix a =
  CscM {
      cscNrows :: {-# UNPACK #-} !Int,
      cscNcols :: {-# UNPACK #-} !Int,
      cscNz :: {-# UNPACK #-} !Int,
      cscRowIx :: V.Vector Int,
      cscColPtr :: V.Vector Int,
      cscVal :: V.Vector a} deriving Eq

instance Functor CscMatrix where
  fmap f (CscM m n nz cc rp x) = CscM m n nz cc rp (fmap f x)

instance Foldable CscMatrix where
  foldr f z cm = foldr f z (cscVal cm)


instance Show a => Show (CscMatrix a) where
  show m'@(CscM m n nz cix rp x) = szs where
    szs = unwords ["CSC (",show m, "x", show n,"),",show nz, "NZ ( sparsity",show (spy m'),"), row indices:",show cix,", column pointers:", show rp,", data:",show x]


-- some instances

instance FiniteDim CscMatrix where
  type FDSize CscMatrix = (Int, Int)
  dim m = (cscNrows m, cscNcols m)

instance HasData CscMatrix a where
  nnz = cscNz
  
instance Sparse CscMatrix a where
  spy m = fromIntegral (nnz m) / fromIntegral (cscNrows m * cscNcols m)

  


-- * Creation
-- | O(N log N) : Copy a Vector containing (row index, column index, entry) into a CSC structure. Sorts the Vector by row columns ( O(log N) ), unzips row indices and data ( O(N) ) and generates the column pointer vector ( 2 O(N) passes )
toCSC :: Int -> Int -> V.Vector (Int, Int, a) -> CscMatrix a
toCSC m n ijxv = CscM m n nz rix crp x where
  nz = V.length x
  (rix, cix, x) = V.unzip3 (sortWith snd3 ijxv)  -- sort by columns
  crp = csPtrV (==) m cix

-- | O(N) : Rebuilds the (row, column, entry) Vector from the CSC representation. 
fromCSC :: CscMatrix a -> V.Vector (Int, Int, a)
fromCSC mc = V.zip3 ii jj xx where (ii,jj,xx) = fromCSC0 mc
fromCSC0 :: CscMatrix a -> (V.Vector Int, V.Vector Int, V.Vector a)
fromCSC0 mc = (rowIx, cols, cscVal mc) where
  (m, n) = dim mc
  rowIx = cscRowIx mc
  l = length rowIx
  cp = cscColPtr mc
  cols = V.create $ do
    rowv <- VM.replicate l 0
    forM_ [0 .. m-1] (\i -> go rowv i 0)
    return rowv
  go vm irow j = when (j <= nj - 1) $ do
                          VM.write vm (j + jmin) irow
                          go vm irow (succ j) where
    jmin = cp V.! irow
    jmax = cp V.! (irow + 1)
    nj = jmax - jmin





-- example data

-- row = np.array([0, 0, 1, 2, 2, 2])
-- col = np.array([0, 2, 2, 0, 1, 2])
-- data = np.array([1, 2, 3, 4, 5, 6])
