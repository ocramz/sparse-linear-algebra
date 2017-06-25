{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses, CPP #-}
module Data.Sparse.Internal.CSR where

-- import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
-- import qualified Data.Vector.Generic as VG (convert)

import Control.Monad (when, forM_)

import Data.Sparse.Types
import Data.Sparse.Internal.SVector
import Data.Sparse.Internal.Utils

import Numeric.LinearAlgebra.Class
-- import Data.Sparse.Common


{-| Compressed Row Storage specification :

   http://netlib.org/utk/people/JackDongarra/etemplates/node373.html

   The compressed row storage (CRS) format puts the subsequent nonzeros of the matrix
   rows in contiguous memory locations. Assuming we have a nonsymmetric sparse matrix
   $A$, we create three vectors: one for floating point numbers (val) and the other
   two for integers (col_ind, row_ptr).

   The val vector stores the values of the nonzero elements of the matrix $A$ as
   they are traversed in a row-wise fashion.
 
   The col_ind vector stores the column indexes of the elements in the val vector,
   that is, if val(k)=a_{i,j}, then  col_ind(k)=j$.

   The row_ptr vector stores the locations in the val vector that start a row;
   that is, if  val(k)=a_{i,j}, then row_ptr(i) <= k < row_ptr(i+1)

-}

data CsrMatrix a =
  CsrM {
      csrNrows :: {-# UNPACK #-} !Int,
      csrNcols :: {-# UNPACK #-} !Int,
      csrNz :: {-# UNPACK #-} !Int,
      csrColIx :: V.Vector Int,
      csrRowPtr :: V.Vector Int,
      csrVal :: V.Vector a} deriving Eq

instance Functor CsrMatrix where
  fmap f (CsrM m n nz cc rp x) = CsrM m n nz cc rp (fmap f x)

instance Foldable CsrMatrix where
  foldr f z cm = foldr f z (csrVal cm)


-- instance Set CsrMatrix where -- TODO: efficiency of intersection and union of sparse matrices depends on internal representation

instance Show a => Show (CsrMatrix a) where
  show m'@(CsrM m n nz cix rp x) = szs where
    szs = unwords ["CSR (",show m, "x", show n,"),",show nz, "NZ ( sparsity",show (spy m'),"), column indices:",show cix,", row pointers:", show rp,", data:",show x]




-- * Creation
-- | O(N log N) : Copy a Vector containing (row index, column index, entry) into a CSR structure. Sorts the Vector by row indices ( O(log N) ), unzips column indices and data ( O(N) ) and generates the row pointer vector ( 2 O(N) passes )
toCSR :: Int -> Int -> V.Vector (Int, Int, a) -> CsrMatrix a
toCSR m n ijxv = CsrM m n nz cix crp x where
  nz = V.length x
  (rix, cix, x) = V.unzip3 (sortWith fst3 ijxv)  -- sort by rows
  crp = csPtrV (==) m rix









-- * Lookup
-- | O(1) : lookup row
lookupRow :: CsrMatrix a -> IxRow -> Maybe (SVector a)
lookupRow cm i | null er = Nothing
               | otherwise = Just er where er = extractRowCSR cm i

-- | O(N) lookup entry by index in a sparse Vector, using a default value in case of missing entry
lookupEntry_WD :: a -> SVector a -> Int -> a
lookupEntry_WD z cr j = maybe z (\j -> svVal cr V.! j) (V.findIndex (== j) (svIx cr))




-- | O(N) : lookup entry by (row, column) indices, using a default value in case of missing entry
lookupCSR_WD :: a -> CsrMatrix a -> (IxRow, IxCol) -> a
lookupCSR_WD z csr (i,j) = maybe z (\r -> lookupEntry_WD z r j) (lookupRow csr i)

-- | O(N) : lookup entry by (row, column) indices
lookupCSR :: a -> CsrMatrix a -> (IxRow, IxCol) -> Maybe a
lookupCSR z csr (i, j) = lookupRow csr i >>= \cr -> return $ lookupEntry_WD z cr j


-- ** Extract a row
-- | O(1) : extract a row from the CSR matrix. Returns an empty Vector if the row is not present.
extractRowCSR :: CsrMatrix a -> IxRow -> SVector a
extractRowCSR (CsrM _ n _ cix rp x) irow = SV n ixs vals where
  imin = rp V.! irow
  imax = rp V.! (irow + 1)
  len = imax - imin
  trimf  = V.slice imin len
  ixs = trimf cix
  vals = trimf x



-- | O(N) : Rebuilds the (row, column, entry) Vector from the CSR representation. 
fromCSR :: CsrMatrix a -> V.Vector (Int, Int, a)
fromCSR mc = V.zip3 ii jj xx where (ii,jj,xx) = fromCSR0 mc
fromCSR0 :: CsrMatrix a -> (V.Vector Int, V.Vector Int, V.Vector a)
fromCSR0 mc = (rows, csrColIx mc, csrVal mc) where
  (m, _) = dim mc
  l = length (csrColIx mc)
  rp = csrRowPtr mc
  rows = V.create $ do
    rowv <- VM.replicate l 0
    forM_ [0 .. m-1] (\i -> go rowv i 0)
    return rowv
  go vm irow j = when (j < nj) $ do
                          VM.write vm (j + jmin) irow
                          go vm irow (succ j) where
    jmin = rp V.! irow
    jmax = rp V.! (irow + 1)
    nj = jmax - jmin

-- | O(N log N) : Transpose CSR matrix
transposeCSR :: CsrMatrix a -> CsrMatrix a
transposeCSR mm = toCSR n m $ V.zip3 jj ii xx where
  (m,n) = dim mm
  (ii, jj, xx) = fromCSR0 mm




-- some instances

instance FiniteDim CsrMatrix where
  type FDSize CsrMatrix = (Int, Int)
  dim m = (csrNrows m, csrNcols m)

instance HasData CsrMatrix a where
  nnz = csrNz
  
instance Sparse CsrMatrix a where
  spy m = fromIntegral (nnz m) / fromIntegral (csrNrows m * csrNcols m) 

instance Elt a => SpContainer CsrMatrix a where
  type ScIx CsrMatrix = (Int, Int)
  -- scInsert = undefined
  (@@) = lookupCSR_WD 0

-- instance Elt a => SparseMatrix CsrMatrix a where








-- test data
v0 :: V.Vector (Int, Double)
v0 = V.fromList [(0, pi), (2, 3), (3, 57)]

v1 = V.fromList [(0,0,pi), (0,1,5), (2,1,exp 1)]

-- sm0 = smFromVector ColsFirst (3, 3) v1 :: SpMatrix1 Double




v2r, v2c :: V.Vector (LexIx, Double)
v2r = V.fromList [(1, 5), (5, 8), (7, 6), (10, 3)]

v2c = V.fromList [(4, 5), (5, 8), (10, 3), (13, 6)]

v2 :: V.Vector (Int, Int, Double)
v2 = V.fromList [(1,0,5), (1,1,8), (2,2,3), (3,1,6)]
-- sm2 = smFromVector RowsFirst (4, 4) v2 :: SpMatrix1 Double
