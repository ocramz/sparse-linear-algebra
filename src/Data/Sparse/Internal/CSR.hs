{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses, CPP #-}
module Data.Sparse.Internal.CSR where

import Control.Applicative
import Control.Monad.Primitive
import Control.Monad.ST

import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Control.Monad
import Data.Maybe

import Data.Complex
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)


import Data.Sparse.Utils
import Data.Sparse.Types
import Data.Sparse.Internal.CSRVector
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
  CM {
      csrNrows :: {-# UNPACK #-} !Int,
      csrNcols :: {-# UNPACK #-} !Int,
      csrNz :: {-# UNPACK #-} !Int,
      csrColIx :: V.Vector Int,
      csrRowPtr :: V.Vector Int,
      csrVal :: V.Vector a} deriving Eq

instance Functor CsrMatrix where
  fmap f (CM m n nz cc rp x) = CM m n nz cc rp (fmap f x)

instance Foldable CsrMatrix where
  foldr f z cm = foldr f z (csrVal cm)


-- instance Set CsrMatrix where -- TODO: efficiency of intersection and union of sparse matrices depends on internal representation

instance Show a => Show (CsrMatrix a) where
  show m'@(CM m n nz cix rp x) = szs where
    szs = unwords ["CSR (",show m, "x", show n,"),",show nz, "NZ ( sparsity",show (spy m'),"), column indices:",show cix,", row pointers:", show rp,", data:",show x]




-- * Creation
-- | Copy a Vector containing (row index, column index, entry) into a CSR structure. Sorts the Vector by row indices ( O(log N) ), unzips column indices and data ( O(N) ) and generates the row pointer vector ( 2 O(N) passes )
toCSR :: Int -> Int -> V.Vector (Int, Int, a) -> CsrMatrix a
toCSR m n ijxv = CM m n nz cix crp x where
  nz = V.length x
  (rix, cix, x) = V.unzip3 (sortByRows ijxv)
  crp = csrPtrV m rix
  sortByRows = V.modify (VA.sortBy f) where
       f a b = compare (fst3 a) (fst3 b)









-- * Lookup
-- | O(1) : lookup row
lookupRow :: CsrMatrix a -> IxRow -> Maybe (CsrVector a)
lookupRow cm i | null er = Nothing
               | otherwise = Just er where er = extractRowCSR cm i

-- | O(N) lookup entry by index in a sparse Vector, using a default value in case of missing entry
lookupEntry_WD :: a -> CsrVector a -> Int -> a
lookupEntry_WD z cr j = maybe z (\j -> cvVal cr V.! j) (V.findIndex (== j) (cvIx cr))




-- | O(N) : lookup entry by (row, column) indices, using a default value in case of missing entry
lookupCSR_WD :: a -> CsrMatrix a -> (IxRow, IxCol) -> a
lookupCSR_WD z csr (i,j) = maybe z (\r -> lookupEntry_WD z r j) (lookupRow csr i)

-- | O(N) : lookup entry by (row, column) indices
lookupCSR :: a -> CsrMatrix a -> (IxRow, IxCol) -> Maybe a
lookupCSR z csr (i, j) = lookupRow csr i >>= \cr -> return $ lookupEntry_WD z cr j


-- ** Extract a row
-- | O(1) : extract a row from the CSR matrix. Returns an empty Vector if the row is not present.
extractRowCSR :: CsrMatrix a -> IxRow -> CsrVector a
extractRowCSR (CM _ n _ cix rp x) irow = CV n ixs vals where
  imin = rp V.! irow
  imax = rp V.! (irow + 1)
  ixs = trimf imin imax cix
  vals = trimf imin imax x
  trimf i1 i2 w = V.force $ V.drop i1 (V.take i2 w)



-- | O(N) : Rebuilds the (row, column, entry) Vector from the CSR representation. 
fromCSR :: CsrMatrix a -> V.Vector (Int, Int, a)
fromCSR mc = V.zip3 ii jj xx where (ii,jj,xx) = fromCSR0 mc

fromCSR0 :: CsrMatrix a -> (V.Vector Int, V.Vector Int, V.Vector a)
fromCSR0 mc = (rows, csrColIx mc, csrVal mc) where
  (m, n) = dim mc
  l = length (csrColIx mc)
  rp = csrRowPtr mc
  rows = V.create $ do
    rowv <- VM.replicate l 0
    forM_ [0 .. m-1] (\i -> go rowv i 0)
    return rowv
  go vm irow j | j <= nj - 1 = do
                          VM.write vm (j + jmin) irow
                          go vm irow (succ j)
               | otherwise = return () where
    jmin = rp V.! irow
    jmax = rp V.! (irow + 1)
    nj = jmax - jmin

-- | O(N log N) : Transpose CSR matrix
transposeCSR :: CsrMatrix a -> CsrMatrix a
transposeCSR mm = toCSR n m $ V.zip3 jj ii xx where
  (m,n) = dim mm
  (ii, jj, xx) = fromCSR0 mm





-- matVec mm v = fmap (\ri -> extractRow mm ri `dot` v) irows_ where
--   nrows = csrNrows mm
--   irows_ = [0 .. nrows - 1]





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





                                                 
           




  







-- * Utilities

tail3 :: (t, t1, t2) -> (t1, t2)
tail3 (_,j,x) = (j,x)

snd3 (_,j,_) = j

fst3 :: (t, t1, t2) -> t
fst3 (i, _, _) = i











-- playground


  


union :: Ord a => [a] -> [a] -> [a]
union u_ v_ = go u_ v_ where
  go [] x = x
  go y [] = y
  go uu@(u:us) vv@(v:vs)
    | u == v =    u : go us vs
    | u < v =     u : go us vv 
    | otherwise = v : go uu vs





  


    

safe :: (a -> Bool) -> (a -> b) -> a -> Maybe b
safe q f v
  | q v = Just (f v)
  | otherwise = Nothing

vectorNonEmpty :: (V.Vector a1 -> a) -> V.Vector a1 -> Maybe a
vectorNonEmpty = safe (not . V.null)

safeHead :: V.Vector a -> Maybe a
safeHead = vectorNonEmpty V.head

safeTail :: V.Vector a -> Maybe (V.Vector a)
safeTail = vectorNonEmpty V.tail













