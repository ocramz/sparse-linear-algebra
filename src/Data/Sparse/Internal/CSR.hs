module Data.Sparse.Internal.CSR where

import Control.Monad.Primitive

import Data.List (group, groupBy)

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
-- import qualified Data.Vector.Unboxed.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA (sort, sortBy)
import qualified Data.Vector.Generic as VG (convert)

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
  CsrMatrix {
      csrNrows :: {-# UNPACK #-} !Int,
      csrNcols :: {-# UNPACK #-} !Int,
      csrNnz :: {-# UNPACK #-} !Int,
      csrColVal :: V.Vector (Int, a),
      csrRowPtr :: V.Vector Int } deriving Eq

instance Show a => Show (CsrMatrix a) where
  show (CsrMatrix m n nz cv rp) = unwords ["CSR", szs, ",",nzs] where
    szs = unwords ["(",show m, "x", show n,")"]
    nzs = unwords [show nz, " NZ"]

toCSR :: Int -> Int -> V.Vector (Int, Int, a) -> CsrMatrix a
toCSR m n ijxv = CsrMatrix m n nz cv crp where
  ijxv' = sortByRows ijxv -- merge sort over row indices ( O(log N) )
  ll = V.fromList .
       groupBy (\x y -> fst3 x == fst3 y) . -- groupBy ( O(N) )
       V.toList $ ijxv'
  cv = V.map tail3 ijxv'  -- map ( O(N) )
  crp = V.scanl (\x vs -> x + length vs) 0 ll -- scanl * length
  nz = V.length ijxv'
  sortByRows = V.modify (VA.sortBy f) where
    f x y = compare (fst3 x) (fst3 y)
  tail3 (_,j,x) = (j,x)
  fst3 (i, _, _) = i

extractRow :: CsrMatrix a -> Int -> V.Vector a
extractRow (CsrMatrix _ _ _ cv rp) irow = V.map snd vals where
  imin = rp V.! irow
  imax = (rp V.! (irow + 1)) - 1
  vals = V.drop imin $ V.take (imax + 1) cv





m0 = toCSR 2 2 $ V.fromList [(0,0, pi), (1,0,3), (1,1,2)]



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
