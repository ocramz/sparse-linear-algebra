module Data.Sparse.Internal.CSR where

import Control.Monad.Primitive

import Data.Foldable (foldl')
import Data.List (group, groupBy)

import qualified Data.Vector as V 
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
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
      csrNz :: {-# UNPACK #-} !Int,
      csrColVal :: V.Vector (Int, a),
      csrRowPtr :: V.Vector Int } deriving Eq

instance Show a => Show (CsrMatrix a) where
  show (CsrMatrix m n nz cv rp) = szs where
    szs = unwords ["CSR (",show m, "x", show n,"),",show nz, "NZ"]

toCSR :: Int -> Int -> V.Vector (Int, Int, a) -> CsrMatrix a
toCSR m n ijxv = CsrMatrix m n nz cv crp where
  ijxv' = sortByRows ijxv -- merge sort over row indices ( O(log N) )
  cv = V.map tail3 ijxv'  -- map ( O(N) )
  crp = csrPtrVM m $ V.map fst3 ijxv'
  sortByRows = V.modify (VA.sortBy f) where
    f x y = compare (fst3 x) (fst3 y)
  nz = V.length ijxv'
  tail3 (_,j,x) = (j,x)
  fst3 (i, _, _) = i

-- csrRowPtr : for every row index in [0 .. m-1], count # NZ elements having matching row number
csrPtrVM :: Int -> V.Vector Int -> V.Vector Int
csrPtrVM nrows xs = V.scanl (+) 0 $ V.modify modf (V.replicate nrows 0) where
 modf vm = do
  let loop v ll i | i == nrows = return ()
                  | otherwise = do
        let lp = V.length $ V.takeWhile (== i) ll
        VM.write v i lp
        loop v (V.drop lp ll) (i + 1)
  loop vm xs 0  


-- csrPtr :: (Num t, Eq t) => t -> V.Vector t -> V.Vector Int
csrPtr nrows xs = scanl (+) 0 $ go xs 0 where
  go ll i | i == nrows = []
          | otherwise = lp : go (drop lp ll) (i + 1) where
              lp = length $ takeWhile (== i) ll
              
csrPtrV :: (Num t, Eq t) => t -> V.Vector t -> V.Vector Int
csrPtrV nrows xs = V.scanl' (+) 0 $ go xs 0 where
 go ll i | i == nrows = V.empty
         | otherwise = V.singleton lp V.++ go (V.drop lp ll) (i + 1) where
               lp = V.length $ V.takeWhile (== i) ll




data CsrVector a = CsrVector { cvDim :: {-# UNPACK #-} !Int,
                               cvIxVal :: V.Vector (Int, a) } deriving Eq

extractRow :: CsrMatrix a -> Int -> V.Vector (Int, a)
extractRow (CsrMatrix _ _ _ cv rp) irow = vals where
  imin = rp V.! irow
  imax = (rp V.! (irow + 1)) - 1
  vals = V.drop imin $ V.take (imax + 1) cv

rows (CsrMatrix m n nz cv crp) = foldl' insf [] ixd where
  ixd = V.zip crp (V.tail crp)
  insf acc (i1, i2) = V.toList (V.drop i1 (V.take (i2 + 1) cv)) : acc



m0 = toCSR 2 2 $ V.fromList [(0,0, pi), (1,0,3), (1,1,2)]
m1 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (1,0,2), (1,1,3), (2,0,4), (2,3,1), (3,2,2)]
m2 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (2,0,4), (2,3,1), (3,2,2)]
m3 = toCSR 4 4 $ V.fromList [(1,0,5), (1,1,8), (2,2,3), (3,1,6)]


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
