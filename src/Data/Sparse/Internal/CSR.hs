module Data.Sparse.Internal.CSR where

import Control.Applicative
import Control.Monad.Primitive
import Control.Monad.ST

import Data.Foldable (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Data.Sparse.Utils

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
  crp = csrPtrVM m $ V.map fst3 ijxv' -- scanl + replicate + takeWhile * map ( O(N) )
  sortByRows = V.modify (VA.sortBy f) where
       f x y = compare (fst3 x) (fst3 y)
  nz = V.length ijxv'
  csrPtrVM nrows xs = V.scanl (+) 0 $ V.modify modf (V.replicate nrows 0) where
   modf vm = do
     let loop v ll i | i == nrows = return ()
                     | otherwise = do
                                     let lp = V.length $ V.takeWhile (== i) ll
                                     VM.write v i lp
                                     loop v (V.drop lp ll) (i + 1)
     loop vm xs 0  


data CsrVector a = CsrVector { cvDim :: {-# UNPACK #-} !Int,
                               cvIxVal :: V.Vector (Int, a) } deriving Eq

extractRow :: CsrMatrix a -> Int -> V.Vector (Int, a)
extractRow (CsrMatrix m n _ cv rp) irow = vals where
  imin = rp V.! irow
  imax = (rp V.! (irow + 1)) - 1
  vals = V.drop imin $ V.take (imax + 1) cv

fromCSR :: CsrMatrix a -> V.Vector (Int, Int, a)
fromCSR mc = mconcat $ map (\i -> withRowIx (extractRow mc i) i) [0 .. csrNrows mc - 1] where
  withRowIx v i = V.zip3 (V.replicate n i) icol_ v_ where
    (icol_, v_) = V.unzip v
    n = V.length v_






-- * Sparse vector

fromDenseV :: Int -> V.Vector a -> CsrVector a
fromDenseV n xs = CsrVector n (V.indexed xs)


v0,v1 :: V.Vector Int
v0 = V.fromList [0,1,2,5,6]
v1 = V.fromList [0,2,3,4,6]

-- | intersection between sorted vectors, in-place updates
isectvM :: (Ord a, Num a) => V.Vector a -> V.Vector a -> V.Vector a
isectvM u_ v_ = V.force $ V.modify (modf u_ v_) (V.replicate n 0) where
  n = max (V.length u_) (V.length v_)
  modf u_ v_ vm = do
    let go u_ v_ i vm | V.null u_ || V.null v_ || i == n = return ()
                      | otherwise =  do
         let (u,us) = (V.head u_, V.tail u_)
             (v,vs) = (V.head v_, V.tail v_)
         if u == v then do VM.write vm i u
                           go us vs (i + 1) vm
                   else if u < v then go us v_ i vm
                                 else go u_ vs i vm
    go u_ v_ 0 vm

isectvM' :: (PrimMonad m, Ord a) => Int -> V.Vector a -> V.Vector a -> m (V.Vector a)
isectvM' n u_ v_ = do
    vm <- VM.new n
    let go u_ v_ i vm | V.null u_ || V.null v_ || i == n = return (vm, i)
                      | otherwise =  do
         let (u,us) = (V.head u_, V.tail u_)
             (v,vs) = (V.head v_, V.tail v_)

         if u == v then do VM.write vm i u
                           go us vs (i + 1) vm
                   else if u < v then go us v_ i vm
                                 else go u_ vs i vm
    (vm', i') <- go u_ v_ 0 vm
    V.take i' <$> V.freeze vm'
 


-- | intersection between sorted Vectors

isectv' u_ v_ = V.create $ do
  let n = max (V.length u_) (V.length v_)
  vm <- VM.new n
  let go u_ v_ i vm | V.null u_ || V.null v_ || i == n = return (vm, i)
                      | otherwise =  do
         let (u,us) = (V.head u_, V.tail u_)
             (v,vs) = (V.head v_, V.tail v_)
         if u == v then do VM.write vm i u
                           go us vs (i + 1) vm
                   else if u < v then go us v_ i vm
                                 else go u_ vs i vm
  (vm', i') <- go u_ v_ 0 vm
  let vm'' = VM.take i' vm'
  return vm''

isectv :: Ord a => V.Vector a -> V.Vector a -> V.Vector a
isectv = isectvWith compare

isectvWith :: (a -> t -> Ordering) -> V.Vector a -> V.Vector t -> V.Vector a
isectvWith cmp u_ v_ | V.null u_ || V.null v_ = V.empty
                     | otherwise = 
  let (u,us) = (V.head u_, V.tail u_)
      (v,vs) = (V.head v_, V.tail v_)
  in case cmp u v of EQ -> V.cons u (isectvWith cmp us vs)
                     LT -> isectvWith cmp us v_
                     GT -> isectvWith cmp u_ vs


-- | intersection between _sorted_ lists (see @isect@ in `data-ordlist`)
intersect :: Ord a => [a] -> [a] -> [a]
intersect (x:xs) (y:ys) | x == y = x : xs `intersect` ys
                        | otherwise = if x < y
                                      then xs `intersect` (y:ys)
                                      else (x:xs) `intersect` ys
intersect _ [] = []
intersect [] _ = []








-- * Utilities

tail3 :: (t, t1, t2) -> (t1, t2)
tail3 (_,j,x) = (j,x)

fst3 :: (t, t1, t2) -> t
fst3 (i, _, _) = i

-- test data

m0,m1,m2,m3 :: CsrMatrix Double
m0 = toCSR 2 2 $ V.fromList [(0,0, pi), (1,0,3), (1,1,2)]
m1 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (1,0,2), (1,1,3), (2,0,4), (2,3,1), (3,2,2)]
m2 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (2,0,4), (2,3,1), (3,2,2)]
m3 = toCSR 4 4 $ V.fromList [(1,0,5), (1,1,8), (2,2,3), (3,1,6)]



-- playground

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
