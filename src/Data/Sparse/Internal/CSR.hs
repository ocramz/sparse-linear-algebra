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

import Data.Complex

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
v1 = V.fromList [0,3,4,6]

-- | intersection between sorted vectors, in-place updates
intersectWith ::
  Ord b => (a -> b) -> (a -> a -> c) -> V.Vector a -> V.Vector a -> V.Vector c
intersectWith f g u_ v_ = V.force $ V.create $ do
  let n = max (V.length u_) (V.length v_)
  vm <- VM.new n
  let go u_ v_ i vm | V.null u_ || V.null v_ || i == n = return (vm, i)
                    | otherwise =  do
         let (u,us) = (V.head u_, V.tail u_)
             (v,vs) = (V.head v_, V.tail v_)
         if f u == f v then do VM.write vm i (g u v)
                               go us vs (i + 1) vm
                   else if f u < f v then go us v_ i vm
                                     else go u_ vs i vm
  (vm', i') <- go u_ v_ 0 vm
  let vm'' = VM.take i' vm'
  return vm''

-- intersectIxVectors :: Ord a => V.Vector (a, c) -> V.Vector (a, c) -> V.Vector c
-- intersectIxVectors = intersectWith fst snd

liftIV ::
  Ord i => (a -> a -> b) -> V.Vector (i, a) -> V.Vector (i, a) -> V.Vector (i, b)
liftIV f = intersectWith fst (\(i, x) (_, y) -> (i, f x y))

dot :: (Ord i, Num b) => V.Vector (i, b) -> V.Vector (i, b) -> b
dot a b = foldlIxV' (+) 0 $ liftIV (*) a b 

dotC ::
  (RealFloat a, Ord t) => V.Vector (t, Complex a) -> V.Vector (t, Complex a) -> a
dotC a b = realPart $ foldlIxV' (+) 0 $ liftIV (\x y -> conjugate x * y) a b 

newtype IxVector a = IxVector { unIxVector :: V.Vector (Int, a)}

-- instance F.Foldable IxVector where
--   foldr f z v = foldr (\(_, x) a -> f a x) z (unIxVector v)

foldrIxV :: (t1 -> t2 -> t1) -> t1 -> V.Vector (t, t2) -> t1
foldrIxV f z v = V.foldr (\(_, x) a -> f a x) z v

foldlIxV' :: (t1 -> t2 -> t1) -> t1 -> V.Vector (t, t2) -> t1
foldlIxV' f z v = V.foldl' (\b (_, x) -> f b x) z v

e1, e2 :: V.Vector (Int, Double)
e1 = V.indexed $ V.fromList [1,0,0]
e2 = V.indexed $ V.fromList [0,1,0]

e1c = V.indexed $ V.fromList [1,0,0] :: V.Vector (Int, Complex Double)


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
