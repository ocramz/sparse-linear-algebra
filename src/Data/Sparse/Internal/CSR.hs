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

import Data.Maybe

import Data.Complex
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)


import Data.Sparse.Utils
import Data.Sparse.Types

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

instance Show a => Show (CsrMatrix a) where
  show m'@(CM m n nz cix rp x) = szs where
    szs = unwords ["CSR (",show m, "x", show n,"),",show nz, "NZ ( sparsity",show (spy m'),"), column indices:",show cix,", row pointers:", show rp,", data:",show x]

-- * Creation
-- | Copy a Vector of (row index, column index, entry) into a CSR structure. Sorts the Vector by row indices ( O(log N) ), unzips column indices and data ( O(N) ) and generates the row pointer vector ( 2 passes of O(N) )
toCSR :: Int -> Int -> V.Vector (Int, Int, a) -> CsrMatrix a
toCSR m n ijxv = CM m n nz cix crp x where
  nz = V.length x
  (rix, cix, x) = V.unzip3 (sortByRows ijxv)
  crp = csrPtrVM m rix
  sortByRows = V.modify (VA.sortBy f) where
       f a b = compare (fst3 a) (fst3 b)
  csrPtrVM nrows xs = V.scanl (+) 0 $ V.create createf where
   createf :: ST s (VM.MVector s Int)
   createf = do
     vm <- VM.new nrows
     let loop v ll i | i == nrows = return ()
                     | otherwise = do
                                     let lp = V.length $ V.takeWhile (== i) ll
                                     VM.write v i lp
                                     loop v (V.drop lp ll) (i + 1)
     loop vm xs 0
     return vm


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
  trimf i1 i2 w = V.drop i1 (V.take i2 w)


-- matVec mm v = fmap (\ri -> extractRow mm ri `dot` v) irows_ where
--   nrows = csrNrows mm
--   irows_ = [0 .. nrows - 1]


-- | Rebuilds the (row, column, entry) Vector from the CSR representation. Not optimized for efficiency.
fromCSR :: CsrMatrix a -> V.Vector (Int, Int, a)
fromCSR mc = mconcat $ map (\i -> withRowIx (extractRowCSR mc i) i) [0 .. csrNrows mc - 1] where
  withRowIx (CV n icol_ v_) i = V.zip3 (V.replicate n i) icol_ v_


-- NOT OPTIMIZED : 
transposeCSR m1@(CM m n _ _ _ _) = toCSR n m $ V.zip3 jj ii xx where
  (ii, jj, xx) = V.unzip3 $ fromCSR m1




-- * Pretty printing
-- toDenseRow :: Elt a => CsrMatrix a -> Int -> [a]
-- toDenseRow cm irow =
--   fmap (\icol -> cm @@ (irow,icol)) [0..ncol-1] where (_, ncol) = dim cm


-- * Some instances

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





-- * Sparse vector

data CsrVector a = CV { cvDim :: {-# UNPACK #-} !Int,
                        cvIx :: V.Vector Int,
                        cvVal :: V.Vector a } deriving Eq

instance Show a => Show (CsrVector a) where
  show (CV n ix v) = unwords ["CV (",show n,"),",show nz,"NZ:",show (V.zip ix v)]
    where nz = V.length ix

instance Functor CsrVector where
  fmap f (CV n ix v) = CV n ix (fmap f v)

instance Foldable CsrVector where
  foldr f z (CV _ _ v) = foldr f z v

instance Traversable CsrVector where
  traverse f (CV n ix v) = CV n ix <$> traverse f v




#define CVType(t) \
  instance AdditiveGroup (CsrVector t) where {zero = CV 0 V.empty V.empty; (^+^) = unionWithCV (+) 0 ; negated = fmap negate};\
  instance VectorSpace (CsrVector t) where {type Scalar (CsrVector t) = (t); n .* x = fmap (* n) x };\
  instance Hilbert (CsrVector t) where {x `dot` y = sum $ intersectWithCV (*) x y };\
  -- instance Normed (CsrVector t) where {norm p v = norm' p v}

#define CVTypeC(t) \
  instance AdditiveGroup (CsrVector (Complex t)) where {zero = CV 0 V.empty V.empty; (^+^) = unionWithCV (+) (0 :+ 0) ; negated = fmap negate};\
  instance VectorSpace (CsrVector (Complex t)) where {type Scalar (CsrVector (Complex t)) = (Complex t); n .* x = fmap (* n) x };\
  instance Hilbert (CsrVector (Complex t)) where {x `dot` y = sum $ intersectWithCV (*) (fmap conjugate x) y};\

#define NormedType(t) \
  instance Normed (CsrVector t) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v};\
  instance Normed (CsrVector (Complex t)) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v}


CVType(Int)
CVType(Integer)
CVType(Float)
CVType(Double)
CVType(CSChar)
CVType(CInt)
CVType(CShort)
CVType(CLong)
CVType(CLLong)
CVType(CIntMax)
CVType(CFloat)
CVType(CDouble)

CVTypeC(Float)
CVTypeC(Double)  
CVTypeC(CFloat)
CVTypeC(CDouble)

NormedType(Float)
NormedType(Double)
NormedType(CFloat)
NormedType(CDouble)

-- ** Construction 

fromDenseCV :: V.Vector a -> CsrVector a
fromDenseCV xs = CV n (V.enumFromTo 0 (n-1)) xs where
  n = V.length xs

fromVectorCV :: Int -> V.Vector (Int, a) -> CsrVector a
fromVectorCV n ixs = CV n ix xs where
  (ix, xs) = V.unzip ixs

fromListCV :: Int -> [(Int, a)] -> CsrVector a
fromListCV n = fromVectorCV n . V.fromList

-- * Query

-- | O(N) Lookup an index in a CsrVector (based on `find` from Data.Foldable)
indexCV :: CsrVector a -> Int -> Maybe a
indexCV cv i =
  case F.find (== i) (cvIx cv) of
    Just i' -> Just $ (V.!) (cvVal cv) i'
    Nothing -> Nothing
      




-- | Intersection between sorted vectors, in-place updates
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


-- | O(N) Applies a function to the index _intersection_ of two CsrVector s. Useful e.g. to compute the inner product of two sparse vectors.
intersectWithCV :: (a -> b -> c) -> CsrVector a -> CsrVector b -> CsrVector c
intersectWithCV g (CV n1 ixu_ uu_) (CV n2 ixv_ vv_) = CV nfin ixf vf where
   nfin = V.length vf
   (ixf, vf) = V.unzip $ V.force $ V.create $ do
    let n = max n1 n2
    vm <- VM.new n
    let go ixu u_ ixv v_ i vm | V.null u_ || V.null v_ || i == n = return (vm, i)
                              | otherwise =  do
           let (u,us) = (V.head u_, V.tail u_)
               (v,vs) = (V.head v_, V.tail v_)
               (ix1,ix1s) = (V.head ixu, V.tail ixu)
               (ix2,ix2s) = (V.head ixv, V.tail ixv)
           if ix1 == ix2 then do VM.write vm i (ix1, g u v)
                                 go ix1s us ix2s vs (i + 1) vm
                     else if ix1 < ix2 then go ix1s us ixv v_ i vm
                                       else go ixu u_ ix2s vs i vm
    (vm', nfin) <- go ixu_ uu_ ixv_ vv_ 0 vm
    let vm'' = VM.take nfin vm'
    return vm''

-- | O(N) Applies a function to the index _union_ of two CsrVector s. Useful e.g. to compute the vector sum of two sparse vectors.
unionWithCV :: (t -> t -> a) -> t -> CsrVector t -> CsrVector t -> CsrVector a
unionWithCV g z (CV n1 ixu_ uu_) (CV n2 ixv_ vv_) = CV n ixf vf where
   n = max n1 n2
   (ixf, vf) = V.unzip $ V.force $ V.create $ do
    vm <- VM.new n
    let go iu u_ iv v_ i vm
          | (V.null u_ && V.null v_) || i == n = return (vm , i)
          | V.null u_ = do
                 VM.write vm i (V.head iv, g z (V.head v_))
                 go iu u_ (V.tail iv) (V.tail v_) (i + 1) vm
          | V.null v_ = do
                 VM.write vm i (V.head iu, g (V.head u_) z)
                 go (V.tail iu) (V.tail u_) iv v_ (i + 1) vm
          | otherwise =  do
           let (u,us) = (V.head u_, V.tail u_)
               (v,vs) = (V.head v_, V.tail v_)
               (iu1,ius) = (V.head iu, V.tail iu)
               (iv1,ivs) = (V.head iv, V.tail iv)
           if iu1 == iv1 then do VM.write vm i (iu1, g u v)
                                 go ius us ivs vs (i + 1) vm
                         else if iu1 < iv1 then do VM.write vm i (iu1, g u z)
                                                   go ius us iv v_ (i + 1) vm
                                           else do VM.write vm i (iv1, g z v)
                                                   go iu u_ ivs vs (i + 1) vm
    (vm', nfin) <- go ixu_ uu_ ixv_ vv_ 0 vm
    let vm'' = VM.take nfin vm'
    return vm''                                                       
           




  







-- * Utilities

tail3 :: (t, t1, t2) -> (t1, t2)
tail3 (_,j,x) = (j,x)

snd3 (_,j,_) = j

fst3 :: (t, t1, t2) -> t
fst3 (i, _, _) = i





-- test data

l0 = [1,2,4,5,8]
l1 = [2,3,6]
l2 = [7]

v0,v1 :: V.Vector Int
v0 = V.fromList [0,1,2,5,6]
v1 = V.fromList [0,3,4,6]

-- e1, e2 :: V.Vector (Int, Double)
-- e1 = V.indexed $ V.fromList [1,0,0]
-- e2 = V.indexed $ V.fromList [0,1,0]

e1, e2:: CsrVector Double
e1 = fromListCV 4 [(0, 1)] 
e2 = fromListCV 4 [(1, 1)]
e3 = fromListCV 4 [(0, 1 :+ 2)] :: CsrVector (Complex Double)

e1c = V.indexed $ V.fromList [1,0,0] :: V.Vector (Int, Complex Double)

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



union :: Ord a => [a] -> [a] -> [a]
union u_ v_ = go u_ v_ where
  go [] x = x
  go y [] = y
  go uu@(u:us) vv@(v:vs)
    | u == v =    u : go us vs
    | u < v =     u : go us vv 
    | otherwise = v : go uu vs
