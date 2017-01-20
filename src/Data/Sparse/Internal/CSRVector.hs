{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses, CPP #-}
module Data.Sparse.Internal.CSRVector where

import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
-- import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Data.Complex
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)

-- import Data.Sparse.Utils
-- import Data.Sparse.Types

import Data.VectorSpace

import Numeric.LinearAlgebra.Class





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

-- | O(N) : Lookup an index in a CsrVector (based on `find` from Data.Foldable)
indexCV :: CsrVector a -> Int -> Maybe a
indexCV cv i =
  case F.find (== i) (cvIx cv) of
    Just i' -> Just $ (V.!) (cvVal cv) i'
    Nothing -> Nothing
      


-- * Lifted binary functions

-- | O(N) : Applies a function to the index _intersection_ of two CsrVector s. Useful e.g. to compute the inner product of two sparse vectors.
intersectWithCV :: (a -> b -> c) -> CsrVector a -> CsrVector b -> CsrVector c
intersectWithCV g (CV n1 ixu_ uu_) (CV n2 ixv_ vv_) = CV nfin ixf vf where
   nfin = V.length vf
   n = min n1 n2
   (ixf, vf) = V.unzip $ V.create $ do
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


      
-- | O(N) : Applies a function to the index _union_ of two CsrVector s. Useful e.g. to compute the vector sum of two sparse vectors.
unionWithCV :: (t -> t -> a) -> t -> CsrVector t -> CsrVector t -> CsrVector a
unionWithCV g z (CV n1 ixu_ uu_) (CV n2 ixv_ vv_) = CV n ixf vf where
   n = max n1 n2
   (ixf, vf) = V.unzip $ V.create $ do
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



#define CVType(t) \
  instance AdditiveGroup (CsrVector t) where {zeroV = CV 0 V.empty V.empty; (^+^) = unionWithCV (+) 0 ; negateV = fmap negate};\
  instance VectorSpace (CsrVector t) where {type Scalar (CsrVector t) = (t); n *^ x = fmap (* n) x };\
  -- instance Hilbert (CsrVector t) where {x `dot` y = sum $ intersectWithCV (*) x y };\
  -- -- instance Normed (CsrVector t) where {norm p v = norm' p v}

#define CVTypeC(t) \
  instance AdditiveGroup (CsrVector (Complex t)) where {zeroV = CV 0 V.empty V.empty; (^+^) = unionWithCV (+) (0 :+ 0) ; negateV = fmap negate};\
  instance VectorSpace (CsrVector (Complex t)) where {type Scalar (CsrVector (Complex t)) = (Complex t); n *^ x = fmap (* n) x };\
  -- instance Hilbert (CsrVector (Complex t)) where {x `dot` y = sum $ intersectWithCV (*) (fmap conjugate x) y};\

-- #define NormedType(t) \
--   instance Normed (CsrVector t) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v};\
--   instance Normed (CsrVector (Complex t)) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v}


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

-- NormedType(Float)
-- NormedType(Double)
-- NormedType(CFloat)
-- NormedType(CDouble)    
