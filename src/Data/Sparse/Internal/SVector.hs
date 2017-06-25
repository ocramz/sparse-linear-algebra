{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses, CPP #-}
module Data.Sparse.Internal.SVector where

import Control.Arrow
import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
-- import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Data.Complex
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)
import Control.Monad.Primitive

-- import Data.Sparse.Utils
-- import Data.Sparse.Types

import Data.VectorSpace

import Numeric.LinearAlgebra.Class





data SVector a = SV { svDim :: {-# UNPACK #-} !Int,
                      svIx :: V.Vector Int,
                      svVal :: V.Vector a } deriving Eq

instance Show a => Show (SVector a) where
  show (SV n ix v) = unwords ["SV (",show n,"),",show nz,"NZ:",show (V.zip ix v)]
    where nz = V.length ix

instance Functor SVector where
  fmap f (SV n ix v) = SV n ix (fmap f v)

instance Foldable SVector where
  foldr f z (SV _ _ v) = foldr f z v

instance Traversable SVector where
  traverse f (SV n ix v) = SV n ix <$> traverse f v




-- ** Construction 

fromDense :: V.Vector a -> SVector a
fromDense xs = SV n (V.enumFromTo 0 (n-1)) xs where
  n = V.length xs

fromVector :: Int -> V.Vector (Int, a) -> SVector a
fromVector n ixs = SV n ix xs where
  (ix, xs) = V.unzip ixs

fromList :: Int -> [(Int, a)] -> SVector a
fromList n = fromVector n . V.fromList

-- * Query

-- | O(N) : Lookup an index in a CsrVector (based on `find` from Data.Foldable)
index :: SVector a -> Int -> Maybe a
index cv i =
  case F.find (== i) (svIx cv) of
    Just i' -> Just $ (V.!) (svVal cv) i'
    Nothing -> Nothing
      


-- * Lifted binary functions

-- | O(N) : Applies a function to the index _intersection_ of two CsrVector s. Useful e.g. to compute the inner product of two sparse vectors.
intersectWith :: (a -> b -> c) -> SVector a -> SVector b -> SVector c
intersectWith g (SV n1 ixu u) (SV n2 ixv v) = SV nfin ixf vf where
   nfin = V.length vf
   n = min n1 n2
   (ixf, vf) = V.unzip $ V.create (intersectWithM g ixu u ixv v)


intersectWithM :: (PrimMonad m, Ord a) =>
     (t1 -> t2 -> t)
     -> V.Vector a
     -> V.Vector t1
     -> V.Vector a
     -> V.Vector t2
     -> m (VM.MVector (PrimState m) (a, t))
intersectWithM g ixu u ixv v = do
  vm <- VM.new n
  (vm', nfin) <- go ixu u ixv v 0 vm
  return $ VM.take nfin vm'
  where
    (n1, n2) = both length (ixu, ixv)
    n = min n1 n2
    go ixu_ u_ ixv_ v_ i vm
      | V.null u_ || V.null v_ || i == n = return (vm, i)
      | otherwise =  do
           let (u, us) = headTail u_ 
               (v, vs) = headTail v_ 
               (ix1, ix1s) = headTail ixu_
               (ix2, ix2s) = headTail ixv_
           if ix1 == ix2 then do VM.write vm i (ix1, g u v)
                                 go ix1s us ix2s vs (i + 1) vm
                     else if ix1 < ix2 then go ix1s us ixv  v_ i vm
                                       else go ixu  u_ ix2s vs i vm





      
-- | O(N) : Applies a function to the index _union_ of two CsrVector s. Useful e.g. to compute the vector sum of two sparse vectors.
unionWith :: (t -> t -> a) -> t -> SVector t -> SVector t -> SVector a
unionWith g z (SV n1 ixu u) (SV n2 ixv v) = SV n ixf vf where
   n = max n1 n2
   (ixf, vf) = V.unzip $ V.create (unionWithM g z ixu u ixv v)

unionWithM :: (PrimMonad m, Ord t) =>
     (a -> a -> t1)
     -> a
     -> V.Vector t
     -> V.Vector a
     -> V.Vector t
     -> V.Vector a
     -> m (VM.MVector (PrimState m) (t, t1))
unionWithM g z ixu u ixv v = do
  vm <- VM.new n
  (vm', nfin) <- go ixu u ixv v 0 vm
  let vm'' = VM.take nfin vm'
  return vm''
  where
    (n1, n2) = both length (u, v)
    n = min n1 n2
    go iu u_ iv v_ i vm
          | (V.null u_ && V.null v_) || i == n = return (vm , i)
          | V.null u_ = do
                 VM.write vm i (V.head iv, g z (V.head v_))
                 go iu u_ (V.tail iv) (V.tail v_) (i + 1) vm
          | V.null v_ = do
                 VM.write vm i (V.head iu, g (V.head u_) z)
                 go (V.tail iu) (V.tail u_) iv v_ (i + 1) vm
          | otherwise =  do
           let (u, us) = headTail u_
               (v, vs) = headTail v_
               (iu1, ius) = headTail iu
               (iv1, ivs) = headTail iv
           if iu1 == iv1 then do VM.write vm i (iu1, g u v)
                                 go ius us ivs vs (i + 1) vm
                         else if iu1 < iv1 then do VM.write vm i (iu1, g u z)
                                                   go ius us iv v_ (i + 1) vm
                                           else do VM.write vm i (iv1, g z v)
                                                   go iu u_ ivs vs (i + 1) vm
  


#define SVType(t) \
  instance AdditiveGroup (SVector t) where {zeroV = SV 0 V.empty V.empty; (^+^) = unionWith (+) 0 ; negateV = fmap negate};\
  instance VectorSpace (SVector t) where {type Scalar (SVector t) = (t); n *^ x = fmap (* n) x };\
  instance InnerSpace (SVector t) where {a <.> b = sum $ intersectWith (*) a b}
  -- instance Normed (CsrVector t) where {norm p v = norm' p v}
  -- instance Hilbert (CsrVector t) where {x `dot` y = sum $ intersectWithCV (*) x y };\
  

#define SVTypeC(t) \
  instance AdditiveGroup (SVector (Complex t)) where {zeroV = SV 0 V.empty V.empty; (^+^) = unionWith (+) (0 :+ 0) ; negateV = fmap negate};\
  instance VectorSpace (SVector (Complex t)) where {type Scalar (SVector (Complex t)) = (Complex t); n *^ x = fmap (* n) x };\
  instance InnerSpace (SVector (Complex t)) where {x <.> y = sum $ intersectWith (*) (conjugate <$> x) y};\

-- #define NormedType(t) \
--   instance Normed (CsrVector t) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v};\
--   instance Normed (CsrVector (Complex t)) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v}


-- CVType(Int)
-- CVType(Integer)
SVType(Float)
SVType(Double)
-- CVType(CSChar)
-- CVType(CInt)
-- CVType(CShort)
-- CVType(CLong)
-- CVType(CLLong)
-- CVType(CIntMax)
SVType(CFloat)
SVType(CDouble)

SVTypeC(Float)
SVTypeC(Double)  
SVTypeC(CFloat)
SVTypeC(CDouble)

-- NormedType(Float)
-- NormedType(Double)
-- NormedType(CFloat)
-- NormedType(CDouble)    




-- * helpers

both :: Arrow arr => arr b c -> arr (b, b) (c, c)
both f = f *** f

headTail :: V.Vector a -> (a, V.Vector a)
headTail = V.head &&& V.tail
