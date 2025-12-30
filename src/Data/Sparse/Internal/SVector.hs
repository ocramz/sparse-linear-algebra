{-# language TypeFamilies, FlexibleInstances, MultiParamTypeClasses, CPP #-}
module Data.Sparse.Internal.SVector where

import Control.Arrow
import Control.Monad (unless)
import qualified Data.Foldable as F -- (foldl')
-- import Data.List (group, groupBy)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
-- import qualified Data.Vector.Algorithms.Merge as VA (sortBy)
-- import qualified Data.Vector.Generic as VG (convert)

import Control.Monad.Primitive

-- import Data.Sparse.Utils
-- import Data.Sparse.Types
import Data.Sparse.Internal.SVector.Mutable hiding (fromList)


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

instance HasData (SVector a) where
  type HDData (SVector a) = V.Vector a
  nnz = length . svIx
  dat (SV _ _ x) = x



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

-- | O(N) : Lookup an index in a SVector (based on `find` from Data.Foldable)
index :: SVector a -> Int -> Maybe a
index cv i =
  case F.find (== i) (svIx cv) of
    Just i' -> Just $ svVal cv V.! i'
    Nothing -> Nothing
      


-- * Lifted binary functions

-- | O(N) : Applies a function to the index _intersection_ of two CsrVector s. Useful e.g. to compute the inner product of two sparse vectors.
intersectWith :: (a -> b -> c) -> SVector a -> SVector b -> SVector c
intersectWith g v1 v2 = SV nfin ixf vf where
   nfin = V.length vf
   (ixf, vf) = V.unzip $ V.create (intersectWithM g v1 v2)

intersectWithM :: PrimMonad m =>
     (a -> b -> c)
     -> SVector a
     -> SVector b
     -> m (VM.MVector (PrimState m) (Int, c))
intersectWithM g (SV n1 ixu u) (SV n2 ixv v) = do
  vm <- VM.new n
  (vm', nfin) <- go ixu u ixv v 0 vm
  return $ VM.take nfin vm'
  where
    n = min n1 n2
    go ixu0 u0 ixv0 v0 i vm0
      | V.null u0 || V.null v0 || i == n = return (vm0, i)
      | otherwise =  do
           let (u', us) = headTail u0 
               (v', vs) = headTail v0 
               (ix1, ix1s) = headTail ixu0
               (ix2, ix2s) = headTail ixv0
           if ix1 == ix2 then do VM.write vm0 i (ix1, g u' v')
                                 go ix1s us ix2s vs (i + 1) vm0
                     else if ix1 < ix2 then go ix1s us ixv0 v0 i vm0
                                       else go ixu0 u0 ix2s vs i vm0

      
-- | O(N) : Applies a function to the index _union_ of two CsrVector s. Useful e.g. to compute the vector sum of two sparse vectors.
unionWith :: (t -> t -> a) -> t -> SVector t -> SVector t -> SVector a
unionWith g z v1 v2 = SV n ixf vf where
   n = max (svDim v1) (svDim v2)
   (ixf, vf) = V.unzip $ V.create (unionWithM g z v1 v2)

unionWithM :: PrimMonad m =>
     (a -> a -> b)
     -> a
     -> SVector a
     -> SVector a
     -> m (VM.MVector (PrimState m) (Int, b))
unionWithM g z (SV n1 ixu u) (SV n2 ixv v) = do
  vm <- VM.new n
  (vm', nfin) <- go ixu u ixv v 0 vm
  let vm'' = VM.take nfin vm'
  return vm''
  where
    n = min n1 n2
    go iu0 u0 iv0 v0 i vm0
          | (V.null u0 && V.null v0) || i == n = return (vm0 , i)
          | V.null u0 = do
                 VM.write vm0 i (V.head iv0, g z (V.head v0))
                 go iu0 u0 (V.tail iv0) (V.tail v0) (i + 1) vm0
          | V.null v0 = do
                 VM.write vm0 i (V.head iu0, g (V.head u0) z)
                 go (V.tail iu0) (V.tail u0) iv0 v0 (i + 1) vm0
          | otherwise =  do
           let (u', us) = headTail u0
               (v', vs) = headTail v0
               (iu1, ius) = headTail iu0
               (iv1, ivs) = headTail iv0
           if iu1 == iv1 then do VM.write vm0 i (iu1, g u' v')
                                 go ius us ivs vs (i + 1) vm0
                         else if iu1 < iv1 then do VM.write vm0 i (iu1, g u' z)
                                                   go ius us iv0 v0 (i + 1) vm0
                                           else do VM.write vm0 i (iv1, g z v')
                                                   go iu0 u0 ivs vs (i + 1) vm0
  


#define SVType(t) \
  instance AdditiveGroup (SVector t) where {zeroV = SV 0 V.empty V.empty; (^+^) = unionWith (+) 0 ; negateV = fmap negate};\
  instance VectorSpace (SVector t) where {type Scalar (SVector t) = (t); n .* x = fmap (* n) x };\
  instance InnerSpace (SVector t) where {a <.> b = sum $ intersectWith (*) a b}
  -- instance Normed (CsrVector t) where {norm p v = norm' p v}
  -- instance Hilbert (CsrVector t) where {x `dot` y = sum $ intersectWithCV (*) x y };\
  

#define SVTypeC(t) \
  instance AdditiveGroup (SVector (Complex t)) where {zeroV = SV 0 V.empty V.empty; (^+^) = unionWith (+) (0 :+ 0) ; negateV = fmap negate};\
  instance VectorSpace (SVector (Complex t)) where {type Scalar (SVector (Complex t)) = (Complex t); n .* x = fmap (* n) x };\
  instance InnerSpace (SVector (Complex t)) where {x <.> y = sum $ intersectWith (*) (conjugate <$> x) y};\

-- #define NormedType(t) \
--   instance Normed (CsrVector t) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v};\
--   instance Normed (CsrVector (Complex t)) where { norm p v | p==1 = norm1 v | otherwise = norm2 v ; normalize p v = v ./ norm p v}


-- -- CVType(Int)
-- -- CVType(Integer)
-- SVType(Float)
-- SVType(Double)
-- -- CVType(CSChar)
-- -- CVType(CInt)
-- -- CVType(CShort)
-- -- CVType(CLong)
-- -- CVType(CLLong)
-- -- CVType(CIntMax)
-- SVType(CFloat)
-- SVType(CDouble)

-- SVTypeC(Float)
-- SVTypeC(Double)  
-- SVTypeC(CFloat)
-- SVTypeC(CDouble)

-- -- NormedType(Float)
-- -- NormedType(Double)
-- -- NormedType(CFloat)
-- -- NormedType(CDouble)    








{-| Modify the mutable vector operand by applying a binary function over the index union of the two.

e.g.

g = (+)
z = 0
u = [(1, a), (2, b)]
v = [(0, d), (2, e)]

unionWithSMV g z u v = [(0, d), (1, a), (2, b + e)]

invariants :

* uu, vv nonzero values
* ixu, ixv nonzero indices
* ixu, ixv sorted in ascending order
* n1 == n2
* length ixu == length uu
* length ixv == length vv

-}
unionWithSMV :: PrimMonad m =>
     (a -> a -> a)
     -> a -> SVector a -> SMVector m a -> m (SMVector m a)
unionWithSMV g z (SV n ixu uu) (SMV n2 ixm_ vm_) = do
  ixmnew <- VM.new nnzero  -- create new mutable vectors
  vmnew <- VM.new nnzero
  unless (n == n2) (error "unionWithSMV : operand vectors must have the same length")
  (ixm, vm, nfin) <- go 0 ixmnew vmnew
  let ixm' = VM.take nfin ixm
      vm' = VM.take nfin vm
  return $ SMV n ixm' vm'
  where
    nnzero = lu + lv
    lu = V.length ixu
    lv = VM.length ixm_
    go i ixm vm
      | i == lu && i == lv || i == n = return (ixm, vm, i)
      | i == lu = do
          v0 <- VM.read vm i
          VM.write ixm i i
          VM.write vm i (g z v0)
          go (i + 1) ixm vm
      | i == lv = do
          let u0 = uu V.! i
          VM.write ixm_ i i
          VM.write vm_ i (g u0 z)
          go (i + 1) ixm vm
      | otherwise = do
          let u = uu V.! i    -- read head elements and indices
              iu = ixu V.! i
          v <- VM.read vm i
          iv <- VM.read ixm i
          if iu == iv then
            do
              VM.write ixm i iu         -- write `iu` at position `i` 
              VM.write vm i (g u v)
              go (i + 1) ixm vm
          else if iu < iv then
            do
              VM.write ixm i iu         -- write `iu` at position `i` 
              VM.write vm i (g u z)
              go (i + 1) ixm vm
              else
                do
                  VM.write ixm i iv     -- write `iv` at position `i` 
                  VM.write vm i (g z v)
                  go (i + 1) ixm vm
                               
                             

-- -- test data
-- testUnionWithSMV :: IO (SVector Double)
-- testUnionWithSMV = do 
--   let v = fromList 4 [(1, 1), (2, 1)]
--   vm <- SMV.fromList 4 [(0, pi), (2, pi)]
--   (vmres, nfin) <- unionWithSMV (+) 0 v vm
--   liftIO $ print nfin
--   freeze vmres




            
          

-- * To/from SMVector

thaw :: PrimMonad m => SVector a -> m (SMVector m a)
thaw (SV n ix v) = do
  vm <- V.thaw v
  ixm <- V.thaw ix
  return $ SMV n ixm vm

freeze :: PrimMonad m => SMVector m a -> m (SVector a)
freeze (SMV n ixm vm) = do
  v <- V.freeze vm
  ix <- V.freeze ixm
  return $ SV n ix v





-- * helpers

both :: Arrow arr => arr b c -> arr (b, b) (c, c)
both f = f *** f

headTail :: V.Vector a -> (a, V.Vector a)
headTail = V.head &&& V.tail
