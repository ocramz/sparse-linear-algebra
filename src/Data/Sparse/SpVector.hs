{-# LANGUAGE FlexibleContexts #-}
{-# language TypeFamilies, MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module Data.Sparse.SpVector where

import Data.Sparse.Utils
import Data.Sparse.Types
import Data.Sparse.Internal.IntMap2

import Numeric.Eps
import Numeric.LinearAlgebra.Class

import Data.Complex
import Data.Maybe

import qualified Data.IntMap as IM
import qualified Data.Foldable as F
import qualified Data.Vector as V

-- * Sparse Vector

data SpVector a = SV { svDim :: {-# UNPACK #-} !Int ,
                       svData :: {-# UNPACK #-} !(IM.IntMap a)} deriving Eq

-- | SpVector sparsity
spySV :: Fractional b => SpVector a -> b
spySV s = fromIntegral (IM.size (dat s)) / fromIntegral (dim s)

-- | Number of nonzeros
nzSV :: SpVector a -> Int
nzSV sv = IM.size (dat sv)


sizeStrSV :: SpVector a -> String
sizeStrSV sv = unwords ["(",show (dim sv),"elements ) , ",show (nzSV sv),"NZ ( sparsity", show (spy sv),")"]



instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

instance Set SpVector where  
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
  
instance Foldable SpVector where
    foldr f d v = F.foldr f d (svData v)

instance Num a => AdditiveGroup (SpVector a) where
  zero = SV 0 IM.empty
  (^+^) = liftU2 (+)


-- | 'SpVector's form a vector space because they can be multiplied by a scalar
instance (Num e , AdditiveGroup e) => VectorSpace (SpVector e) where
  type (Scalar (SpVector e)) = e
  n .* v = scale n v

-- | 'SpVector's are finite-dimensional vectors
instance FiniteDim SpVector where
  type FDSize SpVector = Int
  dim = svDim  

instance HasData SpVector a where
  type HDData SpVector a = IM.IntMap a
  dat = svData

instance Sparse SpVector a where
  spy = spySV


-- | 'SpVector's are sparse containers too, i.e. any specific component may be missing (so it is assumed to be 0)
instance Elt a => SpContainer SpVector a where
  type ScIx SpVector = Int
  scInsert = insertSpVector
  scLookup v i = lookupSV i v
  v @@ i = lookupDenseSV i v


instance SparseVector SpVector Double where
  type SpvIx SpVector = Int
  svFromList = fromListSV
  svFromListDense = fromListDenseSV
  svConcat = foldr concatSV zero

-- instance SparseVector SpVector (Complex Double) where


-- | 'SpVector's form a Hilbert space, in that we can define an inner product over them
instance (AdditiveGroup e, Elt e) => Hilbert (SpVector e) where
  a `dot` b | dim a == dim b = dot (dat a) (dat b)
            | otherwise =
                     error $ "dot : sizes must coincide, instead we got " ++
                           show (dim a, dim b)
                           


-- | Since 'SpVector's form a Hilbert space, we can define a norm for them 
instance Normed (SpVector e) where
  norm p (SV _ v) = norm p v








-- ** Creation

-- | Empty sparse vector (length n, no entries)
zeroSV :: Int -> SpVector a
zeroSV n = SV n IM.empty


-- | Singleton sparse vector (length 1)
singletonSV :: a -> SpVector a
singletonSV x = SV 1 (IM.singleton 0 x)


-- | Canonical basis vector in R^n
ei :: Num a => Int -> IM.Key -> SpVector a
ei n i = SV n (IM.insert (i - 1) 1 IM.empty)



-- | create a sparse vector from an association list while discarding all zero entries
mkSpVector :: Epsilon a => Int -> IM.IntMap a -> SpVector a
mkSpVector d im = SV d $ IM.filterWithKey (\k v -> isNz v && inBounds0 d k) im

-- | ", from logically dense array (consecutive indices)
mkSpVectorD :: Epsilon a => Int -> [a] -> SpVector a
mkSpVectorD d ll = mkSpVector d (IM.fromList $ denseIxArray (take d ll))

-- ", don't filter zero elements
mkSpVector1 :: Int -> IM.IntMap a -> SpVector a
mkSpVector1 d ll = SV d $ IM.filterWithKey (\ k _ -> inBounds0 d k) ll

-- | Create new sparse vector, assumin 0-based, contiguous indexing
fromListDenseSV :: Int -> [a] -> SpVector a
fromListDenseSV d ll = SV d (IM.fromList $ denseIxArray (take d ll))


-- | Map a function over a range of indices and filter the result (indices and values) to fit in a `n`-long SpVector
spVectorDenseIx :: Epsilon a => (Int -> a) -> UB -> [Int] -> SpVector a
spVectorDenseIx f n ix =
  fromListSV n $ filter q $ zip ix $ map f ix where
    q (i, v) = inBounds0 n i && isNz v

-- | ", using just the integer bounds of the interval
spVectorDenseLoHi :: Epsilon a => (Int -> a) -> UB -> Int -> Int -> SpVector a
spVectorDenseLoHi f n lo hi = spVectorDenseIx f n [lo .. hi]    




  


-- | one-hot encoding : `oneHotSV n k` produces a SpVector of length n having 1 at the k-th position
oneHotSVU :: Num a => Int -> IxRow -> SpVector a
oneHotSVU n k = SV n (IM.singleton k 1)

oneHotSV :: Num a => Int -> IxRow -> SpVector a
oneHotSV n k |inBounds0 n k = oneHotSVU n k
             |otherwise = error "`oneHotSV n k` must satisfy 0 <= k <= n"


-- | DENSE vector of `1`s
onesSV :: Num a => Int -> SpVector a
onesSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 1

-- | DENSE vector of `0`s
zerosSV :: Num a => Int -> SpVector a
zerosSV d = SV d $ IM.fromList $ denseIxArray $ replicate d 0





-- *** Vector-related

-- | Populate a SpVector with the contents of a Vector. 
fromVector :: V.Vector a -> SpVector a
fromVector qv = V.ifoldl' ins (zeroSV n) qv where
  n = V.length qv
  ins vv i x = insertSpVector i x vv

-- | Populate a Vector with the entries of a SpVector, discarding the indices (NB: loses sparsity information).
toVector :: SpVector a -> V.Vector a
toVector = V.fromList . snd . unzip . toListSV

-- | -- | Populate a Vector with the entries of a SpVector, replacing the missing entries with 0
toVectorDense :: Num a => SpVector a -> V.Vector a
toVectorDense = V.fromList . toDenseListSV





-- ** Element insertion

-- |insert element `x` at index `i` in a preexisting SpVector
insertSpVector :: Int -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim)
  | inBounds0 d i = SV d (IM.insert i x xim)
  | otherwise = error "insertSpVector : index out of bounds"



-- ** fromList
fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (IM.fromList (filter (inBounds0 d . fst) iix ))

-- fromListSV' d iix = SV d (F.foldl' insf IM.empty iix') where
--   insf im (i, x) = IM.insert i x im
--   iix' = filter (inBounds0 d . fst) iix   -- filtering forces list as only instance
 
-- ** toList
toListSV :: SpVector a -> [(IM.Key, a)]
toListSV sv = IM.toList (dat sv)

-- |To dense list (default = 0)
toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d im) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]




-- | Indexed fold over SpVector
ifoldSV :: (IM.Key -> a -> b -> b) -> b -> SpVector a -> b
ifoldSV f e (SV d im) = IM.foldWithKey f e im







  
instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (IM.toList x)


-- ** Lookup

-- | Lookup an index in a SpVector
lookupSV :: IM.Key -> SpVector a -> Maybe a
lookupSV i (SV _ im) = IM.lookup i im

-- | Lookup an index, return a default value if lookup fails
lookupDefaultSV :: a -> IM.Key -> SpVector a -> a
lookupDefaultSV def i (SV _ im) = IM.findWithDefault def i im

-- |Lookup an index in a SpVector, returns 0 if lookup fails
lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV = lookupDefaultSV 0






-- ** Sub-vectors
-- | Tail elements
tailSV :: SpVector a -> SpVector a
tailSV (SV n sv) = SV (n-1) ta where
  ta = IM.mapKeys (\i -> i - 1) $ IM.delete 0 sv
  
-- | Head element
headSV :: Num a => SpVector a -> a
headSV sv = fromMaybe 0 (IM.lookup 0 (dat sv))

-- | Keep the first n components of the SpVector (like `take` for lists)
takeSV, dropSV :: Int -> SpVector a -> SpVector a
takeSV n (SV _ sv) = SV n $ IM.filterWithKey (\i _ -> i < n) sv
-- | Discard the first n components of the SpVector and rebalance the keys (like `drop` for lists)
dropSV n (SV n0 sv) = SV (n0 - n) $ IM.mapKeys (subtract n) $ IM.filterWithKey (\i _ -> i >= n) sv




-- | Concatenate two sparse vectors
concatSV :: SpVector a -> SpVector a -> SpVector a
concatSV (SV n1 s1) (SV n2 s2) = SV (n1+n2) (IM.union s1 s2') where
  s2' = IM.mapKeys (+ n1) s2


-- | Filter
filterSV :: (a -> Bool) -> SpVector a -> SpVector a
filterSV q sv = SV (dim sv) (IM.filter q (dat sv)) 


-- | Indexed filter
ifilterSV :: (Int -> a -> Bool) -> SpVector a -> SpVector a
ifilterSV q sv = SV (dim sv) (IM.filterWithKey q (dat sv))












-- * Orthogonal vector

-- | Generate an arbitrary (not random) vector `u` such that `v dot u = 0`
-- orthogonalSV :: Fractional a => SpVector a -> SpVector a
orthogonalSV v = u where
  (h, t) = (headSV v, tailSV v)
  n = dim v
  v2 = onesSV (n - 1)
  yn = singletonSV $ - (v2 `dot` t)/h
  u = concatSV yn v2



    











