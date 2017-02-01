{-# LANGUAGE FlexibleContexts #-}
{-# language TypeFamilies, MultiParamTypeClasses, FlexibleInstances #-}
{-# language CPP #-}
{-# language GeneralizedNewtypeDeriving, DeriveFunctor #-}
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

import Control.Exception
import Control.Monad.Catch (MonadThrow (..))
import Control.Exception.Common

import GHC.Exts

import Data.Sparse.Utils
import Data.Sparse.Types
import Data.Sparse.Internal.IntMap2

import Numeric.Eps
import Numeric.LinearAlgebra.Class

import Data.Complex
import Data.Maybe

import qualified Data.IntMap.Strict as IM
import qualified Data.Foldable as F
import qualified Data.Vector as V

import Data.VectorSpace hiding (magnitude)


-- | A synonym for IntMap 
newtype IntM a = IntM {unIM :: IM.IntMap a} deriving (Eq, Show, Functor, Foldable)

empty = IntM IM.empty

size (IntM x) = IM.size x

singleton i x = IntM $ IM.singleton i x

filterWithKey f im = IntM $ IM.filterWithKey f (unIM im)

insert k x (IntM im) = IntM $ IM.insert k x im

filterI f (IntM im) = IntM $ IM.filter f im

instance IsList (IntM a) where
  type Item (IntM a) = (Int, a)
  fromList = IntM . IM.fromList
  toList = IM.toList . unIM

instance Set IntM where
  liftU2 f (IntM a) (IntM b) = IntM $ IM.unionWith f a b
  liftI2 f (IntM a) (IntM b) = IntM $ IM.intersectionWith f a b

instance Num a => AdditiveGroup (IntM a) where
  zeroV = IntM IM.empty
  {-# INLINE zeroV #-}
  (^+^) = liftU2 (+)
  {-# INLINE (^+^) #-}
  (^-^) = liftU2 (-)
  {-# INLINE (^-^) #-}
  negateV = fmap negate
  {-# INLINE negateV #-}




-- | ParamInstance can be used with all types that are instances of Set (which are by construction also instances of Functor)
#define ParamInstance(f, t) \
  instance VectorSpace (f t) where {type (Scalar (f (t))) = (t); n *^ im = fmap (* n) im};\
  instance VectorSpace (f (Complex t)) where {type (Scalar (f (Complex t))) = Complex (t); n *^ im = fmap (* n) im};\
  instance InnerSpace (f t) where {a <.> b = sum $ liftI2 (*) a b};\
  instance InnerSpace (f (Complex t)) where {a <.> b = sum $ liftI2 (*) (conjugate <$> a) b};\
  -- instance Normed (f t) where {type RealScalar (f t) = t ; type Magnitude (f t) = t ; norm1 a = sum (abs <$> a) ; norm2Sq a = sum $ liftI2 (*) a a; normP p v = sum u**(1/p) where u = fmap (**p) v; normalize = normzPR ; normalize2 = normz2R}; \
  -- instance Normed (f (Complex t)) where {type RealScalar (f (Complex t)) = t; type Magnitude (f (Complex t)) = t; norm1 a = realPart $ sum (abs <$> a); norm2Sq a = realPart $ sum $ liftI2 (*) (conjugate <$> a) a; normP p v = realPart $ sum u**(1/(p :+ 0)) where u = fmap (**(p :+ 0)) v; normalize = normzPC; normalize2 = normz2C }


instance Normed (IntM Double) where
  type RealScalar (IntM Double) = Double
  type Magnitude (IntM Double) = Double
  norm1 a = sum (abs <$> a)
  norm2Sq a = sum $ liftI2 (*) a a
  normP p v = sum u**(1/p) where u = fmap (**p) v
  normalize p v = v ./ normP p v 
  normalize2 v = v ./ norm2 v 
  
instance Normed (IntM (Complex Double)) where
  type RealScalar (IntM (Complex Double)) = Double
  type Magnitude (IntM (Complex Double)) = Double
  norm1 a = realPart $ sum (abs <$> a)
  norm2Sq a = realPart $ sum $ liftI2 (*) (conjugate <$> a) a
  normP p v = realPart $ sum u**(1/(p :+ 0)) where u = fmap (**(p :+ 0)) v
  normalize p v = v ./ toC (normP p v)
  normalize2 v = v ./ toC (norm2 v)






-- | IntMap instances
#define IntMapInstance(t) \
  ParamInstance( IntM, t )

IntMapInstance(Double)
-- IntMapInstance(Float)







-- -- | list to IntMap
-- mkIm :: [Double] -> IM.IntMap Double
mkIm xs = fromList $ indexed xs :: IntM Double

-- mkImC :: [Complex Double] -> IM.IntMap (Complex Double)
mkImC xs = fromList $ indexed xs :: IntM (Complex Double)





-- * Sparse Vector

data SpVector a = SV { svDim :: {-# UNPACK #-} !Int ,
                       svData :: !(IntM a)} deriving Eq

instance Show a => Show (SpVector a) where
  show (SV d x) = "SV (" ++ show d ++ ") "++ show (toList x)

-- | SpVector sparsity
spySV :: Fractional b => SpVector a -> b
spySV s = fromIntegral (size (dat s)) / fromIntegral (dim s)

-- | Number of nonzeros
nzSV :: SpVector a -> Int
nzSV sv = size (dat sv)


sizeStrSV :: SpVector a -> String
sizeStrSV sv = unwords ["(",show (dim sv),"elements ) , ",show (nzSV sv),"NZ ( sparsity", show (spy sv),")"]



instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

instance Set SpVector where  
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
  
instance Foldable SpVector where
    foldr f d v = F.foldr f d (svData v)











-- | 'SpVector's form a vector space because they can be multiplied by a scalar


-- | 'SpVector's are finite-dimensional vectors
instance FiniteDim SpVector where
  type FDSize SpVector = Int
  dim = svDim  

instance HasData SpVector a where
  type HDData SpVector a = IntM a
  dat = svData

instance Sparse SpVector a where
  spy = spySV


-- | 'SpVector's are sparse containers too, i.e. any specific component may be missing (so it is assumed to be 0)
instance Elt a => SpContainer SpVector a where
  type ScIx SpVector = Int
  scInsert = insertSpVector
  scLookup v i = lookupSV i v
  scToList = toListSV
  v @@ i = lookupDenseSV i v


-- instance (Elt e, RealFloat e) => SparseVector SpVector e where
--   type SpvIx SpVector = Int
--   svFromList = fromListSV
--   svFromListDense = fromListDenseSV
--   svConcat = foldr concatSV zero

-- instance SparseVector SpVector (Complex Double) where




#define SpVectorInstance(t) \
  instance AdditiveGroup (SpVector (t)) where { zeroV = SV 0 empty; (^+^) = liftU2 (+); negateV = fmap negate };\
  instance AdditiveGroup (SpVector (Complex t)) where { zeroV = SV 0 empty; (^+^) = liftU2 (+); negateV = fmap negate };\
  instance VectorSpace (SpVector t) where { type (Scalar (SpVector t)) = t; n *^ v = scale n v};\
  instance VectorSpace (SpVector (Complex t)) where { type (Scalar (SpVector (Complex t))) = Complex t; n *^ v = scale n v};\
  instance InnerSpace (SpVector (t)) where { (<.>) = dotS };\
  instance InnerSpace (SpVector (Complex (t))) where { (<.>) = dotS };\
  instance Normed (SpVector (t)) where {type RealScalar (SpVector (t)) = t; type Magnitude (SpVector (t)) = t; norm1 (SV _ v) = norm1 v; norm2Sq (SV _ v) = norm2Sq v ; normP p (SV _ v) = normP p v; normalize p (SV n v) = SV n (normalize p v); normalize2 (SV n v) = SV n (normalize2 v)};\
  instance Normed (SpVector (Complex t)) where {type RealScalar (SpVector (Complex t)) = t; type Magnitude (SpVector (Complex t)) = t; norm1 (SV _ v) = norm1 v; norm2Sq (SV _ v) = norm2Sq v ; normP p (SV _ v) = normP p v; normalize p (SV n v) = SV n (normalize p v); normalize2 (SV n v) = SV n (normalize2 v)}



SpVectorInstance(Double)
-- SpVectorInstance(Float)


dotS :: InnerSpace (IntM t) => SpVector t -> SpVector t -> Scalar (IntM t)
(SV m a) `dotS` (SV n b)
  | n == m = a <.> b
  | otherwise = error $ unwords ["<.> : Incompatible dimensions:", show m, show n]

-- dotSSafe :: (MonadThrow m, InnerSpace (IM.IntMap t)) =>
--      SpVector t -> SpVector t -> m (Scalar (IM.IntMap t))
dotSSafe (SV m a) (SV n b)
  | n == m = return $ a <.> b
  | otherwise = throwM (DotSizeMismatch m n)










-- ** Creation

-- | Empty sparse vector (length n, no entries)
zeroSV :: Int -> SpVector a
zeroSV n = SV n empty


-- | Singleton sparse vector (length 1)
singletonSV :: a -> SpVector a
singletonSV x = SV 1 (singleton 0 x)


-- | Canonical basis vector in R^n
-- ei :: Num a => Int -> IM.Key -> SpVector a
ei n i = SV n (insert (i - 1) 1 empty)



-- | Sparse vector from an association list while discarding all zero entries
-- mkSpVector :: Epsilon a => Int -> IntM a -> SpVector a
mkSpVector d im = SV d $ IntM $ IM.filterWithKey (\k v -> isNz v && inBounds0 d k) im


-- ", don't filter zero elements
-- mkSpVector1 :: Int -> IM.IntMap a -> SpVector a
mkSpVector1 d ll = SV d $ IntM $ IM.filterWithKey (\ k _ -> inBounds0 d k) ll


-- | Dense real SpVector (monomorphic Double)
mkSpVR :: Int -> [Double] -> SpVector Double
mkSpVR d ll = SV d $ mkIm ll

-- | Dense complex SpVector (monomorphic Double)
mkSpVC :: Int -> [Complex Double] -> SpVector (Complex Double)
mkSpVC d ll = SV d $ mkImC ll




-- | Create new sparse vector, assumin 0-based, contiguous indexing
fromListDenseSV :: Int -> [a] -> SpVector a
fromListDenseSV d ll = SV d (fromList $ denseIxArray (take d ll))


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
oneHotSVU n k = SV n (singleton k 1)

oneHotSV :: Num a => Int -> IxRow -> SpVector a
oneHotSV n k |inBounds0 n k = oneHotSVU n k
             |otherwise = error "`oneHotSV n k` must satisfy 0 <= k <= n"


-- | DENSE vector of `1`s
onesSV :: Num a => Int -> SpVector a
onesSV d = SV d $ fromList $ denseIxArray $ replicate d 1

-- | DENSE vector of `0`s
zerosSV :: Num a => Int -> SpVector a
zerosSV d = SV d $ fromList $ denseIxArray $ replicate d 0





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

-- |insert element `x` at index `i` in a preexisting SpVector; discards out-of-bounds entries
-- insertSpVector :: IM.Key -> a -> SpVector a -> SpVector a
insertSpVector i x (SV d xim) | inBounds0 d i = SV d (insert i x xim)

insertSpVectorSafe :: MonadThrow m => Int -> a -> SpVector a -> m (SpVector a)
insertSpVectorSafe i x (SV d xim)
  | inBounds0 d i = return $ SV d (insert i x xim)
  | otherwise = throwM (OOBIxError "insertSpVector" i)





-- ** fromList
fromListSV :: Int -> [(Int, a)] -> SpVector a
fromListSV d iix = SV d (fromList (filter (inBounds0 d . fst) iix ))

 
-- ** toList
toListSV :: SpVector a -> [(IM.Key, a)]
toListSV sv = toList (dat sv)

-- |To dense list (default = 0)
toDenseListSV :: Num b => SpVector b -> [b]
toDenseListSV (SV d (IntM im)) = fmap (\i -> IM.findWithDefault 0 i im) [0 .. d-1]




-- | Indexed fold over SpVector
ifoldSV :: (IM.Key -> a -> b -> b) -> b -> SpVector a -> b
ifoldSV f e (SV d (IntM im)) = IM.foldrWithKey f e im







  



-- ** Lookup

-- | Lookup an index in a SpVector
lookupSV :: IM.Key -> SpVector a -> Maybe a
lookupSV i (SV _ (IntM im)) = IM.lookup i im

-- | Lookup an index, return a default value if lookup fails
lookupDefaultSV :: a -> IM.Key -> SpVector a -> a
lookupDefaultSV def i (SV _ (IntM im)) = IM.findWithDefault def i im

-- |Lookup an index in a SpVector, returns 0 if lookup fails
lookupDenseSV :: Num a => IM.Key -> SpVector a -> a
lookupDenseSV = lookupDefaultSV 0






-- ** Sub-vectors
-- | Tail elements
tailSV :: SpVector a -> SpVector a
tailSV (SV n (IntM sv)) = SV (n-1) $ IntM ta where
  ta = IM.mapKeys (\i -> i - 1) $ IM.delete 0 sv
  
-- | Head element
headSV :: Num a => SpVector a -> a
headSV (SV _ (IntM im)) = fromMaybe 0 (IM.lookup 0 im)

-- | Keep the first n components of the SpVector (like `take` for lists)
takeSV, dropSV :: Int -> SpVector a -> SpVector a
takeSV n (SV _ sv) = SV n $ filterWithKey (\i _ -> i < n) sv
-- | Discard the first n components of the SpVector and rebalance the keys (like `drop` for lists)
dropSV n (SV n0 (IntM sv)) = SV (n0 - n) $ IntM $ IM.mapKeys (subtract n) $ IM.filterWithKey (\i _ -> i >= n) sv


-- | Keep a range of entries 
rangeSV :: (IM.Key, IM.Key) -> SpVector a -> SpVector a
rangeSV (rmin, rmax) (SV n (IntM sv))
  | len > 0 && len <= n = SV len $ IntM sv'
  | otherwise = error $ unwords ["rangeSV : invalid bounds", show (rmin, rmax) ] where
  len = rmax - rmin
  sv' = IM.mapKeys (subtract rmin) $ IM.filterWithKey (\i _ -> i >= rmin && i <= rmax) sv





-- | Concatenate two sparse vectors
concatSV :: SpVector a -> SpVector a -> SpVector a
concatSV (SV n1 (IntM s1)) (SV n2 (IntM s2)) = SV (n1+n2) $ IntM (IM.union s1 s2') where
  s2' = IM.mapKeys (+ n1) s2


-- | Filter
filterSV :: (a -> Bool) -> SpVector a -> SpVector a
filterSV q sv = SV (dim sv) $ IntM (IM.filter q (unIM $ dat sv)) 


-- | Indexed filter
ifilterSV :: (Int -> a -> Bool) -> SpVector a -> SpVector a
ifilterSV q sv = SV (dim sv) (filterWithKey q (dat sv))






-- * Sparsify : remove almost-0 elements (|x| < eps)
-- | Sparsify an SpVector
sparsifySV :: Epsilon a => SpVector a -> SpVector a
sparsifySV = filterSV isNz








-- * Orthogonal vector

-- | Generate an arbitrary (not random) vector `u` such that `v dot u = 0`
-- orthogonalSV :: Fractional a => SpVector a -> SpVector a
orthogonalSV v = u where
  (h, t) = (headSV v, tailSV v)
  n = dim v
  v2 = onesSV (n - 1)
  yn = singletonSV $ - (v2 `dot` t)/h
  u = concatSV yn v2



    











