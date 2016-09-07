module Math.Linear.Sparse where

-- import Data.Functor ((<$>))
import qualified Data.IntMap as IM
import qualified Data.Foldable as F

-- additive ring 
class Functor f => Additive f where
  -- | zero element
  zero :: Num a => f a
  
  -- | componentwise operations
  (^+^) :: Num a => f a -> f a -> f a
  (^-^) :: Num a => f a -> f a -> f a

  -- | linear interpolation. Doesn't really belong here though
  
  -- lerp :: Num a => a -> f a -> f a -> f a

  -- union binary lift
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- intersection binary lift
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c


negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate


class Additive f => VectorSpace f where
  (.*) :: Num a => a -> f a -> f a

class VectorSpace f => Normed f where
  dot :: Num a => f a -> f a -> a
  -- norm :: Num a => a -> f a -> a


instance Additive IM.IntMap where
  zero = IM.empty
  {-# INLINE zero #-}
  liftU2 = IM.unionWith
  {-# INLINE liftU2 #-}
  liftI2 = IM.intersectionWith
  {-# INLINE liftI2 #-}
  (^+^) = liftU2 (+)
  {-# INLINE (^+^) #-}
  x ^-^ y = x ^+^ negated y
  {-# INLINE (^-^) #-}

instance VectorSpace IM.IntMap where
  n .* im = IM.map (* n) im
  
-- instance Normed IM.IntMap where
--   im1 `dot` im2 = IM.foldr (*) im1 im2
 

-- | fold functions are applied to non-zero values
instance Foldable SparseVector where
    foldr f d v = F.foldr f d (svData v)
  


data SparseVector a = SV { svDim :: Int ,
                           svData :: IM.IntMap a} deriving Eq

instance Functor SparseVector where
  fmap f (SV n x) = SV n (fmap f x)

instance Additive SparseVector where
  zero = SV 0 IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
                      
instance VectorSpace SparseVector where
  n .* v = fmap (*n) v

-- instance Normed SparseVector where
--   v1 `dot` v2 = sum (IM.foldr (*) (svData v1) (svData v2))
  
