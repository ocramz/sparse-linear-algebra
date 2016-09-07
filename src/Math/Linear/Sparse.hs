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
  
instance Normed IM.IntMap where
   a `dot` b = sum $ liftI2 (*) a b 


   



data SpVector a = SV { svDim :: Int ,
                           svData :: IM.IntMap a} deriving Eq

                      
-- instances for SparseVector
instance Functor SpVector where
  fmap f (SV n x) = SV n (fmap f x)

-- | fold functions are applied to non-zero values
instance Foldable SpVector where
    foldr f d v = F.foldr f d (svData v)

instance Additive SpVector where
  zero = SV 0 IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  liftU2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftU2 f2 x1 x2)
  liftI2 f2 (SV n1 x1) (SV n2 x2) = SV (max n1 n2) (liftI2 f2 x1 x2)
                      
instance VectorSpace SpVector where
  n .* v = fmap (*n) v

instance Normed SpVector where
  sv1 `dot` sv2 = dot (svData sv1) (svData sv2)




-- Sparse Matrices
data SpMatrix a = SM {smDim :: (Int, Int),
                      smData :: IM.IntMap (SpVector a)} deriving Eq

instance Functor SpMatrix where
  fmap f (SM d md) = SM d ((fmap . fmap) f md)

instance Additive SpMatrix where
  zero = SM (0,0) IM.empty
  (^+^) = liftU2 (+)
  (^-^) = liftU2 (-)
  -- liftU2 f2 (SM n1 x1) (SM n2 x2) = SM (maxTup n1 n2) (liftU2 f2 x1 x2)
  -- liftI2 f2 (SM n1 x1) (SM n2 x2) = SM (minTup n1 n2) (liftI2 f2 x1 x2)


maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)



  
-- -- testing testing

v1, v2 :: SpVector Double
v1 = SV 5 (IM.fromList [(0, 4), (1, 3)])
v2 = SV 4 (IM.fromList [(0, 2), (1, 3.2), (3, 15)])
