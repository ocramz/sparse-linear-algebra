{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, KindSignatures #-}
module Numeric.LinearAlgebra.Class where

-- * Additive ring 
class Functor f => Additive f where
  -- | Ring zero element
  zero :: Num a => f a
  
  -- | Ring +
  (^+^) :: Num a => f a -> f a -> f a


  one :: Num a => f a

  (^*^) :: Num a => f a -> f a -> f a



-- | negate the values in a functor
negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate

-- | subtract two Additive objects
(^-^) :: (Additive f, Num a) => f a -> f a -> f a
x ^-^ y = x ^+^ negated y






-- * Vector space
class Additive f => VectorSpace f where
  -- | multiplication by a scalar
  (.*) :: Num a => a -> f a -> f a
  

-- |linear interpolation
lerp :: (VectorSpace f, Num a) => a -> f a -> f a -> f a
lerp a u v = a .* u ^+^ ((1-a) .* v)


-- * Hilbert space (inner product)
class VectorSpace f => Hilbert f where
  -- | inner product
  dot :: Num a => f a -> f a -> a


-- ** Hilbert-space distance function
-- |`hilbertDistSq x y = || x - y ||^2`
hilbertDistSq :: (Hilbert f, Num a) => f a -> f a -> a
hilbertDistSq x y = dot t t where
  t = x ^-^ y

  


-- * Normed vector space
class Hilbert f => Normed f where
  norm :: (Floating a, Eq a) => a -> f a -> a




-- ** Norms and related results

-- | Squared 2-norm
normSq :: (Hilbert f, Num a) => f a -> a
normSq v = v `dot` v


-- |L1 norm
norm1 :: (Foldable t, Num a, Functor t) => t a -> a
norm1 v = sum (fmap abs v)

-- |Euclidean norm
norm2 :: (Hilbert f, Floating a) => f a -> a
norm2 v = sqrt (normSq v)

-- |Lp norm (p > 0)
normP :: (Foldable t, Functor t, Floating a) => a -> t a -> a
normP p v = sum u**(1/p) where
  u = fmap (**p) v

-- |Infinity-norm
normInfty :: (Foldable t, Ord a) => t a -> a
normInfty = maximum



-- |Normalize w.r.t. p-norm (p finite)
normalize :: (Normed f, Floating a, Eq a) => a -> f a -> f a
normalize p v = (1 / norm p v) .* v







-- |Lp inner product (p > 0)
dotLp :: (Set t, Foldable t, Floating a) => a -> t a -> t a ->  a
dotLp p v1 v2 = sum u**(1/p) where
  f a b = (a*b)**p
  u = liftI2 f v1 v2


-- |Reciprocal
reciprocal :: (Functor f, Fractional b) => f b -> f b
reciprocal = fmap recip


-- |Scale
scale :: (Num b, Functor f) => b -> f b -> f b
scale n = fmap (* n)







-- * FiniteDim : finite-dimensional objects

class Additive f => FiniteDim f where
  type FDSize f :: *
  dim :: f a -> FDSize f


-- | unary dimension-checking bracket
withDim :: (FiniteDim f, Show e) =>
     f a
     -> (FDSize f -> f a -> Bool)
     -> (f a -> c)
     -> String
     -> (f a -> e)
     -> c
withDim x p f e ef | p (dim x) x = f x
                   | otherwise = error e' where e' = e ++ show (ef x)

-- | binary dimension-checking bracket
withDim2 :: (FiniteDim f, FiniteDim g, Show e) =>
     f a
     -> g b
     -> (FDSize f -> FDSize g -> f a -> g b -> Bool)
     -> (f a -> g b -> c)
     -> String
     -> (f a -> g b -> e)
     -> c
withDim2 x y p f e ef | p (dim x) (dim y) x y = f x y
                      | otherwise = error e' where e' = e ++ show (ef x y)






-- * HasData : accessing inner data (do not export)

class Additive f => HasData f a where
  type HDData f a :: * 
  dat :: f a -> HDData f a


-- * Sparse : sparse datastructures

class (FiniteDim f, HasData f a) => Sparse f a where
  spy :: Fractional b => f a -> b




-- * Set : types that behave as sets

class Functor f => Set f where
  -- |union binary lift : apply function on _union_ of two Sets
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- |intersection binary lift : apply function on _intersection_ of two Sets
  liftI2 :: (a -> b -> c) -> f a -> f b -> f c  





class Sparse c a => SpContainer c a where
  type ScIx c :: *
  scInsert :: ScIx c -> a -> c a -> c a
  scLookup :: c a -> ScIx c -> Maybe a
  -- -- | Lookup with default, infix form ("safe" : should throw an exception if lookup is outside matrix bounds)
  (@@) :: c a -> ScIx c -> a





-- * IxContainer : indexed container types

-- class IxContainer (c :: * -> *) a where
--   type Ix c :: *
--   type IxSz c :: *
--   ixcLookup :: c a -> Ix c -> Maybe a
--   ixcIfilter :: (Ix c -> a -> Bool) -> c a -> c a
--   ixcInsert :: Ix c -> a -> c a -> c a
--   ixcFromList :: Foldable t => IxSz c -> t (Ix c, a) -> c a
--   ixcToList :: c a -> [(Ix c, a)]

-- newtype IM_ a = IM (IM.IntMap a)

-- instance IxContainer IM_ a where
--   type Ix IM_  = Int
-- --   -- ixcLookupDefault = lookupDefaultSV
-- --   -- ixcFilter = filterSV


-- newtype IM2 a = IM2 { unIM2 :: IM.IntMap (IM.IntMap a)}

-- instance IxContainer IM2 a where
--   type Ix IM2 = (Int, Int)
--   ixcIfilter f im2 = IM2 $ ifilterIM2 (curry f) (unIM2 im2)



-- class Rank2 (c :: * -> *) a where
--   type R2IxRow c :: *
--   type R2IxCol c :: *
--   type R2V c :: * -> *
--   rank2Act :: c a -> R2V c a -> R2V c a
--   extractRow :: c a -> R2IxRow c -> Maybe (R2V c a)
--   extractCol :: c a -> R2IxCol c -> Maybe (R2V c a)
--   -- extractRows :: c a -> [R2IxRow c] -> [R2V c]
--   -- extractCols :: c a -> [R2IxCol c] -> [R2V c]
  




-- * SMatrix : sparse matrix types

-- class (IxContainer c a, Sparse c a, Additive c) => SMatrix c a where
