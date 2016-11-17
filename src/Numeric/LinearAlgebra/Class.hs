{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, KindSignatures, FlexibleContexts, FlexibleInstances #-}
module Numeric.LinearAlgebra.Class where

import Data.Complex

-- * Matrix and vector elements (possibly Complex)
class (Eq e , Fractional e) => Elt e where
  conj :: e -> e
  conj = id
  
instance Elt Double
instance Elt Float
instance RealFloat e => Elt (Complex e) where
  conj = conjugate
  


-- * Additive monoid
class Functor f => Additive f where
  -- | Monoid neutral element
  zero :: Num e => f e
  
  -- | Monoid +
  (^+^) :: Num e => f e -> f e -> f e

  -- one :: Num a => f a
  -- (^*^) :: Num a => f a -> f a -> f a



-- | negate the values in a functor
negated :: (Num a, Functor f) => f a -> f a
negated = fmap negate

-- | subtract two Additive objects
(^-^) :: (Num e , Additive f) => f e -> f e -> f e
x ^-^ y = x ^+^ negated y






-- * Vector space
class (Elt e , Additive f) => VectorSpace f e where
  -- | multiplication by a scalar
  (.*) :: e -> f e -> f e
  

-- |linear interpolation
lerp :: VectorSpace f e => e -> f e -> f e -> f e
lerp a u v = a .* u ^+^ ((1-a) .* v)









-- * Hilbert space (inner product)
class VectorSpace f e => Hilbert f e where
  type HT e :: *
  type instance HT e = Double
  -- | inner product
  dot :: f e -> f e -> HT e



-- ** Hilbert-space distance function
-- |`hilbertDistSq x y = || x - y ||^2`
hilbertDistSq :: Hilbert f e => f e -> f e -> HT e
hilbertDistSq x y = dot t t where
  t = x ^-^ y








-- * Normed vector space
class Hilbert f e => Normed f e where
  -- |p-norm (p finite)
  norm :: RealFloat a => a -> f e -> HT e
  -- |Normalize w.r.t. p-norm
  normalize :: RealFloat a => a -> f e -> f e




-- ** Norms and related results

-- | Squared 2-norm
normSq :: Hilbert f e => f e -> HT e
normSq v = v `dot` v


-- |L1 norm
norm1 :: (Foldable t, Num a, Functor t) => t a -> a
norm1 v = sum (fmap abs v)

-- |Euclidean norm
norm2 :: (Hilbert f e, Floating (HT e)) => f e -> HT e
norm2 v = sqrt (normSq v)

-- |Lp norm (p > 0)
-- normP :: (Hilbert f e, Floating a) => a -> f e -> HT e
-- normP p v = sum u**(1/p) where
--   u = fmap (**p) v

-- |Infinity-norm
normInfty :: (Foldable t, Ord a) => t a -> a
normInfty = maximum










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
  dim :: f e -> FDSize f


-- | unary dimension-checking bracket
withDim :: (FiniteDim f, Show s) =>
     f e
     -> (FDSize f -> f e -> Bool)
     -> (f e -> c)
     -> String
     -> (f e -> s)
     -> c
withDim x p f e ef | p (dim x) x = f x
                   | otherwise = error e' where e' = e ++ show (ef x)

-- | binary dimension-checking bracket
withDim2 :: (FiniteDim f, FiniteDim g, Show s) =>
     f e
     -> g e
     -> (FDSize f -> FDSize g -> f e -> g e -> Bool)
     -> (f e -> g e -> c)
     -> String
     -> (f e -> g e -> s)
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
  scToList :: c a -> [a]
  -- -- | Lookup with default, infix form ("safe" : should throw an exception if lookup is outside matrix bounds)
  (@@) :: c a -> ScIx c -> a









-- * SparseVector

class (SpContainer v e, Hilbert v e) => SparseVector v e where
  type SpvIx v :: *
  svFromList :: Int -> [(SpvIx v, e)] -> v e
  svFromListDense :: Int -> [e] -> v e
  svConcat :: Foldable t => t (v e) -> v e
  -- svZipWith :: (e -> e -> e) -> v e -> v e -> v e

-- * SparseMatrix

class (SpContainer m e, Additive m) => SparseMatrix m e where
  type SpmIxRow m :: *
  type SpmIxCol m :: *  
  smFromFoldable :: Foldable t => (Int, Int) -> t (SpmIxRow m, SpmIxCol m, e) -> m e
  smFromFoldableDense :: Foldable t => t e -> m e  
  smTranspose :: m e -> m e
  smExtractSubmatrix ::
    m e -> (SpmIxRow m, SpmIxRow m) -> (SpmIxCol m, SpmIxCol m) -> m e


class (SparseMatrix m e, SparseVector v e) => LinearSpace m v e where
  lsInsertRow :: m e -> v e -> SpmIxRow m -> m e
  lsInsertCol :: m e -> v e -> SpmIxCol m -> m e
  lsExtractRow :: m e -> SpmIxRow m -> v e
  lsExtractCol :: m e -> SpmIxCol m -> v e  






-- class (Elt e, SparseVector v e) => DOT s v e where
--   dot' :: StratCxt s -> v e -> v e -> ResS s e

class (Elt e, SparseVector v e) => VectorSpace' v e where
  vscale :: e -> v e -> v e

class (VectorSpace' v e) => Hilbert' v e where
  vdot :: v e -> v e -> e

class Hilbert' v e => Normed' v e where
  vnorm :: Floating a => Int -> v e -> a
