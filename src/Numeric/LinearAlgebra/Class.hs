{-# language TypeFamilies, MultiParamTypeClasses, KindSignatures, FlexibleContexts, FlexibleInstances #-}
{-# language CPP #-}
module Numeric.LinearAlgebra.Class where

import Control.Applicative
import Data.Complex
import Data.Ratio
import Foreign.C.Types (CSChar, CInt, CShort, CLong, CLLong, CIntMax, CFloat, CDouble)

-- * Matrix and vector elements (possibly Complex)
class (Eq e , Fractional e) => Elt e where
  conj :: e -> e
  conj = id
instance Elt Double
instance Elt Float
instance RealFloat e => Elt (Complex e) where
  conj = conjugate
  

-- * Additive group
class AdditiveGroup e where
  -- | Identity element
  zero :: e
  -- | Group action
  (^+^) :: e -> e -> e
  -- | Inverse element
  negated :: e -> e
  -- | Inverse group action
  (^-^) :: e -> e -> e
  x ^-^ y = x ^+^ negated y

-- | Numeric isntances for AdditiveGroup

-- For 'Num' types:
-- 
-- instance AdditiveGroup n where {zeroV=0; (^+^) = (+); negateV = negate}




instance Integral a => AdditiveGroup (Ratio a) where
  {zero=0; (^+^) = (+); negated = negate}
instance (RealFloat v, AdditiveGroup v) => AdditiveGroup (Complex v) where
  zero   = zero :+ zero
  (^+^)   = (+)
  negated = negate

-- | Standard instance for an applicative functor applied to a vector space.
instance AdditiveGroup v => AdditiveGroup (a -> v) where
  zero   = pure   zero
  (^+^)   = liftA2 (^+^)
  negated = fmap   negated



-- * Vector space
class AdditiveGroup v => VectorSpace v where
  type Scalar v :: *
  -- | Scale a vector
  (.*) :: Scalar v -> v -> v  

-- |Linear interpolation
lerp :: (VectorSpace e, Num (Scalar e)) => Scalar e -> e -> e -> e
lerp a u v = a .* u ^+^ ((1-a) .* v)


-- linearCombination :: (VectorSpace v , Foldable t) => t (Scalar v, v) -> v
-- linearCombination  =  foldr (\(a, x) (b, y) -> (a .* x) ^+^ (b .* y)) 

-- | numerical instances for VectorSpace
instance (RealFloat v, VectorSpace v) => VectorSpace (Complex v) where
  type Scalar (Complex v) = Scalar v
  s .* (u :+ v) = s .* u :+ s .* v




-- * Hilbert space (inner product)
infixr 7 `dot`
class (VectorSpace v, AdditiveGroup (Scalar v)) => Hilbert v where
  dot :: v -> v -> Scalar v
  
-- infixr 7 <.>
-- (<.>) = dot  





-- ** Hilbert-space distance function
-- |`hilbertDistSq x y = || x - y ||^2`
hilbertDistSq :: Hilbert v => v -> v -> Scalar v
hilbertDistSq x y = dot t t where
  t = x ^-^ y








-- * Normed vector space
class Hilbert e => Normed e where
  -- |p-norm (p finite)
  norm :: RealFloat a => a -> e -> Scalar e
  -- |Normalize w.r.t. p-norm
  normalize :: RealFloat a => a -> e -> e




-- ** Norms and related results

-- | Squared 2-norm
normSq :: Hilbert v => v -> Scalar v
normSq v = v `dot` v


-- |L1 norm
norm1 :: (Foldable t, Num a, Functor t) => t a -> a
norm1 v = sum (fmap abs v)

-- |Euclidean norm
norm2 :: (Hilbert v, Floating (Scalar v)) => v -> Scalar v
norm2 v = sqrt (normSq v)

-- |Lp norm (p > 0)
-- normP :: (Hilbert f e, Floating a) => a -> f e -> HT e
normP :: (Foldable t, Functor t, Floating a) => a -> t a -> a
normP p v = sum u**(1/p) where
  u = fmap (**p) v

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

class Functor f => FiniteDim f where
  type FDSize f :: *
  dim :: f a -> FDSize f


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

class HasData f a where
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
  liftI2 :: (a -> a -> b) -> f a -> f a -> f b





class Sparse c a => SpContainer c a where
  type ScIx c :: *
  scInsert :: ScIx c -> a -> c a -> c a
  scLookup :: c a -> ScIx c -> Maybe a
  scToList :: c a -> [a]
  -- -- | Lookup with default, infix form ("safe" : should throw an exception if lookup is outside matrix bounds)
  (@@) :: c a -> ScIx c -> a









-- * SparseVector

class (SpContainer v e, Hilbert e) => SparseVector v e where
  type SpvIx v :: *
  svFromList :: Int -> [(SpvIx v, e)] -> v e
  svFromListDense :: Int -> [e] -> v e
  svConcat :: Foldable t => t (v e) -> v e
  -- svZipWith :: (e -> e -> e) -> v e -> v e -> v e

-- * SparseMatrix

class (SpContainer m e, AdditiveGroup (m e)) => SparseMatrix m e where
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









#define ScalarType(t) \
  instance AdditiveGroup (t) where {zero = 0; (^+^) = (+); negated = negate};\
  instance VectorSpace (t) where {type Scalar (t) = (t); (.*) = (*) };\
  instance Hilbert (t) where dot = (*)

ScalarType(Int)
ScalarType(Integer)
ScalarType(Float)
ScalarType(Double)
ScalarType(CSChar)
ScalarType(CInt)
ScalarType(CShort)
ScalarType(CLong)
ScalarType(CLLong)
ScalarType(CIntMax)
ScalarType(CFloat)
ScalarType(CDouble)
