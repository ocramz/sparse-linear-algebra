{-# language TypeFamilies, MultiParamTypeClasses, KindSignatures, FlexibleContexts, FlexibleInstances, ConstraintKinds #-}
{-# language AllowAmbiguousTypes #-}
{-# language CPP #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Class
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  GPL-3 (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- Typeclasses for linear algebra and related concepts
--
-----------------------------------------------------------------------------
module Numeric.LinearAlgebra.Class where

-- import Control.Applicative
import Data.Complex

-- import Control.Exception
-- import Control.Exception.Common
import Control.Monad.Catch
import Control.Monad.IO.Class

-- import Data.Typeable (Typeable)

import qualified Data.Vector as V (Vector)

import Data.Sparse.Types
import Numeric.Eps



-- * Matrix and vector elements (optionally Complex)
class (Eq e , Fractional e, Floating e, Num (EltMag e), Ord (EltMag e)) => Elt e where
  type EltMag e :: *
  -- | Complex conjugate, or identity function if its input is real-valued
  conj :: e -> e
  conj = id
  -- | Magnitude
  mag :: e -> EltMag e

instance Elt Double where {type EltMag Double = Double ; mag = id}
instance Elt Float where {type EltMag Float = Float; mag = id}
instance (RealFloat e) => Elt (Complex e) where
  type EltMag (Complex e) = e
  conj = conjugate
  mag = magnitude
  



infixl 6 ^+^, ^-^

-- * Additive group
class AdditiveGroup v where
  -- | The zero element: identity for '(^+^)'
  zeroV :: v
  -- | Add vectors
  (^+^) :: v -> v -> v
  -- | Additive inverse
  negateV :: v -> v
  -- | Group subtraction
  (^-^) :: v -> v -> v
  (^-^) x y = x ^+^ negateV y


infixr 7 .*

-- * Vector space @v@.
class (AdditiveGroup v, Num (Scalar v)) => VectorSpace v where
  type Scalar v :: *
  -- | Scale a vector
  (.*) :: Scalar v -> v -> v

-- | Adds inner (dot) products.
class (VectorSpace v, AdditiveGroup (Scalar v)) => InnerSpace v where
  -- | Inner/dot product
  (<.>) :: v -> v -> Scalar v

-- | Inner product
dot :: InnerSpace v => v -> v -> Scalar v
dot = (<.>)
  

infixr 7 ./
infixl 7 *.

-- | Scale a vector by the reciprocal of a number (e.g. for normalization)
(./) :: (VectorSpace v, s ~ Scalar v, Fractional s) => v -> s -> v
v ./ s = (recip s) .* v

-- | Vector multiplied by scalar
(*.) :: (VectorSpace v, s ~ Scalar v) => v -> s -> v
(*.) = flip (.*)  



-- | Convex combination of two vectors (NB: 0 <= `a` <= 1). 
cvx :: VectorSpace v => Scalar v -> v -> v -> v
cvx a u v = a .* u ^+^ ((1-a) .* v)




-- ** Hilbert-space distance function

-- |`hilbertDistSq x y = || x - y ||^2` computes the squared L2 distance between two vectors
hilbertDistSq :: InnerSpace v => v -> v -> Scalar v
hilbertDistSq x y = t <.> t where
  t = x ^-^ y








-- * Normed vector spaces

class (InnerSpace v, Num (RealScalar v), Eq (RealScalar v), Epsilon (Magnitude v), Show (Magnitude v), Ord (Magnitude v)) => Normed v where
  type Magnitude v :: *
  type RealScalar v :: *
  -- | L1 norm
  norm1 :: v -> Magnitude v
  -- | Euclidean (L2) norm squared
  norm2Sq :: v -> Magnitude v
  -- | Lp norm (p > 0)
  normP :: RealScalar v -> v -> Magnitude v
  -- | Normalize w.r.t. Lp norm
  normalize :: RealScalar v -> v -> v
  -- | Normalize w.r.t. L2 norm
  normalize2 :: v -> v
  -- | Normalize w.r.t. norm2' instead of norm2
  normalize2' :: Floating (Scalar v) => v -> v 
  normalize2' x = x ./ norm2' x
  -- | Euclidean (L2) norm
  norm2 :: Floating (Magnitude v) => v -> Magnitude v 
  norm2 x = sqrt (norm2Sq x)
  -- | Euclidean (L2) norm; returns a Complex (norm :+ 0) for Complex-valued vectors
  norm2' :: Floating (Scalar v) => v -> Scalar v 
  norm2' x = sqrt $ x <.> x
  -- | Lp norm (p > 0)
  norm :: Floating (Magnitude v) => RealScalar v -> v -> Magnitude v 
  norm p v
    | p == 1 = norm1 v
    | p == 2 = norm2 v
    | otherwise = normP p v



-- | Infinity-norm (Real)
normInftyR :: (Foldable t, Ord a) => t a -> a
normInftyR x = maximum x

-- | Infinity-norm (Complex)
normInftyC :: (Foldable t, RealFloat a, Functor t) => t (Complex a) -> a
normInftyC x = maximum (magnitude <$> x)





-- instance Normed Double where
--   type Magnitude Double = Double
--   type RealScalar Double = Double
--   norm1 = abs
--   norm2Sq = abs
--   normP _ = abs
--   normalize _ _ = 1
--   normalize2 _ = 1
--   norm2 = abs
--   norm2' = abs

-- instance Normed (Complex Double) where
--   type Magnitude (Complex Double) = Double
--   type RealScalar (Complex Double) = Double
--   norm1 (r :+ i) = abs r + abs i
--   norm2Sq = (**2) . magnitude
--   normP p (r :+ i) = (r**p + i**p)**(1/p)
--   normalize p c = c ./ normP p c
--   normalize2 c = c ./ magnitude c
--   norm2 = magnitude
--   norm2' = magnitude
  




    






-- | Lp inner product (p > 0)
dotLp :: (Set t, Foldable t, Floating a) => a -> t a -> t a ->  a
dotLp p v1 v2 = sum u**(1/p) where
  f a b = (a*b)**p
  u = liftI2 f v1 v2


-- | Reciprocal
reciprocal :: (Functor f, Fractional b) => f b -> f b
reciprocal = fmap recip


-- |Scale
scale :: (Num b, Functor f) => b -> f b -> f b
scale n = fmap (* n)










-- * Matrix ring

-- | A matrix ring is any collection of matrices over some ring R that form a ring under matrix addition and matrix multiplication

class (AdditiveGroup m, Epsilon (MatrixNorm m)) => MatrixRing m where
  type MatrixNorm m :: *
  -- | Matrix-matrix product
  (##) :: m -> m -> m
  -- | Matrix times matrix transpose (A B^T)
  (##^) :: m -> m -> m
  -- | Matrix transpose times matrix (A^T B)
  (#^#) :: m -> m -> m
  a #^# b = transpose a ## b
  -- | Matrix transpose (Hermitian conjugate in the Complex case)
  transpose :: m -> m
  -- | Frobenius norm
  normFrobenius :: m -> MatrixNorm m



-- a "sparse matrix ring" ?

-- class MatrixRing m a => SparseMatrixRing m a where
--   (#~#) :: Epsilon a => Matrix m a -> Matrix m a -> Matrix m a







-- * Linear vector space

class (VectorSpace v, MatrixRing (MatrixType v)) => LinearVectorSpace v where
  type MatrixType v :: *
  -- | Matrix-vector action
  (#>) :: MatrixType v -> v -> v
  -- | Dual matrix-vector action
  (<#) :: v -> MatrixType v -> v





-- ** LinearVectorSpace + Normed

type V v = (LinearVectorSpace v, Normed v) 




-- ** Linear systems
  
class LinearVectorSpace v => LinearSystem v where
  -- | Solve a linear system; uses GMRES internally as default method
  (<\>) :: (MonadIO m, MonadThrow m) =>
           MatrixType v   -- ^ System matrix
        -> v              -- ^ Right-hand side
        -> m v            -- ^ Result










-- * FiniteDim : finite-dimensional objects

class FiniteDim f where
  type FDSize f
  -- | Dimension (i.e. Int for SpVector, (Int, Int) for SpMatrix)
  dim :: f -> FDSize f




-- * HasData : accessing inner data (do not export)

class HasData f where
  type HDData f 
  -- | Number of nonzeros
  nnz :: f -> Int
  dat :: f -> HDData f



-- * Sparse : sparse datastructures

class (FiniteDim f, HasData f) => Sparse f where
  -- | Sparsity (fraction of nonzero elements)
  spy :: Fractional b => f -> b



-- * Set : types that behave as sets

class Functor f => Set f where
  -- | Union binary lift : apply function on _union_ of two "sets"
  liftU2 :: (a -> a -> a) -> f a -> f a -> f a

  -- | Intersection binary lift : apply function on _intersection_ of two "sets"
  liftI2 :: (a -> a -> b) -> f a -> f a -> f b




-- * SpContainer : sparse container datastructures. Insertion, lookup, toList, lookup with 0 default
class Sparse c => SpContainer c where
  type ScIx c :: *
  type ScElem c
  scInsert :: ScIx c -> ScElem c -> c -> c
  scLookup :: c -> ScIx c -> Maybe (ScElem c)
  scToList :: c -> [(ScIx c, ScElem c)]
  -- -- | Lookup with default, infix form ("safe" : should throw an exception if lookup is outside matrix bounds)
  (@@) :: c -> ScIx c -> ScElem c







-- * SparseVector

class SpContainer v => SparseVector v where
  type SpvIx v :: *
  svFromList :: Int -> [(SpvIx v, ScElem v)] -> v
  svFromListDense :: Int -> [ScElem v] -> v
  svConcat :: Foldable t => t v -> v
  -- svZipWith :: (e -> e -> e) -> v e -> v e -> v e








-- * SparseMatrix

class SpContainer m => SparseMatrix m where
  smFromVector :: LexOrd -> (Int, Int) -> V.Vector (IxRow, IxCol, ScElem m) -> m
  -- smFromFoldableDense :: Foldable t => t e -> m e  
  smTranspose :: m -> m
  -- smExtractSubmatrix ::
  --   m e -> (IxRow, IxRow) -> (IxCol, IxCol) -> m e
  encodeIx :: m -> LexOrd -> (IxRow, IxCol) -> LexIx
  decodeIx :: m -> LexOrd -> LexIx -> (IxRow, IxCol)


-- data RowsFirst = RowsFirst
-- data ColsFirst = ColsFirst








-- * SparseMatVec

-- | Combining functions for relating (structurally) matrices and vectors, e.g. extracting/inserting rows/columns/submatrices

-- class (SparseMatrix m o e, SparseVector v e) => SparseMatVec m o v e where
--   smvInsertRow :: m e -> v e -> IxRow -> m e
--   smvInsertCol :: m e -> v e -> IxCol -> m e
--   smvExtractRow :: m e -> IxRow -> v e
--   smvExtractCol :: m e -> IxCol -> v e  



-- * Utilities


-- | Lift a real number onto the complex plane
toC :: Num a => a -> Complex a
toC r = r :+ 0















-- -- | Instances for AdditiveGroup
-- instance Integral a => AdditiveGroup (Ratio a) where
--   {zero=0; (^+^) = (+); negated = negate}

-- instance (RealFloat v, AdditiveGroup v) => AdditiveGroup (Complex v) where
--   zero    = zero :+ zero
--   (^+^)   = (+)
--   negated = negate

-- -- | Standard instance for an applicative functor applied to a vector space.
-- instance AdditiveGroup v => AdditiveGroup (a -> v) where
--   zero    = pure   zero
--   (^+^)   = liftA2 (^+^)
--   negated = fmap   negated


-- -- | Instances for VectorSpace
-- instance (RealFloat v, VectorSpace v) => VectorSpace (Complex v) where
--   type Scalar (Complex v) = Scalar v
--   s .* (u :+ v) = s .* u :+ s .* v



-- #define ScalarType(t) \
--   instance AdditiveGroup (t) where {zero = 0; (^+^) = (+); negated = negate};\
--   instance VectorSpace (t) where {type Scalar (t) = (t); (.*) = (*) };\
--   instance Hilbert (t) where dot = (*)

-- ScalarType(Int)
-- ScalarType(Integer)
-- ScalarType(Float)
-- ScalarType(Double)
-- ScalarType(CSChar)
-- ScalarType(CInt)
-- ScalarType(CShort)
-- ScalarType(CLong)
-- ScalarType(CLLong)
-- ScalarType(CIntMax)
-- ScalarType(CFloat)
-- ScalarType(CDouble)
