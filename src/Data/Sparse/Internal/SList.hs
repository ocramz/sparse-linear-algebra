{-# language TypeFamilies, FlexibleInstances, DeriveFunctor #-}
module Data.Sparse.Internal.SList where


import Data.Sparse.Utils
import Numeric.LinearAlgebra.Class

import Data.VectorSpace




-- | Sparse list
newtype SList a = SL {unSV :: [(Int, a)]} deriving (Eq, Show, Functor)

emptySL :: SList a
emptySL = SL []

consSL :: (Int, a) -> SList a -> SList a
consSL x (SL xs) = SL (x : xs)

headSL :: SList a -> Maybe (Int, a)
headSL (SL (x:_)) = Just x
headSL (SL []) = Nothing

fromList :: [(Int, a)] -> SList a
fromList = SL

toList :: SList a -> [(Int, a)]
toList = unSV


{-|
NB : unionWith and intersectWith work only if the indices are _sorted_
-}


-- | Inner product between sparse lists
sldot :: (Elt a, Ord i) => [(i, a)] -> [(i, a)] -> a
sldot u v = sum $ intersectWith pf u v where
  pf x y = conj x * y

-- | Vector sum of sparse lists
slsum :: (Ord i, Elt a) => [(i, a)] -> [(i, a)] -> [(i, a)]
slsum = unionWith (+) 0




-- | `vector-space` instances

instance Elt a => AdditiveGroup (SList a) where
  zeroV = SL []
  negateV = fmap (* (-1))
  u ^+^ v = SL $ slsum (unSV u) (unSV v)


instance Elt a => VectorSpace (SList a) where
  type Scalar (SList a) = a
  a *^ v = fmap (* a) v


instance (AdditiveGroup a, Elt a) => InnerSpace (SList a) where
  u <.> v = sldot (unSV u) (unSV v)  


-- instance InnerSpace (SList Double) where
--   u <.> v = inner (unSV u) (unSV v)

-- instance InnerSpace (SList (Complex Double)) where
--   u <.> v = innerC (unSV u) (unSV v)
