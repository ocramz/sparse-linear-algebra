{-# language TypeFamilies, FlexibleInstances, DeriveFunctor #-}
module Data.Sparse.Internal.SList where


import Data.Sparse.Utils
import Numeric.LinearAlgebra.Class




-- | Sparse list
newtype SList a = SL {unSL :: [(Int, a)]} deriving (Eq, Show, Functor)

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
toList = unSL


{-|
NB : unionWith and intersectWith work only if the indices are _sorted_

NB2 : we use the _descending_ order comparison
-}


-- | Inner product between sparse lists
sldot :: (Elt a, Ord i) => [(i, a)] -> [(i, a)] -> a
sldot u v = sum $ intersectWithD pf u v where
  pf x y = conj x * y

-- | Vector sum of sparse lists
slsum :: (Ord i, Elt a) => [(i, a)] -> [(i, a)] -> [(i, a)]
slsum = unionWithD (+) 0









-- | `vector-space` instances

instance Elt a => AdditiveGroup (SList a) where
  zeroV = SL []
  negateV = fmap (* (-1))
  u ^+^ v = SL $ slsum (unSL u) (unSL v)


instance Elt a => VectorSpace (SList a) where
  type Scalar (SList a) = a
  a .* v = fmap (* a) v


instance (AdditiveGroup a, Elt a) => InnerSpace (SList a) where
  u <.> v = sldot (unSL u) (unSL v)  


-- instance InnerSpace (SList Double) where
--   u <.> v = inner (unSV u) (unSV v)

-- instance InnerSpace (SList (Complex Double)) where
--   u <.> v = innerC (unSV u) (unSV v)




-- test data

l1, l2 :: [(Int, Double)]
l1 = [(0, pi), (2, pi),  (3, 5.4) ]

l2 = [(1, exp 1), (2, 3.4)]

-- l1c :: [(Int, Complex Double)]
-- l1c = zip ii $ zipWith (:+) [1..3] [3,2,1] where
--   ii = [0, 2, 5]

-- sl1c = SL l1c


-- helpers
-- sortIndices :: [(IM.Key, a)] -> [(IM.Key, a)]
-- sortIndices = IM.toList . IM.fromList
