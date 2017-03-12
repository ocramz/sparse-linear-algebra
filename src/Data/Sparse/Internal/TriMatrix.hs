{-# language TypeFamilies, FlexibleInstances, DeriveFunctor #-}
module Data.Sparse.Internal.TriMatrix where

-- import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as IM
import Data.IntMap.Strict ((!))
-- import qualified Data.Set as S
-- import qualified Data.Vector as V

import Data.Maybe (fromMaybe)
-- import Data.Monoid
import Data.Ord (comparing)
import Data.Complex

import Numeric.Eps
import Data.Sparse.Types
import Data.Sparse.Utils

import Data.VectorSpace
import Numeric.LinearAlgebra.Class

{- | triangular sparse matrix, row-major order

Intmap-of-lists
* fast random access of rows
* fast consing of row elements
-}
-- newtype TriMatrix a = TM { unTM :: IM.IntMap [(Int, a)]} deriving (Show, Functor)

newtype TriMatrix a = TM { unTM :: IM.IntMap (SList a)} deriving (Show, Functor)

emptyTM :: Int -> TriMatrix a
emptyTM n = TM $ IM.fromList [(i, emptySL) | i <- [0 .. n-1]]

-- | `appendIM i x im` appends an element `x` to the i'th list in an IntMap-of-lists structure

appendIM :: IM.Key -> (Int, a) -> IM.IntMap (SList a) -> IM.IntMap (SList a)
appendIM i x im = IM.insert i (x `consSL` e) im where
  e = fromMaybe emptySL (IM.lookup i im)


-- appendToRowTM :: IxRow -> IxCol -> a -> TriMatrix a -> TriMatrix a
-- appendToRowTM i j x tm = TM $ appendIM i (j, x) (unTM tm)





-- -- innerMaybe :: Epsilon a => SVector a -> SVector a -> Maybe a
-- innerMaybe :: (Epsilon a, Ord i) => [(i, a)] -> [(i, a)] -> Maybe a
-- innerMaybe u v | isNz x = Just x
--                | otherwise = Nothing where x = inner u v

 















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





{- | LU factorization : store L and U^T in TriMatrix format

* 

-}



-- test data

l1, l2 :: [(Int, Double)]
l1 = [(0, pi), (2, pi),  (3, 5.4) ]

l2 = [(1, exp 1), (2, 3.4)]

l1c :: [(Int, Complex Double)]
l1c = zip ii $ zipWith (:+) [1..3] [3,2,1] where
  ii = [0, 2, 5]

sl1c = SL l1c


-- helpers
sortIndices :: [(IM.Key, a)] -> [(IM.Key, a)]
sortIndices = IM.toList . IM.fromList

