{-# language TypeFamilies, FlexibleInstances, DeriveFunctor #-}
module Data.Sparse.Internal.TriMatrix where

-- import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as IM
import Data.IntMap.Strict ((!))
-- import qualified Data.Set as S
-- import qualified Data.Vector as V

import Data.Maybe (fromMaybe)
-- import Data.Monoid
-- import Data.Complex

import Numeric.Eps
import Data.Sparse.Types
-- import Data.Sparse.Utils
import Data.Sparse.Internal.SList

import Data.VectorSpace
import Numeric.LinearAlgebra.Class

{- | triangular sparse matrix, row-major order

Intmap-of-sparse lists
* fast random access of rows
* fast consing of row elements
-}

newtype TriMatrix a = TM { unTM :: IM.IntMap (SList a)} deriving (Show, Functor)

emptyTM :: Int -> TriMatrix a
emptyTM n = TM $ IM.fromList [(i, emptySL) | i <- [0 .. n-1]]

-- | `appendIM i x im` appends an element `x` to the i'th SList in an IntMap-of-SLists structure
appendIM :: IM.Key -> (Int, a) -> IM.IntMap (SList a) -> IM.IntMap (SList a)
appendIM i x im = IM.insert i (x `consSL` e) im where
  e = fromMaybe emptySL (IM.lookup i im)

-- -- | Appends a column to a TriMatrix from an IntMap
-- appendColTM j z imx tm = IM.foldlWithKey ins z imx where
--   ins acc i x = appendIM tm (j, x) acc
  
solveUij i j aij lmat umat = aij - li <.> uj
  where
    li = unTM lmat ! i
    uj = unTM umat ! j

solveLij i j aij lmat umat = (aij - li <.> uj)/utt
  where
    li = unTM lmat ! i
    uj = unTM umat ! j
    utt = snd $ head $ unSL uj




-- -- innerMaybe :: Epsilon a => SVector a -> SVector a -> Maybe a
-- innerMaybe :: (Epsilon a, Ord i) => [(i, a)] -> [(i, a)] -> Maybe a
-- innerMaybe u v | isNz x = Just x
--                | otherwise = Nothing where x = inner u v

 










{- | LU factorization : store L and U^T in TriMatrix format

* 

-}





