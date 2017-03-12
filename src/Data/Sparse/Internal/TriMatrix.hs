{-# language TypeFamilies, FlexibleInstances, DeriveFunctor #-}
module Data.Sparse.Internal.TriMatrix where

-- import qualified Data.Map.Strict as M
import qualified Data.IntMap.Strict as IM
import Data.IntMap.Strict ((!))
-- import qualified Data.Set as S
-- import qualified Data.Vector as V

import Data.Foldable (foldrM)
import Data.Maybe (fromMaybe)
-- import Data.Monoid
-- import Data.Complex

import Numeric.Eps
import Data.Sparse.Types
-- import Data.Sparse.Utils
import Data.Sparse.Internal.SList

import Data.VectorSpace
import Numeric.LinearAlgebra.Class
-- import Data.Sparse.SpMatrix ((@@!))
import qualified Data.Sparse.Internal.IntM as IntM
import Data.Sparse.Common ((@@!), nrows, ncols, lookupSM, SpMatrix)

import Control.Monad.Catch
import Control.Exception.Common

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
  




-- | Nested lookup with default value = 0
lookupWD :: Num a =>
     (irow -> mat -> Maybe row)    -- ^ row lookup
     -> (jcol -> row -> Maybe a)   -- ^ in-row lookup
     -> mat                
     -> irow
     -> jcol
     -> a
lookupWD rlu clu aa i j = fromMaybe 0 (rlu i aa >>= clu j)



-- -- innerMaybe :: Epsilon a => SVector a -> SVector a -> Maybe a
-- innerMaybe :: (Epsilon a, Ord i) => [(i, a)] -> [(i, a)] -> Maybe a
-- innerMaybe u v | isNz x = Just x
--                | otherwise = Nothing where x = inner u v

 










{- | LU factorization : store L and U^T in TriMatrix format -}


-- luStep amat (lmat, umat, t) = 
--   let umat' = solveUrow amat lmat umat t
--       lmat' = solveLcol amat lmat umat' t
--   in (lmat', umat', succ t)



-- solveUrowM amat lmat umat i = foldrM ins umat [i .. n-1] where
--   n = ncols amat
--   li = lmat ! i
--   ins j acc
--        | nearZero x && i == j =
--            throwM (NeedsPivoting "solveUrowM" "bla")
--             -- throwM (NeedsPivoting "solveUrowM" (unwords ["U",show (i,i :: Int)]))
--        | isNz x = pure $ appendIM j (i, x) acc
--        | otherwise = pure acc where
--     x = aij - li <.> uj 
--     aij = amat @@! (i,j)
--     uj = umat ! j

solveUrow
  :: (Elt a, Epsilon a) =>
     SpMatrix a
     -> IM.IntMap (SList a)
     -> IM.IntMap (SList a)
     -> IM.Key
     -> IM.IntMap (SList a)   -- ^ updated U
solveUrow amat lmat umat i = foldr ins umat [i .. n-1] where
  n = ncols amat
  li = lmat ! i
  ins j acc | isNz x = appendIM j (i, x) acc
            | otherwise = acc where
    x = aij - li <.> uj 
    aij = amat @@! (i,j)
    uj = umat ! j

solveLcol
  :: (Elt a, Epsilon a) =>
     SpMatrix a
     -> IM.IntMap (SList a)
     -> IM.IntMap (SList a)
     -> a                   -- ^ diagonal element of U (must be nonzero)
     -> IM.Key
     -> IM.IntMap (SList a) -- ^ updated L
solveLcol amat lmat umat utt j = foldr ins lmat [j .. m-1] where
  m = nrows amat
  uj = umat ! j
  ins i acc
    | i == j = appendIM i (j, 1) acc  -- write 1 on the diagonal 
    | isNz x = appendIM i (j, x) acc
    | otherwise = acc where
    x = (aij - li <.> uj)/utt 
    aij = amat @@! (i,j)
    li = lmat ! i






