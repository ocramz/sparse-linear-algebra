{-# LANGUAGE FlexibleContexts #-}
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
import Data.Sparse.Common ((@@!), nrows, ncols, lookupSM, extractRow, extractCol, SpVector, SpMatrix, foldlWithKeySV)

import Control.Monad.Catch
import Control.Exception.Common

import Control.Monad (when)
import Control.Monad.Trans.State

{- | triangular sparse matrix, row-major order

Intmap-of-sparse lists
* fast random access of rows
* fast consing of row elements
-}

newtype TriMatrix a = TM { unTM :: IM.IntMap (SList a)} deriving (Show, Functor)

emptyIMSL :: Int -> IM.IntMap (SList a)
emptyIMSL n = IM.fromList [(i, emptySL) | i <- [0 .. n-1]]

emptyTM :: Int -> TriMatrix a
emptyTM n = TM (emptyIMSL n)

-- | `appendIM i x im` appends an element `x` to the i'th SList in an IntMap-of-SLists structure
appendIM :: IM.Key -> (Int, a) -> IM.IntMap (SList a) -> IM.IntMap (SList a)
appendIM i x im = IM.insert i (x `consSL` e) im where
  e = fromMaybe emptySL (IM.lookup i im)


-- | Nested lookup with default value = 0
lookupWD :: Num a =>
     (irow -> mat -> Maybe row)    -- ^ row lookup
     -> (jcol -> row -> Maybe a)   -- ^ in-row lookup
     -> mat                
     -> irow
     -> jcol
     -> a
lookupWD rlu clu aa i j = fromMaybe 0 (rlu i aa >>= clu j)

 







{- | LU factorization : store L and U^T in TriMatrix format -}

lu :: (Scalar (SpVector t) ~ t, Elt t, VectorSpace (SpVector t),
      MonadThrow m, Epsilon t) =>
     SpMatrix t
     -> m (IM.IntMap (SList t), IM.IntMap (SList t))
lu amat = do
  (lfin, ufin, ifin) <- execStateT (luStep amat) (luInit amat)
  let (ufin', _) = uStep amat lfin ufin ifin
  return (lfin, ufin')

luInit amat = (lmat0, umat0, 1)
  where
    (m,n) = (nrows amat, ncols amat)
    urow0 = extractRow amat 0
    lcol0 = extractCol amat 0 ./ (urow0 @@ 0) 
    umat0 = foldlWithKeySV ins (emptyIMSL n) lcol0 
    lmat0 = IM.insert 0 (SL [(0, 1)]) l0 where
      l0 = foldlWithKeySV ins (emptyIMSL m) urow0
    ins acc i x = appendIM i (0, x) acc  


luStep :: (Elt a, MonadThrow m, Epsilon a) =>
     SpMatrix a
     -> StateT (IM.IntMap (SList a), IM.IntMap (SList a), IM.Key) m ()
luStep amat = do
  (lmat, umat, t) <- get
  let (umat', utt) = uStep amat lmat umat t
  when (nearZero utt) $
       throwM (NeedsPivoting "LU" (unwords ["U", show (t,t)]) :: MatrixException Double)
  let lmat' = lStep amat lmat umat' utt t
  put (lmat', umat', t + 1)



uStep :: (Elt a, Epsilon a) =>
     SpMatrix a
     -> IM.IntMap (SList a)
     -> IM.IntMap (SList a)
     -> IM.Key
     -> (IM.IntMap (SList a), a)   -- ^ updated U, i'th diagonal element Uii
uStep amat lmat umat i = (umat', udiag) where
  n = ncols amat
  udiag = amat@@!(i,i) - (li <.> umat ! i)
  li = lmat ! i
  umat' = foldr ins umat [i .. n-1]
  ins j acc
      | i == j   = appendIM j (i, udiag) acc
      | isNz uij = appendIM j (i, uij) acc
      | otherwise = acc where
    uij = aij - li <.> uj 
    aij = amat @@! (i,j)
    uj = umat ! j
  

lStep :: (Elt a, Epsilon a) =>
     SpMatrix a
     -> IM.IntMap (SList a)
     -> IM.IntMap (SList a)
     -> a                   -- ^ diagonal element of U (must be nonzero)
     -> IM.Key
     -> IM.IntMap (SList a) -- ^ updated L
lStep amat lmat umat udiag j = foldr ins lmat [j .. m-1] where
  m = nrows amat
  uj = umat ! j
  ins i acc
    | i == j   = appendIM i (j, 1) acc  -- write 1 on the diagonal 
    | isNz lij = appendIM i (j, lij) acc
    | otherwise = acc where
    lij = (aij - li <.> uj)/udiag
    aij = amat @@! (i,j)
    li = lmat ! i






