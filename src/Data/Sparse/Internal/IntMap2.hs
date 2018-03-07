module Data.Sparse.Internal.IntMap2 where

import qualified Data.Sparse.Internal.IntM as I

-- import Numeric.LinearAlgebra.Class
import qualified Data.IntMap.Strict as IM
import Data.Sparse.Types

import Data.Maybe
import GHC.Exts






-- * Insertion

-- | Insert an element
-- insertIM2 ::
--   IM.Key -> IM.Key -> a -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
insertIM2
  :: IM.Key -> IM.Key -> a -> I.IntM (I.IntM a) -> I.IntM (I.IntM a)
insertIM2 i j x imm = I.insert i ro imm where
  ro = maybe (I.singleton j x) (I.insert j x) (I.lookup i imm)
{-# inline insertIM2 #-}  

-- * Lookup

-- |Lookup a key
-- lookupIM2 ::
--   IM.Key -> IM.Key -> IM.IntMap (IM.IntMap a) -> Maybe a
lookupIM2 i j imm = I.lookup i imm >>= I.lookup j
{-# inline lookupIM2 #-}  

-- | Lookup with default 0
lookupWD_IM :: Num a => IM.IntMap (IM.IntMap a) -> (IxRow, IxCol) -> a
lookupWD_IM im (i,j) = fromMaybe 0 (IM.lookup i im >>= IM.lookup j)




-- |Populate an IM2 from a list of (row index, column index, value)  
-- fromListIM2 ::
--   Foldable t =>
--      t (IM.Key, IM.Key, a) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
fromListIM2 iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertIM2 i j x t


-- * Folding

-- |Indexed left fold over an IM2, with general accumulator
-- ifoldlIM2' :: (IM.Key -> IM.Key -> a -> b -> b) -> b -> IM.IntMap (IM.IntMap a) -> b
ifoldlIM2' f empty mm = I.foldlWithKey' accRow empty mm where
  accRow acc i r = I.foldlWithKey' (accElem i) acc r
  accElem i acc j x = f i j x acc
{-# inline ifoldlIM2' #-}

-- |Indexed left fold over an IM2
-- ifoldlIM2 ::
--   (IM.Key -> IM.Key -> t -> IM.IntMap a -> IM.IntMap a) ->
--   IM.IntMap (IM.IntMap t) ->  
--   IM.IntMap a
ifoldlIM2 f m         = I.foldlWithKey' accRow I.empty m where
  accRow    acc i row = I.foldlWithKey' (accElem i) acc row
  accElem i acc j x   = f i j x acc
{-# inline ifoldlIM2 #-}  

-- |Left fold over an IM2, with general accumulator
-- foldlIM2 :: (a -> b -> b) -> b -> IM.IntMap (IM.IntMap a) -> b
foldlIM2 f empty mm = foldl accRow empty mm where
  accRow acc r = foldl accElem acc r
  accElem acc x = f x acc
{-# inline foldlIM2 #-}


-- | Inner indices become outer ones and vice versa. No loss of information because both inner and outer IntMaps are nubbed.
-- transposeIM2 :: IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
transposeIM2 = ifoldlIM2 (flip insertIM2)
{-# inline transposeIM2 #-}

-- specialized folds

-- -- extract diagonal elements
-- extractDiagonalIM2 :: IM.IntMap (IM.IntMap a) -> [a]
-- extractDiagonalIM2 = ifoldlIM2' (\i j x xs -> if i==j then x : xs else xs) []




-- * Filtering

-- |Map over outer IM and filter all inner IM's
-- ifilterIM2 ::
--   (IM.Key -> IM.Key -> a -> Bool) ->
--   IM.IntMap (IM.IntMap a) ->
--   IM.IntMap (IM.IntMap a)
ifilterIM2 f  =
  I.mapWithKey (\irow row -> I.filterWithKey (f irow) row)
{-# inline ifilterIM2 #-}

-- |Specialized filtering : keep only sub-diagonal elements
-- filterSubdiag :: IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
filterSubdiag = ifilterIM2 (\i j _ -> i>j)

-- countSubdiagonalNZ :: IM.IntMap (IM.IntMap a) -> Int
countSubdiagonalNZ im =
  I.size $ I.filterI (not . null) (filterSubdiag im)

-- |List of (row, col) indices of (nonzero) subdiagonal elements
-- subdiagIndices :: IM.IntMap (IM.IntMap a) -> [(IM.Key, IM.Key)]
subdiagIndices im = concatMap rpairs $ toList (I.keys <$> im') where
  im' = filterSubdiag im

rpairs :: (a, [b]) -> [(a, b)]
rpairs (i, jj@(_:_)) = zip (replicate (length jj) i) jj
rpairs (_, []) = []

-- -- list of (row, col) indices of elements that satisfy a criterion
-- indicesThatIM2 ::
--   (IM.Key -> IM.IntMap a -> Bool) -> IM.IntMap (IM.IntMap a) -> [(IM.Key, IM.Key)]
-- indicesThatIM2 f im = concatMap rpairs $ IM.toList (IM.map IM.keys im') where
--   im' = IM.filterWithKey f im

  


-- * Mapping

-- |Map over IM2
mapIM2 :: (a -> b) -> I.IntM (I.IntM a) -> I.IntM (I.IntM b)
mapIM2 = fmap . fmap





-- |Indexed map over IM2
-- imapIM2 ::
--   (IM.Key -> IM.Key -> a -> b) ->
--   IM.IntMap (IM.IntMap a) ->
--   IM.IntMap (IM.IntMap b)
imapIM2 f im = I.mapWithKey ff im where
  ff j x = I.mapWithKey (`f` j) x



-- |Mapping keys
-- mapKeysIM2 ::
--   (IM.Key -> IM.Key) -> (IM.Key -> IM.Key) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
mapKeysIM2 fi fj im = adjCols <$> adjRows where
  adjRows = I.mapKeys fi im
  adjCols = I.mapKeys fj




-- map over a single `column`

-- mapColumnIM2 :: (b -> b) -> IM.IntMap (IM.IntMap b) -> Int -> IM.IntMap (IM.IntMap b)
mapColumnIM2 f im jj = imapIM2 (\i j x -> if j == jj then f x else x) im




-- | utilities




