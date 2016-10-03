module Math.Linear.Sparse.IntMap where

import qualified Data.IntMap.Strict as IM




-- | ========= IntMap-of-IntMap (IM2) stuff


-- insert an element
insertIM2 ::
  IM.Key -> IM.Key -> a -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
insertIM2 i j x imm = IM.insert i (IM.insert j x ro) imm where
  ro = maybe (IM.singleton j x) (IM.insert j x) (IM.lookup i imm)

-- lookup a key
lookupIM2 ::
  IM.Key -> IM.Key -> IM.IntMap (IM.IntMap a) -> Maybe a
lookupIM2 i j imm = IM.lookup i imm >>= IM.lookup j

-- populate an IM2 from a list of (row index, column index, value)  
fromListIM2 ::
  Foldable t =>
     t (IM.Key, IM.Key, a) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
fromListIM2 iix sm = foldl ins sm iix where
  ins t (i,j,x) = insertIM2 i j x t


-- | folding

-- indexed fold over an IM2
ifoldlIM2' :: (IM.Key -> IM.Key -> a -> b -> b) -> b -> IM.IntMap (IM.IntMap a) -> b
ifoldlIM2' f empty mm = IM.foldlWithKey' accRow empty mm where
  accRow acc i r = IM.foldlWithKey' (accElem i) acc r
  accElem i acc j x = f i j x acc
{-# inline ifoldlIM2' #-}

ifoldlIM2 ::
  (IM.Key -> IM.Key -> t -> IM.IntMap a -> IM.IntMap a) ->
  IM.IntMap (IM.IntMap t) ->  
  IM.IntMap a
ifoldlIM2 f m         = IM.foldlWithKey' accRow IM.empty m where
  accRow    acc i row = IM.foldlWithKey' (accElem i) acc row
  accElem i acc j x   = f i j x acc
{-# inline ifoldlIM2 #-}  

foldlIM2 :: (a -> b -> b) -> b -> IM.IntMap (IM.IntMap a) -> b
foldlIM2 f empty mm = IM.foldl accRow empty mm where
  accRow acc r = IM.foldl accElem acc r
  accElem acc x = f x acc
{-# inline foldlIM2 #-}


-- transposeIM2 : inner indices become outer ones and vice versa. No loss of information because both inner and outer IntMaps are nubbed.
transposeIM2 :: IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
transposeIM2 = ifoldlIM2 (flip insertIM2)
{-# inline transposeIM2 #-}

-- specialized folds

-- -- extract diagonal elements
-- extractDiagonalIM2 :: IM.IntMap (IM.IntMap a) -> [a]
-- extractDiagonalIM2 = ifoldlIM2' (\i j x xs -> if i==j then x : xs else xs) []




-- | filtering

-- map over outer IM and filter all inner IM's
ifilterIM2 ::
  (IM.Key -> IM.Key -> a -> Bool) ->
  IM.IntMap (IM.IntMap a) ->
  IM.IntMap (IM.IntMap a)
ifilterIM2 f  =
  IM.mapWithKey (\irow row -> IM.filterWithKey (f irow) row) 
{-# inline ifilterIM2 #-}

-- specialized filtering function

-- keep only sub-diagonal elements
filterSubdiag :: IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
filterSubdiag = ifilterIM2 (\i j _ -> i>j)

countSubdiagonalNZ :: IM.IntMap (IM.IntMap a) -> Int
countSubdiagonalNZ im =
  IM.size $ IM.filter (not . IM.null) (filterSubdiag im)

-- list of (row, col) indices of (nonzero) subdiagonal elements
subdiagIndices :: IM.IntMap (IM.IntMap a) -> [(IM.Key, IM.Key)]
subdiagIndices im = concatMap rpairs $ IM.toList (IM.map IM.keys im') where
  im' = filterSubdiag im

rpairs :: (a, [b]) -> [(a, b)]
rpairs (i, jj@(_:_)) = zip (replicate (length jj) i) jj
rpairs (_, []) = []

-- -- list of (row, col) indices of elements that satisfy a criterion
-- indicesThatIM2 ::
--   (IM.Key -> IM.IntMap a -> Bool) -> IM.IntMap (IM.IntMap a) -> [(IM.Key, IM.Key)]
-- indicesThatIM2 f im = concatMap rpairs $ IM.toList (IM.map IM.keys im') where
--   im' = IM.filterWithKey f im

  


-- | mapping

-- map over IM2

mapIM2 :: (a -> b) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap b)
mapIM2 = IM.map . IM.map   -- imapIM2 (\_ _ x -> f x)


-- indexed map over IM2
imapIM2 ::
  (IM.Key -> IM.Key -> a -> b) ->
  IM.IntMap (IM.IntMap a) ->
  IM.IntMap (IM.IntMap b)
imapIM2 f im = IM.mapWithKey ff im where
  ff j x = IM.mapWithKey (`f` j) x



-- mapping keys

mapKeysIM2 ::
  (IM.Key -> IM.Key) -> (IM.Key -> IM.Key) -> IM.IntMap (IM.IntMap a) -> IM.IntMap (IM.IntMap a)
mapKeysIM2 fi fj im = IM.map adjCols adjRows where
  adjRows = IM.mapKeys fi im
  adjCols = IM.mapKeys fj 




-- map over a single `column`

mapColumnIM2 :: (b -> b) -> IM.IntMap (IM.IntMap b) -> Int -> IM.IntMap (IM.IntMap b)
mapColumnIM2 f im jj = imapIM2 (\i j x -> if j == jj then f x else x) im





-- sparsification :

