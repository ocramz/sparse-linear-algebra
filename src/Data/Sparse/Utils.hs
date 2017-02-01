module Data.Sparse.Utils where

import qualified Data.Vector as V

-- * Misc. utilities


-- | Wrap a function with a null check, returning in Maybe
harness :: (t -> Bool) -> (t -> a) -> t -> Maybe a
harness q f v | q v = Nothing
              | otherwise = Just $ f v


-- | Componentwise tuple operations
-- TODO : use semilattice properties instead
maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)




-- | integer-indexed ziplist
denseIxArray :: [b] -> [(Int, b)]
denseIxArray xs = zip [0..length xs-1] xs 

-- | ", 2d arrays
denseIxArray2 :: Int -> [c] -> [(Int, Int, c)]
denseIxArray2 m xs = zip3 (concat $ replicate n ii_) jj_ xs where
  ii_ = [0 .. m-1]
  jj_ = concatMap (replicate m) [0 .. n-1]
  ln = length xs
  n = ln `div` m


-- indexed :: [b] -> [(Int, b)]
-- indexed x = zip [0 .. length x - 1] x



-- folds

-- | foldr over the results of a fmap
foldrMap :: (Foldable t, Functor t) => (a -> b) -> (b -> c -> c) -> c -> t a -> c
foldrMap ff gg x0 = foldr gg x0 . fmap ff

-- | strict left fold
foldlStrict :: (a -> b -> a) -> a -> [b] -> a
foldlStrict f = go
  where
    go z []     = z
    go z (x:xs) = let z' = f z x in z' `seq` go z' xs

-- | indexed right fold
ifoldr :: Num i =>
     (a -> b -> b) -> b -> (i -> c -> d -> a) -> c -> [d] -> b  
ifoldr mjoin mneutral f  = go 0 where
  go i z (x:xs) = mjoin (f i z x) (go (i+1) z xs)
  go _ _ [] = mneutral


-- ** Bounds checking
type LB = Int
type UB = Int

inBounds :: LB -> UB -> Int -> Bool
inBounds ibl ibu i = i>= ibl && i<ibu

inBounds2 :: (LB, UB) -> (Int, Int) -> Bool
inBounds2 (ibl,ibu) (ix,iy) = inBounds ibl ibu ix && inBounds ibl ibu iy


-- ", lower bound = 0
inBounds0 :: UB -> Int -> Bool
inBounds0 = inBounds 0

inBounds02 :: (UB, UB) -> (Int, Int) -> Bool
inBounds02 (bx,by) (i,j) = inBounds0 bx i && inBounds0 by j




-- ** Safe indexing



head' :: V.Vector a -> Maybe a
head' = harness V.null V.head

tail' :: V.Vector a -> Maybe (V.Vector a)
tail' = harness V.null V.tail




-- | a cons-based moving-window datatype of length at least 3

data W3 a = W3 a a a [a] deriving (Eq, Show)
initW3 :: a -> a -> a -> [a] -> W3 a
initW3 = W3

pushW3 :: a -> W3 a -> W3 a
pushW3 i (W3 i0 i1 i2 is) = W3 i i0 i1 (i2 : initSafe is)

fstW3 (W3 i _ _ _) = i
sndW3 (W3 _ i _ _) = i
thirdW3 (W3 _ _ i _) = i

withInitW3, withTailW3 :: W3 t -> (t -> t -> t1) -> t1
withInitW3 (W3 a b _ _) f = f a b
withTailW3 (W3 _ b c _) f = f b c

initSafe (x:xs) = x : initSafe xs
initSafe [] = []
