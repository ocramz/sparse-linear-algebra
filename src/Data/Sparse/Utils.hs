module Data.Sparse.Utils where

-- import Control.Arrow (first, second)
import Data.Ord
import qualified Data.Vector as V

-- * Misc. utilities






-- | Intersection and union of sparse lists having indices in _ascending_ order
intersectWith :: Ord i => (a -> a -> b) -> [(i, a)] -> [(i, a)] -> [b]
intersectWith f = intersectWith0 (comparing fst) (lift2snd f)

unionWith :: Ord i =>
     (a -> a -> a) -> a -> [(i, a)] -> [(i, a)] -> [(i, a)]
unionWith = unionWith0 compare



-- | Intersection of sparse lists having indices in _descending_ order
intersectWithD :: Ord i => (a -> a -> b) -> [(i, a)] -> [(i, a)] -> [b]
intersectWithD f = intersectWith0 (comparing (Down . fst)) (lift2snd f)

-- | Union of sparse lists having indices in _descending_ order
unionWithD :: Ord i =>
     (a -> a -> a) -> a -> [(i, a)] -> [(i, a)] -> [(i, a)]
unionWithD = unionWith0 (comparing Down)


--
  
intersectWith0 :: (a -> b -> Ordering) -> (a -> b -> c) -> [a] -> [b] -> [c]
intersectWith0 q f = go [] where
  go acc ls@(x : xs) rs@(y : ys) =
    case q x y of EQ -> go (f x y : acc) xs ys
                  LT -> go acc xs rs
                  _  -> go acc ls ys
  go acc [] _ = acc
  go acc _ [] = acc



-- | Lift a binary function onto the second entry of a tuple
lift2snd :: (t -> t1 -> t2) -> (a, t) -> (a1, t1) -> t2
lift2snd f a b = f (snd a) (snd b)



unionWith0 :: (i -> i -> Ordering) -> (a -> a -> a) -> a -> [(i, a)] -> [(i, a)] -> [(i, a)]
unionWith0 q f z = go [] where
  go acc ls@((ix, x) : xs) rs@((iy, y) : ys) =
    case q ix iy of EQ -> go ((ix, f x y) : acc) xs ys
                    LT -> go ((ix, f x z) : acc) xs rs
                    _  -> go ((iy, f z y) : acc) ls ys
  go acc [] r = acc ++ r
  go acc l [] = acc ++ l


-- unionWith0 :: (a -> a -> Ordering) -> (a -> a -> a) -> a -> [a] -> [a] -> [a]
-- unionWith0 q f z = go [] where
--   go acc ls@(x : xs) rs@(y : ys) =
--     case q x y of EQ -> go (f x y : acc) xs ys
--                   LT -> go (f x z : acc) xs rs
--                   _  -> go (f z y : acc) ls ys
--   go acc [] r = acc ++ r
--   go acc l [] = acc ++ l

-- union :: Ord a => [a] -> [a] -> [a]
-- union u_ v_ = go u_ v_ where
--   go [] x = x
--   go y [] = y
--   go uu@(u:us) vv@(v:vs)
--     | u == v =    u : go us vs
--     | u < v =     u : go us vv 
--     | otherwise = v : go uu vs


  






-- | Wrap a function with a null check, returning in Maybe
safe :: (t -> Bool) -> (t -> a) -> t -> Maybe a
safe q f v | q v = Nothing
              | otherwise = Just $ f v


-- | Componentwise tuple operations
-- TODO : use semilattice properties instead
maxTup, minTup :: Ord t => (t, t) -> (t, t) -> (t, t)
maxTup (x1,y1) (x2,y2) = (max x1 x2, max y1 y2)
minTup (x1,y1) (x2,y2) = (min x1 x2, min y1 y2)




-- | integer-indexed ziplist
indexed :: [b] -> [(Int, b)]
indexed xs = indexed' (length xs) xs

indexed' :: Int -> [a] -> [(Int, a)]
indexed' n xs = zip [0 .. n-1] xs




  

-- | ", 2d arrays
indexed2 :: Int -> [c] -> [(Int, Int, c)]
indexed2 m xs = zip3 (concat $ replicate n ii_) jj_ xs where
  ii_ = [0 .. m-1]
  jj_ = concatMap (replicate m) [0 .. n-1]
  ln = length xs
  n = ln `div` m




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
head' = safe V.null V.head

tail' :: V.Vector a -> Maybe (V.Vector a)
tail' = safe V.null V.tail





