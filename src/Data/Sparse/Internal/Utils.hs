{-# language ExistentialQuantification, TypeFamilies, FlexibleInstances, MultiParamTypeClasses #-}
module Data.Sparse.Internal.Utils where

import Control.Monad (unless)
import Control.Monad.State
import Control.Monad.ST

import Data.Ord (comparing)

import qualified Data.Vector as V 
-- import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Mutable as VM
import qualified Data.Vector.Algorithms.Merge as VA

-- import Data.Ix
-- import Data.Maybe

-- import Data.Sparse.Types
-- import Numeric.LinearAlgebra.Class





-- | Given a number of rows(resp. columns) `n` and a _sorted_ Vector of Integers in increasing order (containing the row(col) indices of nonzero entries), return the cumulative vector of nonzero entries of length `n + 1` (the "row(col) pointer" of the CSR(CSC) format). NB: Fused count-and-accumulate
-- E.g.:
-- > csrPtrV (==) 4 (V.fromList [1,1,2,3])
-- [0,0,2,3,4]
csrPtrV :: (a -> Int -> Bool) -> Int -> V.Vector a -> V.Vector Int
csrPtrV eqf n xs = V.create createf where
  createf :: ST s (VM.MVector s Int)
  createf = do
    let c = 0
    vm <- VM.new (n + 1)
    VM.write vm 0 0  -- write `0` at position 0
    let loop v ll i count | i == n = return ()
                          | otherwise = do
                              let lp = V.length $ V.takeWhile (`eqf` i) ll
                                  count' = count + lp
                              VM.write v (i + 1) count'
                              loop v (V.drop lp ll) (succ i) count'
    loop vm xs 0 c
    return vm


-- csrPtrV' eqf n xs = V.create createf where
--   createf :: ST s (VM.MVector s Int)
--   createf = do
--     let c = 0
--     vm <- VM.new (n + 1)
--     VM.write vm 0 0
--     numLoop fw 0 where
--       fw ix = let lp = V.length $ V.takeWhile (`eqf` ix) ll
--                   count
--     return vm

numLoop :: Monad m => (Int -> m a) -> Int -> m ()
numLoop fm n = go 0 where
  go i | i == n = return ()
       | otherwise = do
           _ <- fm i
           go (succ i)

numLoopST' :: Monad m => (Int -> s -> m a) -> Int -> (s -> s) -> s -> m ()
numLoopST' fm n fs s0 = go 0 s0 where
  go i s | i == n = return ()
         | otherwise = do
             _ <- fm i s
             go (succ i) (fs s)

numLoopST'' ::
  MonadState s m => (Int -> s -> m a) -> Int -> (s -> s) -> m ()
numLoopST'' fm n fs = go 0 where
  go i = unless (i == n) $ do
    s <- get
    _ <- fm i s -- ignore result of `fm`
    put $ fs s
    go (succ i)

  
             




  
-- | O(N) : Intersection between sorted vectors, in-place updates
intersectWith :: Ord o =>
    (a -> o) -> (a -> a -> c) -> V.Vector a -> V.Vector a -> V.Vector c
intersectWith f = intersectWithCompare (comparing f)


intersectWithCompare ::
  (a1 -> a2 -> Ordering) ->
  (a1 -> a2 -> a) ->
  V.Vector a1 ->
  V.Vector a2 ->
  V.Vector a
intersectWithCompare fcomp g u_ v_ = V.create $ do
  let n = min (V.length u_) (V.length v_)
  vm <- VM.new n
  let go u_ v_ i vm | V.null u_ || V.null v_ || i == n = return (vm, i)
                    | otherwise =  do
         let (u,us) = (V.head u_, V.tail u_)
             (v,vs) = (V.head v_, V.tail v_)
         case fcomp u v of EQ -> do VM.write vm i (g u v)
                                    go us vs (i + 1) vm
                           LT -> go us v_ i vm
                           GT -> go u_ vs i vm
  (vm', i') <- go u_ v_ 0 vm
  let vm'' = VM.take i' vm'
  return vm''







-- | O(N) : Union between sorted vectors, in-place updates
unionWith :: Ord o =>
  (t -> o) -> (t -> t -> a) -> t -> V.Vector t -> V.Vector t -> V.Vector a
unionWith f = unionWithCompare (comparing f)


unionWithCompare ::
  (t -> t -> Ordering) -> (t -> t -> a) -> t -> V.Vector t -> V.Vector t -> V.Vector a
unionWithCompare fcomp g z u_ v_ = V.create $ do
  let n = (V.length u_) + (V.length v_)
  vm <- VM.new n
  let go u_ v_ i vm
        | (V.null u_ && V.null v_) || i==n = return (vm, i)
        | V.null u_ = do
            VM.write vm i (g z (V.head v_))
            go u_ (V.tail v_) (i+1) vm
        | V.null v_ = do
            VM.write vm i (g (V.head u_) z)
            go (V.tail u_) v_ (i+1) vm
        | otherwise =  do
           let (u,us) = (V.head u_, V.tail u_)
               (v,vs) = (V.head v_, V.tail v_)
           case fcomp u v of EQ -> do VM.write vm i (g u v)
                                      go us vs (i + 1) vm
                             LT -> do VM.write vm i (g u z)
                                      go us v_ (i + 1) vm
                             GT -> do VM.write vm i (g z v)
                                      go u_ vs (i + 1) vm
  (vm', nfin) <- go u_ v_ 0 vm
  let vm'' = VM.take nfin vm'
  return vm''





-- * Sorting
sortWith :: Ord o => (t -> o) -> V.Vector t -> V.Vector t
sortWith f v = V.modify (VA.sortBy (comparing f)) v

sortWith3 :: Ord o =>
   ((a, b, c) -> o) ->
   V.Vector a ->
   V.Vector b ->
   V.Vector c ->
   V.Vector (a, b, c)
sortWith3 f x y z = sortWith f $ V.zip3 x y z


sortByFst3 :: Ord a => V.Vector a -> V.Vector b -> V.Vector c -> V.Vector (a, b, c)
sortByFst3 = sortWith3 fst3

sortBySnd3 :: Ord b => V.Vector a -> V.Vector b -> V.Vector c -> V.Vector (a, b, c)
sortBySnd3 = sortWith3 snd3












-- * Utilities

-- ** 3-tuples
fst3 :: (a, b, c) -> a
fst3 (a, _, _) = a
snd3 :: (a, b, c) -> b
snd3 (_, b, _) = b
third3 :: (a, b, c) -> c
third3 (_, _, c) = c


tail3 :: (t, t1, t2) -> (t1, t2)
tail3 (_,j,x) = (j,x)



mapFst3 :: (a -> b) -> (a, y, z) -> (b, y, z)
mapFst3 f (a, b, c) = (f a, b, c)

mapSnd3 :: (a -> b) -> (x, a, z) -> (x, b, z)
mapSnd3 f (a, b, c) = (a, f b, c)

mapThird3 :: (a -> b) -> (x, y, a) -> (x, y, b)
mapThird3 f (a, b, c) = (a, b, f c)



lift2 :: (a -> b) -> (b -> b -> c) -> a -> a -> c
lift2 p f t1 t2 = f (p t1) (p t2)






-- | Stream fusion based version of the above, from [1]
-- [1] : https://www.schoolofhaskell.com/user/edwardk/revisiting-matrix-multiplication/part-3

data Stream m a = forall s . Stream (s -> m (Step s a)) s 

data Step s a = Yield a s | Skip s | Done  

data MergeState sa sb i a
  = MergeL sa sb i a
  | MergeR sa sb i a
  | MergeLeftEnded sb
  | MergeRightEnded sa
  | MergeStart sa sb

mergeStreamsWith :: (Ord i, Monad m) =>
     (a -> a -> Maybe a) -> Stream m (i, a) -> Stream m (i, a) -> Stream m (i, a)
mergeStreamsWith f (Stream stepa sa0) (Stream stepb sb0)
  = Stream step (MergeStart sa0 sb0)  where
 step (MergeStart sa sb) = do
    r <- stepa sa
    return $ case r of
      Yield (i, a) sa' -> Skip (MergeL sa' sb i a)
      Skip sa'         -> Skip (MergeStart sa' sb)
      Done             -> Skip (MergeLeftEnded sb)
 step (MergeL sa sb i a) = do
    r <- stepb sb
    return $ case r of
      Yield (j, b) sb' -> case compare i j of
        LT -> Yield (i, a)     (MergeR sa sb' j b)
        EQ -> case f a b of
           Just c  -> Yield (i, c) (MergeStart sa sb')
           Nothing -> Skip (MergeStart sa sb')
        GT -> Yield (j, b)     (MergeL sa sb' i a)
      Skip sb' -> Skip (MergeL sa sb' i a)
      Done     -> Yield (i, a) (MergeRightEnded sa)
 step (MergeR sa sb j b) = do
    r <- stepa sa
    return $ case r of
      Yield (i, a) sa' -> case compare i j of
        LT -> Yield (i, a)     (MergeR sa' sb j b)
        EQ -> case f a b of
          Just c  -> Yield (i, c) (MergeStart sa' sb)
          Nothing -> Skip (MergeStart sa' sb)
        GT -> Yield (j, b)     (MergeL sa' sb i a)
      Skip sa' -> Skip (MergeR sa' sb j b)
      Done     -> Yield (j, b) (MergeLeftEnded sb)
 step (MergeLeftEnded sb) = do
    r <- stepb sb
    return $ case r of
      Yield (j, b) sb' -> Yield (j, b) (MergeLeftEnded sb')
      Skip sb'         -> Skip (MergeLeftEnded sb')
      Done             -> Done
 step (MergeRightEnded sa) = do
    r <- stepa sa
    return $ case r of
      Yield (i, a) sa' -> Yield (i, a) (MergeRightEnded sa')
      Skip sa'         -> Skip (MergeRightEnded sa')
      Done             -> Done
 {-# INLINE [0] step #-}
{-# INLINE [1] mergeStreamsWith #-} 


-- test data

-- m0 = V.fromList [O (0,0,1), O(0,1,2)]
-- m1 = V.fromList [O (0,0,1), O(0,2,3)]


isOrderedV :: Ord a => V.Vector a -> Bool
isOrderedV l = V.all (== True) $ V.zipWith (<) l (V.tail l)











