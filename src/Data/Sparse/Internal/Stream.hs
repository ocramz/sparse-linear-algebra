{-# language ExistentialQuantification #-}
module Data.Sparse.Internal.Stream where

import Data.Vector.Fusion.Stream.Monadic (Stream(..), Step(..))
import Data.Vector.Fusion.Bundle.Size (Size(..))

-- from Data.Vector.Fusion.Stream.Monadic

-- -- | Result of taking a single step in a stream
-- data Step s a = Yield a s  -- ^ a new element and a new seed
--               | Skip    s  -- ^ just a new seed
--               | Done       -- ^ end of stream

-- -- | Result of taking a single step in a stream
-- data Step s a where
--   Yield :: a -> s -> Step s a
--   Skip  :: s -> Step s a
--   Done  :: Step s a

-- instance Functor (Step s) where
--   {-# INLINE fmap #-}
--   fmap f (Yield x s) = Yield (f x) s
--   fmap _ (Skip s) = Skip s
--   fmap _ Done = Done

-- -- | Monadic streams
-- data Stream m a = forall s. Stream (s -> m (Step s a)) s


-- . -- 

data MergeState sa sb i a
  = MergeL sa sb i a
  | MergeR sa sb i a
  | MergeLeftEnded sb
  | MergeRightEnded sa
  | MergeStart sa sb

mergeStreamsWith0 f (Stream stepa sa0) (Stream stepb sb0) =
  Stream step (MergeStart sa0 sb0) -- (toMax na + toMax nb)
  where
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
