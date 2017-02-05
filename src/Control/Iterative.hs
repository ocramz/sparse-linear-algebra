{-# language FlexibleContexts, GeneralizedNewtypeDeriving, DeriveFunctor #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Control.Iterative
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  GPL-style (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- Combinators and helper functions for iterative algorithms, with support for monitoring and exceptions.
--
-----------------------------------------------------------------------------
module Control.Iterative where

import Control.Exception.Common
import Numeric.LinearAlgebra.Class
import Numeric.Eps

import Control.Monad.Catch
import Data.Typeable

import Control.Monad (unless)
import Control.Monad.State.Strict
import Control.Monad.Trans.Writer.CPS
import Control.Monad.Trans.Class (lift)
import qualified Control.Monad.Trans.State.Strict  as MTS -- (runStateT)
-- import Control.Monad.Trans.Writer (runWriterT)

import Data.VectorSpace


-- * Control primitives for bounded iteration with convergence check

-- | transform state until a condition is met
modifyUntil :: MonadState s m => (s -> Bool) -> (s -> s) -> m s
modifyUntil q f = do
  x <- get
  let y = f x
  put y
  if q y then return y
         else modifyUntil q f
    
modifyUntilM :: MonadState s m => (s -> Bool) -> (s -> m s) -> m s
modifyUntilM q f = do
  x <- get
  y <- f x
  put y
  if q y then return y
         else modifyUntilM q f   



-- -- | iterate until convergence is verified or we run out of a fixed iteration budget
-- -- untilConverged :: (MonadState s m, Epsilon a) => (s -> SpVector a) -> (s -> s) -> m s
-- untilConverged :: (MonadState s m, Normed v, Epsilon (Magnitude v)) =>
--      (s -> v) -> (s -> s) -> m s
-- untilConverged fproj = modifyInspectN 200 (normDiffConverged fproj)



-- -- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached (i.e. same thing as `loopUntilAcc` but this one runs in the State monad)
-- modifyInspectN ::
--   MonadState s m =>
--     Int ->           -- iteration budget
--     ([s] -> Bool) -> -- convergence criterion
--     (s -> s) ->      -- state stepping function
--     m s
-- modifyInspectN nitermax q f 
--   | nitermax > 0 = go 0 []
--   | otherwise = error "modifyInspectN : n must be > 0" where
--       go i ll = do
--         x <- get
--         let y = f x
--         if length ll < 2
--           then do put y
--                   go (i + 1) (y : ll)
--           else if q ll || i == nitermax
--                then do put y
--                        return y
--                else do put y
--                        go (i + 1) (take 2 $ y : ll)







-- | `untilConvergedG0` is a special case of `untilConvergedG` that assesses convergence based on the L2 distance to a known solution `xKnown`
untilConvergedG0 ::
  (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable s, Show (Magnitude v), Show s, Ord (Magnitude v)) =>
        String 
     -> Int
     -> (s -> v)
     -> v        -- ^ Known solution
     -> (s -> s)
     -> s
     -> m s
untilConvergedG0 fname nitermax fproj xKnown =
  untilConvergedG fname nitermax fproj (\s -> nearZero (norm2 $ fproj s ^-^ xKnown))
  


-- | This function makes some default choices on the `modifyInspectGuarded` machinery: convergence is assessed using the squared L2 distance between consecutive states, and divergence is detected when this function is increasing between pairs of measurements.
untilConvergedG :: (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable s,
      Show (Magnitude v), Show s, Ord (Magnitude v)) =>
        String
     -> Int 
     -> (s -> v) 
     -> (s -> Bool) 
     -> (s -> s)               -- ^ state update _function_
     -> s 
     -> m s
untilConvergedG fname nitermax fproj qfinal =
  modifyInspectGuarded fname nitermax (convergf fproj) nearZero qdiverg qfinal

-- | ", monadic version
untilConvergedGM ::
  (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable s, Show s) =>
     String
     -> Int
     -> (s -> v)
     -> (s -> Bool)
     -> (s -> m s)  -- ^ state update _arrow_
     -> s
     -> m s
untilConvergedGM fname nitermax fproj qfinal =
  modifyInspectGuardedM fname nitermax (convergf fproj) nearZero qdiverg qfinal     

qdiverg :: Ord a => a -> a -> Bool
qdiverg = (>)

convergf :: Normed v => (t -> v) -> [t] -> Magnitude v
convergf fp [s1, s0] = norm2 (fp s1 ^-^ fp s0)
convergf _ _ = 1/0



-- | `modifyInspectGuarded` is a high-order abstraction of a numerical iterative process. It accumulates a rolling window of 3 states and compares a summary `q` of the latest 2 with that of the previous two in order to assess divergence (e.g. if `q latest2 > q prev2` then it). The process ends when either we hit an iteration budget or relative convergence is verified. The function then assesses the final state with a predicate `qfinal` (e.g. against a known solution; if this is not known, the user can just supply `const True`)
modifyInspectGuarded ::
  (MonadThrow m, Typeable s, Typeable a, Show s, Show a) =>
        String             -- ^ Calling function name
     -> Int                -- ^ Iteration budget
     -> ([s] -> a)         -- ^ State array projection
     -> (a -> Bool)        -- ^ Convergence criterion
     -> (a -> a -> Bool)   -- ^ Divergence criterion
     -> (s -> Bool)        -- ^ Final state evaluation
     -> (s -> s)           -- ^ State evolution
     -> s                  -- ^ Initial state
     -> m s                -- ^ Final state
modifyInspectGuarded fname nits q qc qd qf f x0 =
  modifyInspectGuardedM fname nits q qc qd qf (pure . f) x0

  


-- | ", monadic version
modifyInspectGuardedM
  :: (MonadThrow m, Typeable s, Typeable a, Show s, Show a) =>
     String
     -> Int
     -> ([s] -> a)
     -> (a -> Bool)
     -> (a -> a -> Bool)
     -> (s -> Bool)
     -> (s -> m s)
     -> s
     -> m s
modifyInspectGuardedM fname nitermax q qconverg qdiverg qfinal f x0
  | nitermax > 0 = checkFinal 
  | otherwise = throwM (NonNegError fname nitermax)
  where
    lwindow = 3
    checkFinal = do
      xfinal <- MTS.execStateT (go 0 []) x0
      if qfinal xfinal
        then return xfinal
        else throwM (NotConverged fname nitermax xfinal)
    go i ll = do
      x <- MTS.get
      y <- lift $ f x
      if length ll < lwindow
      then do MTS.put y
              go (i + 1) (y : ll) -- accumulate a rolling state window to observe
      else do
         let qi = q (init ll)     -- summary of latest 2 states
             qt = q (tail ll)     -- "       "  previous 2 states
         if qconverg qi           -- relative convergence  
         then do MTS.put y
                 return ()
         else if i == nitermax    -- end of iterations w/o convergence
              then do
                MTS.put y
                throwM (NotConverged fname nitermax y)
              else
                if qdiverg qi qt  -- diverging
                then throwM (Diverging fname i qi qt)
                else do MTS.put y -- not diverging, keep iterating
                        let ll' = take lwindow $ y : ll
                        go (i + 1) ll' 
             

          
       

          





instance MonadThrow m => MonadThrow (WriterT w m) where
  throwM = lift . throwM


 

-- | iter1 prints output at every iteration until the loop terminates OR is interrupted by an exception, whichever happens first
iter1 :: (MonadThrow m, MonadIO m, Typeable t, Show t) =>
     (t -> m t) -> (t -> String) -> (t -> Bool) -> (t -> Bool) -> t -> m t
iter1 f wf qe qx x0 = execStateT (go 0) x0 where
 go i = do
   x <- get
   y <- lift $ f x
   liftIO $ putStrLn $ wf y
   when (qx y) $ throwM (NotConverged "bla" (i+1)  y)
   put y
   unless (qe y) $ go (i + 1) 


-- | iter2 concatenates output with WriterT but does NOT show
iter2 :: (MonadThrow m, Monoid w, Typeable t, Show t) => (t -> m t)
     -> (t -> w) -> (t -> Bool) -> (t -> Bool) -> t -> m (t, w)
iter2 f wf qe qx x0 = runWriterT $ execStateT (go 0) x0 where
 go i = do
   x <- get
   y <- lift . lift $ f x
   lift $ tell $ wf y
   when (qx y) $ throwM (NotConverged "bla" (i+1)  y)
   put y
   unless (qe y) $ go (i + 1) 



-- test :: IO (Int, String)
test :: IO Int
-- test :: IO ()
test = do
  -- (yt, w ) <- iter2 f wf qe qexc x0
  -- putStrLn w
  -- return yt
  iter1 f wf qe qexc x0
  -- iter2 f wf qe qexc x0
  where
    f = pure . (+ 1)
    wf v = unwords ["state =", show v,"\n"]
    qe = (== 5)
    qexc = (== 6)
    x0 = 0 :: Int

    



-- | Helpers


-- meanl :: (Foldable t, Fractional a) => t a -> a
-- meanl xx = 1/fromIntegral (length xx) * sum xx

-- norm2l :: (Foldable t, Functor t, Floating a) => t a -> a
-- norm2l xx = sqrt $ sum (fmap (**2) xx)

-- | Squared difference of a 2-element list.
-- | NB: unsafe !
diffSqL :: Floating a => [a] -> a
diffSqL xx = (x1 - x2)**2 where [x1, x2] = [head xx, xx!!1]


-- | Relative tolerance :
-- relTol a b := ||a - b|| / (1 + min (||norm2 a||, ||norm2 b||))
relTol :: Normed v => v -> v -> Magnitude v
relTol a b = norm2 (a ^-^ b) / m where
  m = 1 + min (norm2 a) (norm2 b)


-- | convergence test
normDiff :: Normed v => (t -> v) -> [t] -> Magnitude v
normDiff fp [x1, x0] = norm2Sq (fp x0 ^-^ fp x1)
normDiff _ _ = 1/0   --- ugggggh

normDiffConverged :: Normed v => (s -> v) -> [s] -> Bool
normDiffConverged fp ll = nearZero (normDiff fp ll)
