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

import Control.Monad (when)
-- import Control.Monad.Trans.Reader
import Control.Monad.State.Strict
-- import Control.Monad.Trans.Writer.CPS
import Control.Monad.Trans.Class (lift)
import qualified Control.Monad.Trans.State.Strict  as MTS -- (runStateT)

import Data.VectorSpace


data ConvergenceStatus a = BufferNotReady
                         | Converging
                         | Converged a
                         | Diverging a a
                         | NotConverged deriving (Eq, Show)

     
data IterationConfig a b =
  IterConf { numIterationsMax :: Int,
             printDebugInfo :: Bool,
             iterationView :: a -> b, 
             printDebugIO :: b -> IO ()}
instance Show (IterationConfig a b) where
  show (IterConf n qd _ _) = unwords ["Max. # of iterations:",show n,", print debug information:", show qd]
  


-- * Control primitives for bounded iteration with convergence check

-- -- | transform state until a condition is met
modifyUntil :: MonadState s m => (s -> Bool) -> (s -> s) -> m s
modifyUntil q f = modifyUntilM q (pure . f)
       
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

-- execState (modifyInspectN 2 (nearZero . diffSqL) (/2)) (1 :: Double) `shouldBe` 1/8


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

data S = S {unS1 :: Double, unS2 :: String} deriving (Eq, Show)
liftS1 f (S x i) = S (f x) i
s0 = S 1 "blah"
ic1 = IterConf 2 True unS1 print





-- | `untilConvergedG0` is a special case of `untilConvergedG` that assesses convergence based on the L2 distance to a known solution `xKnown`
untilConvergedG0 ::
  (Normed v, MonadThrow m, MonadIO m, Typeable (Magnitude v), Typeable s, Show s) =>
     String
     -> IterationConfig s v
     -> v                    -- ^ Known value
     -> (s -> s)
     -> s
     -> m s
untilConvergedG0 fname config xKnown =
  modifyInspectGuarded fname config norm2Diff nearZero qdiverg qfin where
    qfin s = nearZero $ norm2 (xKnown ^-^ s)
  


-- | This function makes some default choices on the `modifyInspectGuarded` machinery: convergence is assessed using the squared L2 distance between consecutive states, and divergence is detected when this function is increasing between pairs of measurements.
untilConvergedG :: (Normed v, MonadThrow m, MonadIO m, Typeable (Magnitude v), Typeable s, Show s) =>
        String
     -> IterationConfig s v
     -> (v -> Bool)
     -> (s -> s)               
     -> s 
     -> m s
untilConvergedG fname config =
  modifyInspectGuarded fname config norm2Diff nearZero qdiverg


-- | ", monadic version
untilConvergedGM ::
  (Normed v, MonadThrow m, MonadIO m, Typeable (Magnitude v), Typeable s, Show s) =>
     String
     -> IterationConfig s v
     -> (v -> Bool)
     -> (s -> m s)          
     -> s
     -> m s
untilConvergedGM fname config =
  modifyInspectGuardedM fname config norm2Diff nearZero qdiverg





-- | `modifyInspectGuarded` is a high-order abstraction of a numerical iterative process. It accumulates a rolling window of 3 states and compares a summary `q` of the latest 2 with that of the previous two in order to assess divergence (e.g. if `q latest2 > q prev2` then it). The process ends when either we hit an iteration budget or relative convergence is verified. The function then assesses the final state with a predicate `qfinal` (e.g. against a known solution; if this is not known, the user can just supply `const True`)
modifyInspectGuarded ::
  (MonadThrow m, MonadIO m, Typeable s, Typeable a, Show s, Show a) =>
        String              -- ^ Calling function name
     -> IterationConfig s v -- ^ Configuration
     -> ([v] -> a)          -- ^ State summary array projection
     -> (a -> Bool)         -- ^ Convergence criterion
     -> (a -> a -> Bool)    -- ^ Divergence criterion
     -> (v -> Bool)         -- ^ Final state acceptance criterion
     -> (s -> s)            -- ^ State evolution
     -> s                   -- ^ Initial state
     -> m s                 -- ^ Final state
modifyInspectGuarded fname config q qc qd qfin f x0 =
  modifyInspectGuardedM fname config q qc qd qfin (pure . f) x0

  


-- | ", monadic version
modifyInspectGuardedM ::
  (MonadThrow m, MonadIO m, Typeable s, Show s, Typeable a, Show a) =>
     String
     -> IterationConfig s v
     -> ([v] -> a)
     -> (a -> Bool)
     -> (a -> a -> Bool)
     -> (v -> Bool)
     -> (s -> m s)
     -> s
     -> m s
modifyInspectGuardedM fname config sf qconverg qdiverg qfinal f x0 
  | nitermax > 0 = xout
  | otherwise = throwM (NonNegError fname nitermax)
  where
    lwindow = 3
    nitermax = numIterationsMax config
    pf = iterationView config
    xout = do
      x <- MTS.execStateT (go 0 []) x0
      if qfinal (pf x) then return x
                       else throwM (NotConvergedE fname nitermax x)
    checkConvergStatus i ll
      | i == nitermax - 1 && not (qconverg qi) = NotConverged    
      | length ll < lwindow = BufferNotReady
      | qdiverg qi qt && not (qconverg qi) = Diverging qi qt        
      | qconverg qi = Converged qi
      | otherwise = Converging
      where llf = pf <$> ll
            qi = sf $ init llf         -- summary of latest 2 states
            qt = sf $ tail llf         -- "       "  previous 2 states
    go i ll = do
      x <- MTS.get
      y <- lift $ f x
      when (printDebugInfo config) $ liftIO $ printDebugIO config (pf y) 
      case checkConvergStatus i ll of
        BufferNotReady -> do  
          MTS.put y
          let ll' = y : ll    -- cons current state to buffer
          go (i + 1) ll'
        Converged qi -> MTS.put y
        Diverging qi qt -> do
          MTS.put y
          throwM (DivergingE fname i qi qt)
        Converging -> do
          MTS.put y
          let ll' = take lwindow (y : ll) -- rolling state window
          go (i + 1) ll'
        NotConverged -> do
          MTS.put y
          throwM (NotConvergedE fname nitermax y)
             






    



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


qdiverg :: Ord a => a -> a -> Bool
qdiverg = (>)

norm2Diff [s1, s0] = norm2 (s1 ^-^ s0)
norm2Diff _ = 1/0


-- norm2Diff :: Normed v => (t -> v) -> [t] -> Magnitude v
-- norm2Diff fp [s1, s0] = norm2 (fp s1 ^-^ fp s0)
-- norm2Diff _ _ = 1/0







-- playground

-- instance MonadThrow m => MonadThrow (WriterT w m) where
--   throwM = lift . throwM

-- -- | iter0 also accepts a configuration, e.g. for optional printing of debug info
-- -- iter0 :: MonadIO m =>
-- --      Int -> (s -> m s) -> (s -> String) -> IterationConfig s -> s -> m s
-- iter0 nmax f sf config x0 = flip runReaderT config $ MTS.execStateT (go (0 :: Int)) x0
--  where
--   go i = do
--     x <- get
--     c <- lift $ asks printDebugInfo  -- neat
--     y <- lift . lift $ f x           -- not neat
--     when c $ liftIO $ putStrLn $ sf y 
--     put y
--     unless (i >= nmax) (go $ i + 1)
 

-- -- | iter1 prints output at every iteration until the loop terminates OR is interrupted by an exception, whichever happens first
-- -- iter1 :: (MonadThrow m, MonadIO m, Typeable t, Show t) =>
-- --      (t -> m t) -> (t -> String) -> (t -> Bool) -> (t -> Bool) -> t -> m t
-- iter1 f wf qe qx x0 = execStateT (go 0) x0 where
--  go i = do
--    x <- get
--    y <- lift $ f x
--    _ <- liftIO $ wf y
--    when (qx y) $ throwM (NotConvergedE "bla" (i+1)  y)
--    put y
--    unless (qe y) $ go (i + 1) 


-- -- | iter2 concatenates output with WriterT but does NOT `tell` any output if an exception is raised before the end of the loop
-- iter2 :: (MonadThrow m, Monoid w, Typeable t, Show t) => (t -> m t)
--      -> (t -> w) -> (t -> Bool) -> (t -> Bool) -> t -> m (t, w)
-- iter2 f wf qe qx x0 = runWriterT $ execStateT (go 0) x0 where
--  go i = do
--    x <- get
--    y <- lift . lift $ f x
--    lift $ tell $ wf y
--    when (qx y) $ throwM (NotConvergedE "bla" (i+1)  y)
--    put y
--    unless (qe y) $ go (i + 1) 


-- -- test :: IO (Int, [String])
-- test :: IO Int
-- -- test :: IO ()
-- test = do
--   (yt, w ) <- iter2 f wf qe qexc x0
--   putStrLn w
--   return yt
--   -- iter1 f wf qe qexc x0
--   where
--     f = pure . (+ 1)
--     wf v = unwords ["state =", show v]
--     qe = (== 5)
--     qexc = (== 3)
--     x0 = 0 :: Int
