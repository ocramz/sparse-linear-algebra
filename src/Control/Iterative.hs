{-# language FlexibleContexts, GeneralizedNewtypeDeriving, DeriveFunctor, DeriveGeneric #-}
{-# language OverloadedStrings #-}
{-# language MultiParamTypeClasses #-}
{-# language FlexibleInstances #-}
{-# language ScopedTypeVariables #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Control.Iterative
-- Copyright   :  (c) Marco Zocca 2017-2018
-- License     :  GPL-style (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- Combinators and helper functions for iterative algorithms, with support for monitoring and exceptions.
--
-----------------------------------------------------------------------------
module Control.Iterative (
  -- * Iteration types
  ConvergenceStatus(..), IterConfig(..), ConvergConfig(..), LoopState(..),
  convergenceL2,
  -- * Iteration combinators
  modifyInspectGuardedM, modifyUntil, modifyUntilM, modifyUntilM_,
  modifyUntilM',
  -- * Helpers
  getBuffers, mkLoopState,
  onRangeSparse, onRangeSparseM, unfoldZipM0, unfoldZipM, combx,
  sqDiffPairs, sqDiff, relRes, diffSqL, relTol, norm2Diff
) where

import Control.Applicative

import Control.Monad (when, replicateM)
import Control.Monad.Reader (MonadReader(..), asks)
import Control.Monad.State.Strict (MonadState(..), get, put, gets)
import Control.Monad.Writer.Class (MonadWriter)
import Control.Monad.Writer.Strict (WriterT, runWriterT)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT, execStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Catch (Exception(..), MonadThrow(..), throwM)

-- import Data.Bool (bool)
import Data.Char (toUpper)
import Data.Semigroup
import Data.Monoid (Sum(..), Product(..))

import Data.Typeable
import qualified Control.Exception as E (Exception, Handler)

import Data.Foldable (foldrM)

import GHC.Generics

import Control.Exception.Common
import Control.Iterative.Internal

import Numeric.LinearAlgebra.Class
import Numeric.Eps



-- * ITERATION

-- | Iteration state machine
data ConvergenceStatus s a =
  BufferNotReady
  | Converging
  | Converged s -- ^ Final state
  | Diverging a a 
  | NotConverged s  -- ^ Final state
  deriving (Eq, Show)

-- | Configuration data for the iterative process
data IterConfig s t = IterConfig {
    icFunctionName :: String -- ^ Name of calling function, for logging purposes
  , icNumIterationsMax :: Int -- ^ Max # of iterations
  , icStateWindowLength :: Int -- ^ # of states used to assess convergence/divergence
  , icStateProj :: s -> t     -- ^ Project the state
  -- Note: Logging is now done via MonadWriter - no handler needed
    } deriving Generic

-- | Configuration for numerical convergence
--
-- This can be used to specify convenient defaults for convergence in e.g. L2
data ConvergConfig t a = ConvergConfig {
    ccStateSummary :: [t] -> a     -- ^ Produce a summary from a list of state projections
  , ccStateConverging :: a -> Bool  -- ^ Are we converging ?
  , ccStateDiverging :: a -> a -> Bool -- ^ Are we diverging ?
  , ccStateFinal :: t -> Bool -- ^ Has the state converged ?  
                                   }

convergenceL2 :: Normed v => (v -> Bool) -> ConvergConfig v (Magnitude v)
convergenceL2 = ConvergConfig norm2Diff nearZero (>) 

-- -- | Build an 'IterConfig'
-- mkIterConfig :: String
--              -> Int
--              -> Int
--              -> (s -> t)
--              -> ([t] -> a)
--              -> (a -> Bool)
--              -> (a -> a -> Bool)
--              -> (t -> Bool)
--              -> Handler m (WithSeverity msg)             
--              -> (s -> (Severity, msg))
--              -> IterConfig s t msg m a
-- mkIterConfig = IterConfig


-- class MonadState s m => MonadStateBuffer s m where
--   -- getBuffer :: s -> m (Maybe [a])
--   -- getBuffer :: s -> m (LoopState a)
--   getBuffer :: m (LoopState s)
--   putBuffer :: LoopState s -> m a

-- baz n f = do
--   -- s <- get
--   -- sb <- getBuffer s
--   sb <- getBuffer
--   let sb' = f $ getBuffers n sb
--   putBuffer sb'




updBuffer n snew (LoopState i ls s)
  | n <= 0 = Nothing
  | length ls < n = Just $ LoopState (i+1) (s : ls) snew
  | otherwise = Just $ LoopState (i+1) (s : take n ls) snew

data StateBuffer s = StateBuffer { sbPrevStates :: [s]
                                 , sbCurrentState :: s } deriving (Eq, Show)

-- instance Semigroup (StateBuffer s) where

initStateBuffer :: s -> StateBuffer s
initStateBuffer = StateBuffer []

-- reconstructStateBuffer n (StateBuffer _ ls s)
--   | length ls < n || n <= 0 = Nothing
--   | otherwise = Just buffer where
--       buffer = s : take n ls

-- mkStateBuffer :: Int -> [s] -> s -> Maybe (StateBuffer s)
-- mkStateBuffer n ls s
--   | length ls < n = Nothing
--   | length ls == n = Just $ StateBuffer ls s
--   | otherwise = Just $ StateBuffer (take n ls) s

  
-- | A record to keep track of the current iteration, a list of the prior states and the current state.
data LoopState s = LoopState { lsCounter :: !Int
                             , lsPrevStates :: [s]
                             , lsCurrentState :: s } deriving (Eq, Show)

-- | Reconstruct state buffer for convergence/divergence estimation
getBuffers :: Int -> LoopState a -> Maybe [a]
getBuffers n (LoopState _ ls s)
  | length ls < n || n <= 0 = Nothing
  | otherwise = Just buffer where
      buffer = s : take n ls
      
-- | Construct the initial LoopState with 'lsCounter' = 0, 'lsPrevStates' = []
mkLoopState :: s -> LoopState s
mkLoopState = LoopState 0 []

-- | Configurable iteration combinator, with convergence monitoring and logging via MonadWriter
modifyInspectGuardedM :: forall m a t s w. (MonadThrow m, MonadWriter w m, Show a, Typeable a, Show t, Typeable t) =>
                         ConvergConfig t a 
                      -> IterConfig s t
                      -> (s -> m s)
                      -> s
                      -> m s
modifyInspectGuardedM (ConvergConfig sf qconverg qdiverg qfinal) r f x0
  | nitermax > 0 = run
  | otherwise = throwM (NonNegError fname nitermax)
  where
    (IterConfig fname nitermax lwindow pf) = r
    updState snew (LoopState i lss s) = LoopState (i + 1) lssUpd snew
      where
        lss' = s : lss
        lssUpd | length lss < lwindow = lss'
               | otherwise            = take lwindow lss'    
    run :: m s
    run = do
      let s0 = mkLoopState x0 
      -- Run in WriterT to collect logs, then discard them
      ((aLast, sLast), _logs) <- runWriterT $ runIterativeT r s0 (loop :: IterativeT (IterConfig s t) (LoopState s) (WriterT w m) (Either (ConvergenceStatus s a) s))
      let i = lsCounter sLast
      case aLast of
        Left (NotConverged y) -> throwM $ NotConvergedE fname nitermax (pf y)
        Left (Diverging qi qt) -> throwM $ DivergingE fname i qi qt
        Right x -> pure x
        -- _ -> throwM $ IterE fname "asdf"
    loop :: IterativeT (IterConfig s t) (LoopState s) (WriterT w m) (Either (ConvergenceStatus s a) s)
    loop = do
      s@(LoopState i _ x) <- get
      y <- lift $ lift $ f x 
      let
        s' = updState y s        
        status = case getBuffers lwindow s' of
          Nothing -> BufferNotReady
          Just buffer -> stat
            where
              llf = pf `map` buffer
              qi = sf $ init llf  -- summary of [lwindow + 1 .. 0] states
              qt = sf $ tail llf  -- "       "  [lwindow     .. 1] states
              stat | qdiverg qi qt && not (qconverg qi) = Diverging qi qt
                   | qconverg qi || qfinal (pf y) =       Converged y
                   | i == nitermax - 1 =                  NotConverged y
                   | otherwise =                          Converging 
      case status of
        BufferNotReady -> do
          put s'
          loop 
        Converging -> do
          -- Could log here via MonadWriter tell if needed
          put s'
          loop 
        Diverging qi qt -> 
          pure $ Left (Diverging qi qt) 
        Converged qi -> 
          pure $ Right qi
        NotConverged yy -> 
          pure $ Left (NotConverged yy)
    



-- -- baz :: StateBuffer s m => (s -> Maybe [a] -> s) -> m ()
-- baz f = do
--   s <- get
--   sb <- getBuffer s
--   let s' = f sb
--   put s'




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

modifyUntilM_ :: MonadState s m => (s -> Bool) -> (s -> m s) -> m s
modifyUntilM_ q f = do
  x <- get
  y <- f x
  if q y
    then pure y
    else do
        put y
        modifyUntilM_ q f
    
-- -- | modifyUntil with optional iteration logging to stdout
-- -- modifyUntil' :: MonadLog String m =>
-- --    IterationConfig a b -> (a -> Bool) -> (a -> a) -> a -> m a
-- modifyUntil' config q f x0 = modifyUntilM' config q (pure . f) x0


-- modifyUntilM' :: MonadLog String m =>
--    IterationConfig a b -> (a -> Bool) -> (a -> m a) -> a -> m a
modifyUntilM' config q f x0 = execStateT (go 0) x0 where
  -- logf ii = (Informational, unwords ["Iteration", show ii, "\n"])
  go i = do
   x <- get
   y <- lift $ f x
   -- logWith (icLogWith config) i
   put y
   if q y
     then return y
     else go (i + 1)




-- -- | `untilConvergedG0` is a special case of `untilConvergedG` that assesses convergence based on the L2 distance to a known solution `xKnown`
-- -- untilConvergedG0 :: (Show p, MonadThrow m, Typeable v, Typeable (Magnitude v), Normed v) =>
-- --                     IterConfig s v msg m a
-- --                  -> v -> (s -> s) -> s -> m s
-- untilConvergedG0 config xKnown f x0 = 
--   modifyInspectGuarded config' f x0
--    where
--     config' = config {
--         icStateSummary = norm2Diff
--       , icStateConverging = nearZero
--       , icStateDiverging = (>)
--       , icStateFinal = \s -> nearZero $ norm2 (xKnown ^-^ s)
--       }
  


-- | This function makes some default choices on the `modifyInspectGuarded` machinery: convergence is assessed using the squared L2 distance between consecutive states, and divergence is detected when this function is increasing between pairs of measurements.
-- untilConvergedG :: (Show v, Epsilon v, MonadThrow m, Typeable v, Typeable (Magnitude v), Normed v) =>
--                    Handler m (WithSeverity msg)
--                 -> String
--                 -> Int
--                 -> Int
--                 -> (s -> v)
--                 -> (s -> (Severity, msg))
--                 -> (s -> s)
--                 -> s
--                 -> m s
-- untilConvergedG fh fname nitermax lwindow fp flog =
--   modifyInspectGuarded config
--   where
--     config = mkL2ConvergenceIterConf fname nitermax lwindow fp qfinal fh flog
--     qfinal = nearZero



-- untilConvergedGM fname config =
--   modifyInspectGuardedM fname config norm2Diff nearZero (>)



-- -- | Create a configuration for L2 convergence:
-- --
-- -- state summary : squared distance of vector sequence
-- -- convergence criterion : " ~= zero
-- -- divergence criterion : current summary is > previous one 
-- mkL2ConvergenceIterConf :: Normed v =>
--                            String      -- ^ Function name
--                         -> Int         -- ^ Max # iterations
--                         -> Int         -- ^ Buffer size
--                         -> (s -> v)    -- ^ State projection 
--                         -> (v -> Bool) -- ^ Termination criterion
--                         -> Handler m (WithSeverity msg)  -- ^ Logging handler
--                         -> (s -> (Severity, msg)) -- ^ Log formatting
--                         -> IterConfig s v msg m (Magnitude v)
-- mkL2ConvergenceIterConf fname nitermax lwindow fp =
--   mkIterConfig fname nitermax lwindow fp norm2Diff nearZero (>)


-- -- | Pure version of 'modifyInspectGuardedM'
-- modifyInspectGuarded :: (MonadThrow m, Show t, Show a, Typeable t, Typeable a) =>
--                         Handler m (WithSeverity msg) -- ^ Logging handler
--                      -> IterConfig s t msg a         -- ^ Configuration
--                      -> (s -> s)                     -- ^ State evolution
--                      -> s                            -- ^ Initial state
--                      -> m s                          -- ^ Final state
-- modifyInspectGuarded config f x0 = modifyInspectGuardedM config (pure . f) x0
  



-- * REMOVED LOGGING FUNCTIONS
-- Logging levels (Severity, WithSeverity) have been removed.
-- Use MonadWriter directly with tell for logging.



-- | Some useful combinators


-- | Apply a function over a range of integer indices, zip the result with it and filter out the almost-zero entries
onRangeSparse :: Epsilon b => (Int -> b) -> [Int] -> [(Int, b)]
onRangeSparse f ixs = foldr ins [] ixs where
  ins x xr | isNz (f x) = (x, f x) : xr
           | otherwise = xr

-- | ", monadic version
onRangeSparseM :: (Epsilon b, Foldable t, Monad m) =>
     (a -> m b) -> t a -> m [(a, b)]
onRangeSparseM f ixs = unfoldZipM mf f ixs where
  mf x = isNz <$> f x
  


unfoldZipM0 :: (Foldable t, Monad m) =>
     (a -> Bool) -> (a -> b) -> t a -> m [(a, b)]
unfoldZipM0 q f = unfoldZipM (pure . q) (pure . f)


unfoldZipM :: (Foldable t, Monad m) =>
     (a -> m Bool) -> (a -> m b) -> t a -> m [(a, b)]
unfoldZipM q f ixs = foldrM insf [] ixs where
  insf x xr = do
    qx <- q x
    if qx
    then do
      y <- f x
      pure $ (x, y) : xr
    else pure xr

-- | A combinator I don't know how to call
combx :: Functor f => (a -> b) -> (t -> f a) -> t -> f b
combx g f x = g <$> f x 


          
    



-- | Helpers

sqDiffPairs :: Num a => (v -> v -> a) -> [v] -> a
sqDiffPairs f uu = sqDiff f (init uu) (tail uu)

sqDiff :: Num a => (u -> v -> a) -> [u] -> [v] -> a
sqDiff f uu vv = sum $ zipWith f uu vv



-- | Relative residual
relRes :: (Normed t, LinearVectorSpace t) =>
     MatrixType t -> t -> t -> Magnitude t
relRes aa b x = n / d where
  n = norm2 $ (aa #> x) ^-^ b
  d = norm2 b


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



-- | NB: use it on a list of >= 2 elements !!
norm2Diff :: Normed v => [v] -> Magnitude v 
norm2Diff v = sum $ zipWith f va vb where
  f v1 v2 = norm2 $ v1 ^-^ v2
  va = init v
  vb = tail v

