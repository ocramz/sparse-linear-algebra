{-# language FlexibleContexts, GeneralizedNewtypeDeriving, DeriveFunctor #-}
{-# language OverloadedStrings #-}
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
module Control.Iterative where

import Control.Applicative

import Control.Monad (when, replicateM)
import Control.Monad.Reader (MonadReader(..), asks)
import Control.Monad.State.Strict (MonadState(..), get, put, gets)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT, execStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Catch (Exception(..), MonadThrow(..), throwM)
import Control.Monad.Log (MonadLog(..), WithSeverity(..), Severity(..), renderWithSeverity, LoggingT(..), runLoggingT, Handler, logMessage, logError, logDebug, logInfo, logNotice)

-- import Data.Bool (bool)
import Data.Char (toUpper)

import Data.Typeable
import qualified Control.Exception as E (Exception, Handler)

import Data.Foldable (foldrM)

import Control.Exception.Common
import Numeric.LinearAlgebra.Class
import Numeric.Eps


data ConvergenceStatus a = BufferNotReady
                         | Converging
                         | Converged a
                         | Diverging a a
                         | NotConverged
                           deriving (Eq, Show, Typeable)
instance (Typeable a, Show a) => Exception (ConvergenceStatus a)

-- FIXME I don't like this; logging shouldn't be in IO, `iterationView` is too general, etc.
data IterationConfig a b =
  IterConf { numIterationsMax :: Int     -- ^ Max.# of iterations
           , printDebugInfo :: Bool      -- ^ Print iteration info to stdout 
           , iterationView :: a -> b     -- ^ Project state to a type `b`
           , printDebugIO :: b -> IO ()} -- ^ print function for type `b`
instance Show (IterationConfig a b) where
  show (IterConf n qd _ _) = unwords ["Max. # of iterations:",show n,", print debug information:", show qd]


-- | Iterative algorithms need configuration, state and logging; here we use a transformer stack of ReaderT / StateT / LoggingT (from `logging-effect`)
newtype IterativeT c msg s m a =
  IterativeT { unIterativeT :: ReaderT c (StateT s (LoggingT msg m)) a } deriving (Functor, Applicative, Monad, MonadReader c, MonadState s, MonadLog msg)

instance MonadTrans (IterativeT c msg s) where
  lift = liftIterativeT

liftIterativeT :: Monad m => m a -> IterativeT c msg s m a
liftIterativeT m = IterativeT . lift . lift $ LoggingT mlog
  where mlog = ReaderT (const m)
    

  
-- instance MonadThrow m => MonadThrow (IterativeT c msg s m) where
--   throwM e = lift $ throwM e

runIterativeT :: Handler m message
              -> r
              -> s
              -> IterativeT r message s m a
              -> m (a, s)
runIterativeT lh c x0 m =
  runLoggingT (runStateT (runReaderT (unIterativeT m) c) x0) lh

execIterativeT :: Functor m =>
                  Handler m message
               -> r
               -> s
               -> IterativeT r message s m a
               -> m s
execIterativeT lh c x0 m = snd <$> runIterativeT lh c x0 m

-- | Log with a function that computes a severity and a message from the input
logWith :: MonadLog (WithSeverity a) m => (p -> (Severity, a)) -> p -> m ()
logWith f x = logMessage (WithSeverity sev sevMsg) where
  (sev, sevMsg) = f x
      

mkIterativeT :: Monad m =>
          (a -> s -> message)
       -> (s -> r -> (a, s))
       -> IterativeT r message s m a
mkIterativeT flog fs = IterativeT $ do
  s <- get
  c <- ask
  let (a, s') = fs s c
  logMessage $ flog a s'
  put s'
  return a  










sqDiffPairs :: Num a => (v -> v -> a) -> [v] -> a
sqDiffPairs f uu = sqDiff f (init uu) (tail uu)

sqDiff :: Num a => (u -> v -> a) -> [u] -> [v] -> a
sqDiff f uu vv = sum $ zipWith f uu vv






data LoopState s = LoopState { lsCounter :: !Int
                             , lsPrevStates :: [s]
                             , lsCurrentState :: s } deriving (Eq, Show)

-- | Reconstruct state buffers for convergence/divergence estimation
getBuffers :: Int -> LoopState a -> Maybe ([a], [a])
getBuffers n (LoopState _ ls s)
  | length ls < n = Nothing
  | otherwise = Just (curr, prev) where
      curr = s : take (n - 1) ls
      prev = take n ls

mkLoopState :: s -> LoopState s
mkLoopState = LoopState 0 []



data ConvergenceStatus' a b =
  BufferNotReady'
  | Converging'
  | Converged' a
  | Diverging' a a
  | NotConverged' b
  deriving (Eq, Show, Typeable)
instance (Typeable a, Show a, Typeable b, Show b) => Exception (ConvergenceStatus' a b)

-- | Configuration data for the iterative process
data IterConfig s t msg a = IterConfig {
    icFunctionName :: String -- ^ Name of calling function, for logging purposes
  , icNumIterationsMax :: Int -- ^ Max # of iterations
  , icStateWindowLength :: Int -- ^ # of states used to assess convergence/divergence
  , icStateProj :: s -> t          -- ^ Project the state
  , icStateSummary :: [t] -> a     -- ^ Produce a summary from a list of states
  , icStateConverging :: a -> Bool  -- ^ Are we converging ?
  , icStateDiverging :: a -> a -> Bool -- ^ Are we diverging ?
  , icStateFinal :: t -> Bool -- ^ Has the state converged ?
  , icLogWith :: s -> (Severity, msg) -- ^ Compute log severity and message
    }


-- TODO : remove Show, Typeable constraints from 's'
modifyInspectGuardedM_Iter :: (MonadThrow m, Show s, Show a, Typeable s, Typeable a) =>
                              Handler m (WithSeverity msg)
                           -> IterConfig s t msg a
                           -> (s -> m s)
                           -> s
                           -> m a
modifyInspectGuardedM_Iter lh r@(IterConfig fname nitermax lwindow pf sf qconverg qdiverg qfinal lwf) f x0
  | nitermax > 0 = run
  | otherwise = throwM (NonNegError fname nitermax)
  where
    updStateN snew (LoopState i lss s) = LoopState (i + 1) lssUpd snew
      where
        lss' = s : lss
        lssUpd | length lss < lwindow = lss'
               | otherwise = take lwindow lss'    
    run = do
      let s0 = mkLoopState x0 
      (aLast, sLast) <- runIterativeT lh r s0 loop -- (loop 0 [])
      let i = lsCounter sLast
      case aLast of
        Left (NotConverged' y) -> throwM $ NotConvergedE fname nitermax y
        Left (Diverging' qi qt) -> throwM $ DivergingE fname i qi qt
        Right x -> pure x
    loop = do
      s@(LoopState i ll x) <- get
      y <- lift $ f x 
      let
        llf = pf `map` ll
        qi = sf $ init llf  -- summary of [lwindow + 1 .. 0] states
        qt = sf $ tail llf  -- "       "  [lwindow     .. 1] states
        status | length ll < lwindow =                BufferNotReady'
               | qdiverg qi qt && not (qconverg qi) = Diverging' qi qt        
               | qconverg qi || qfinal (pf y) =       Converged' qi
               | i == nitermax - 1 =                  NotConverged' y
               | otherwise =                          Converging'         
      case status of
        BufferNotReady' -> do
          let s' = updStateN y s
          -- logMessage $ WithSeverity Debug ("Buffer not ready") 
          put s'
          loop 
        Converging' -> do
          let s' = updStateN y s
          -- logMessage $ WithSeverity Informational "Converging"
          logWith lwf y          
          put s'
          loop 
        Diverging' qi qt -> 
          pure $ Left (Diverging' qi qt) -- (DivergingE fname i qi qt)
        Converged' qi -> 
          pure $ Right qi
        NotConverged' y -> 
          pure $ Left (NotConverged' y) -- (NotConvergedE fname nitermax y)           
    






-- modifyInspectGuardedM_Iter lh config fname sf qconverg qdiverg qfinal f x0 
--   | nitermax > 0 = execIterativeT lh r x0 (go 0 [])
--   | otherwise = throwM (NonNegError fname nitermax)
--   where
--     lwindow = 3
--     nitermax = numIterationsMax config
--     -- pf = iterationView config
--     checkConvergStatus y i ll
--       | length ll < lwindow = BufferNotReady
--       | qdiverg qi qt && not (qconverg qi) = Diverging qi qt        
--       | qconverg qi || qfinal (pf y) = Converged qi
--       | i == nitermax - 1 = NotConverged         
--       | otherwise = Converging
--       where llf = pf <$> ll
--             qi = sf $ init llf         -- summary of latest 2 states
--             qt = sf $ tail llf         -- "       "  previous 2 states
--     go i ll = do
--       x <- get
--       y <- lift $ f x
--       -- when (printDebugInfo config) $ do
--       --   logMessage $ unwords ["Iteration", show i]
--       -- when (printDebugInfo config) $ liftIO $ do
--       --   putStrLn $ unwords ["Iteration", show i]
--       --   printDebugIO config (pf y) 
--       case checkConvergStatus y i ll of
--         BufferNotReady -> do  
--           put y
--           let ll' = y : ll    -- cons current state to buffer
--           go (i + 1) ll'
--         Converged qi -> put y
--         Diverging qi qt -> do
--           put y
--           throwM (DivergingE fname i qi qt)
--         Converging -> do
--           put y
--           let ll' = init (y : ll) -- rolling state window
--           go (i + 1) ll'
--         NotConverged -> do
--           put y
--           throwM (NotConvergedE fname nitermax y)



  





bracketsUpp :: Show a => a -> String
bracketsUpp p = unwords ["[", map toUpper (show p), "]"]

withSeverity :: (t -> String) -> WithSeverity t -> String
withSeverity k (WithSeverity u a ) = unwords [bracketsUpp u, k a]



-- -- >>> renderWithSeverity id (WithSeverity Informational "Flux capacitor is functional")
-- -- [Informational] Flux capacitor is functional
-- renderWithSeverity
--   :: (a -> PP.Doc) -> (WithSeverity a -> PP.Doc)
-- renderWithSeverity k (WithSeverity u a) =
--   PP.brackets (PP.pretty u) PP.<+> PP.align (k a)



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
    
-- | modifyUntil with optional iteration logging to stdout
modifyUntil' :: MonadLog String m =>
   IterationConfig a b -> (a -> Bool) -> (a -> a) -> a -> m a
modifyUntil' config q f x0 = modifyUntilM' config q (pure . f) x0


modifyUntilM' :: MonadLog String m =>
   IterationConfig a b -> (a -> Bool) -> (a -> m a) -> a -> m a
modifyUntilM' config q f x0 = execStateT (go 0) x0 where
  pf = iterationView config
  go i = do
   x <- get
   y <- lift $ f x
   when (printDebugInfo config) $ do
     logMessage $ unwords ["Iteration", show i, "\n"]
   -- when (printDebugInfo config) $ liftIO $ do
   --   putStrLn $ unwords ["Iteration", show i, "\n"]
   --   printDebugIO config (pf y) 
   put y
   if q y
     then return y
     else go (i + 1)




-- | `untilConvergedG0` is a special case of `untilConvergedG` that assesses convergence based on the L2 distance to a known solution `xKnown`
untilConvergedG0 :: (Normed v, MonadThrow m, MonadLog String m, Typeable (Magnitude v), Typeable s, Show s) =>
     String
     -> IterationConfig s v
     -> v                    -- ^ Known value
     -> (s -> s)
     -> s
     -> m s
untilConvergedG0 fname config xKnown f x0 = 
  modifyInspectGuarded fname config norm2Diff nearZero (>) qfin f x0
   where
    qfin s = nearZero $ norm2 (xKnown ^-^ s)
  


-- | This function makes some default choices on the `modifyInspectGuarded` machinery: convergence is assessed using the squared L2 distance between consecutive states, and divergence is detected when this function is increasing between pairs of measurements.
untilConvergedG :: (Normed v, MonadThrow m, MonadLog String m, Typeable (Magnitude v), Typeable s, Show s) =>
        String
     -> IterationConfig s v
     -> (v -> Bool)
     -> (s -> s)               
     -> s 
     -> m s
untilConvergedG fname config =
  modifyInspectGuarded fname config norm2Diff nearZero (>)


-- | ", monadic version
untilConvergedGM ::
  (Normed v, MonadThrow m, MonadLog String m, Typeable (Magnitude v), Typeable s, Show s) =>
     String
     -> IterationConfig s v
     -> (v -> Bool)
     -> (s -> m s)          
     -> s
     -> m s
untilConvergedGM fname config =
  modifyInspectGuardedM fname config norm2Diff nearZero (>)





-- | `modifyInspectGuarded` is a high-order abstraction of a numerical iterative process. It accumulates a rolling window of 3 states and compares a summary `q` of the latest 2 with that of the previous two in order to assess divergence (e.g. if `q latest2 > q prev2` then the function throws an exception and terminates). The process ends by either hitting an iteration budget or by relative convergence, whichever happens first. After the iterations stop, the function then assesses the final state with a predicate `qfinal` (e.g. for comparing the final state with a known one; if this is not available, the user can just supply `const True`)
modifyInspectGuarded ::
  (MonadThrow m, MonadLog String m, Typeable s, Typeable a, Show s, Show a) =>
        String              -- ^ Calling function name
     -> IterationConfig s v -- ^ Configuration
     -> ([v] -> a)          -- ^ State summary array projection
     -> (a -> Bool)         -- ^ Convergence criterion
     -> (a -> a -> Bool)    -- ^ Divergence criterion
     -> (v -> Bool)         -- ^ Final state acceptance criterion
     -> (s -> s)            -- ^ State evolution
     -> s                   -- ^ Initial state
     -> m s                 -- ^ Final state
modifyInspectGuarded fname config sf qc qd qfin f x0 =
  modifyInspectGuardedM fname config sf qc qd qfin (pure . f) x0

  


-- | ", monadic version
modifyInspectGuardedM ::
  (MonadThrow m, MonadLog String m, Typeable s, Show s, Typeable a, Show a) =>
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
  | nitermax > 0 = execStateT (go 0 []) x0
  | otherwise = throwM (NonNegError fname nitermax)
  where
    lwindow = 3
    nitermax = numIterationsMax config
    pf = iterationView config
    checkConvergStatus y i ll
      | length ll < lwindow = BufferNotReady
      | qdiverg qi qt && not (qconverg qi) = Diverging qi qt        
      | qconverg qi || qfinal (pf y) = Converged qi
      | i == nitermax - 1 = NotConverged         
      | otherwise = Converging
      where llf = pf <$> ll
            qi = sf $ init llf         -- summary of latest 2 states
            qt = sf $ tail llf         -- "       "  previous 2 states
    go i ll = do
      x <- get
      y <- lift $ f x
      when (printDebugInfo config) $ do
        logMessage $ unwords ["Iteration", show i]
      -- when (printDebugInfo config) $ liftIO $ do
      --   putStrLn $ unwords ["Iteration", show i]
      --   printDebugIO config (pf y) 
      case checkConvergStatus y i ll of
        BufferNotReady -> do  
          put y
          let ll' = y : ll    -- cons current state to buffer
          go (i + 1) ll'
        Converged qi -> put y
        Diverging qi qt -> do
          put y
          throwM (DivergingE fname i qi qt)
        Converging -> do
          put y
          let ll' = init (y : ll) -- rolling state window
          go (i + 1) ll'
        NotConverged -> do
          put y
          throwM (NotConvergedE fname nitermax y)
             


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


-- qdiverg :: Ord a => a -> a -> Bool
-- qdiverg = (>)

norm2Diff [s1, s0] = norm2 (s1 ^-^ s0)
norm2Diff _ = 1/0






