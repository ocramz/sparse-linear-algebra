{-# language DeriveFunctor, GeneralizedNewtypeDeriving #-}
module Control.Iterative.Internal (IterativeT(..), runIterativeT) where

import Control.Monad.Reader (MonadReader(..))
import Control.Monad.State.Strict (MonadState(..), get, put)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Log (MonadLog(..), LoggingT(..), runLoggingT, Handler, logMessage)
import Control.Monad.Catch (MonadThrow(..), throwM)


-- | Iterative algorithms need configuration, state and logging; here we use a transformer stack of ReaderT + StateT + LoggingT (from `logging-effect`).
--
-- The idea is to compose the plumbing of iterative programs from 'MonadState', 'MonadReader' and 'MonadLog' instructions (the "specification" of the effects to be used), which is usually referred to as the "@mtl@ style".
--
-- Example usage:
--
-- @
-- mkIterativeT :: Monad m =>
--           (a -> s -> message)
--        -> (s -> r -> (a, s))
--        -> IterativeT r message s m a
-- mkIterativeT flog fs = 'IterativeT' $ do
--   s <- 'get'
--   c <- 'ask'
--   let (a, s') = fs s c
--   'logMessage' $ flog a s'
--   'put' s'
--   return a
-- @
--
newtype IterativeT c msg s m a =
  IterativeT { unIterativeT :: ReaderT c (StateT s (LoggingT msg m)) a } deriving (Functor, Applicative, Monad, MonadReader c, MonadState s, MonadLog msg)

instance MonadTrans (IterativeT c msg s) where
  lift = liftIterativeT

liftIterativeT :: Monad m => m a -> IterativeT c msg s m a
liftIterativeT m = IterativeT . lift . lift $ LoggingT mlog
  where mlog = ReaderT (const m)
    
instance MonadThrow m => MonadThrow (IterativeT c msg s m) where
  throwM e = lift $ throwM e

-- | Run an 'IterativeT' computation, return result and final state
runIterativeT :: Handler m message -- ^ Logging handler
              -> r  -- ^ Configuration
              -> s  -- ^ Initial state
              -> IterativeT r message s m a 
              -> m (a, s) -- ^ (result, final state)
runIterativeT lh c x0 m =
  runLoggingT (runStateT (runReaderT (unIterativeT m) c) x0) lh

execIterativeT :: Functor m =>
                  Handler m message
               -> r
               -> s
               -> IterativeT r message s m a
               -> m s
execIterativeT lh c x0 m = snd <$> runIterativeT lh c x0 m


      

