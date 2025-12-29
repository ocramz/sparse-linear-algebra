{-# language DeriveFunctor, GeneralizedNewtypeDeriving, MultiParamTypeClasses, FlexibleInstances, UndecidableInstances, FunctionalDependencies #-}
module Control.Iterative.Internal (IterativeT(..), runIterativeT) where

import Control.Monad.Reader (MonadReader(..))
import Control.Monad.State.Strict (MonadState(..), get, put)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Writer.Class (MonadWriter(..))
import Control.Monad.Catch (MonadThrow(..), throwM)


-- | Iterative algorithms need configuration, state and logging; here we use a transformer stack of ReaderT + StateT.
-- Logging is done via any MonadWriter instance from mtl.
--
-- The idea is to compose the plumbing of iterative programs from 'MonadState', 'MonadReader' and 'MonadWriter' instructions (the "specification" of the effects to be used), which is usually referred to as the "@mtl@ style".
--
-- Example usage:
--
-- @
-- mkIterativeT :: (Monad m, MonadWriter [msg] m) =>
--           (a -> s -> msg)
--        -> (s -> r -> (a, s))
--        -> IterativeT r s m a
-- mkIterativeT flog fs = 'IterativeT' $ do
--   s <- 'get'
--   c <- 'ask'
--   let (a, s') = fs s c
--   'tell' [flog a s']
--   'put' s'
--   return a
-- @
--
newtype IterativeT c s m a =
  IterativeT { unIterativeT :: ReaderT c (StateT s m) a } deriving (Functor, Applicative, Monad, MonadReader c, MonadState s)

instance MonadTrans (IterativeT c s) where
  lift = liftIterativeT

liftIterativeT :: Monad m => m a -> IterativeT c s m a
liftIterativeT = IterativeT . lift . lift
    
instance MonadThrow m => MonadThrow (IterativeT c s m) where
  throwM e = lift $ throwM e

instance MonadWriter w m => MonadWriter w (IterativeT c s m) where
  tell w = IterativeT $ lift $ lift $ tell w
  listen m = IterativeT $ ReaderT $ \r -> StateT $ \s -> do
    ((a, s'), w) <- listen $ runStateT (runReaderT (unIterativeT m) r) s
    return ((a, w), s')
  pass m = IterativeT $ ReaderT $ \r -> StateT $ \s -> pass $ do
    ((a, f), s') <- runStateT (runReaderT (unIterativeT m) r) s
    return ((a, s'), f)

-- | Run an 'IterativeT' computation in any monad with MonadWriter support
runIterativeT :: (Monad m, MonadWriter w m)
              => r  -- ^ Configuration
              -> s  -- ^ Initial state
              -> IterativeT r s m a 
              -> m (a, s) -- ^ (result, final state) - logs are accumulated via MonadWriter
runIterativeT c x0 m =
  runStateT (runReaderT (unIterativeT m) c) x0

execIterativeT :: (Monad m, Functor m, MonadWriter w m) =>
                  r
               -> s
               -> IterativeT r s m a
               -> m s
execIterativeT c x0 m = snd <$> runIterativeT c x0 m


      

