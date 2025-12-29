{-# language DeriveFunctor, GeneralizedNewtypeDeriving, MultiParamTypeClasses, FlexibleInstances, UndecidableInstances, FunctionalDependencies #-}
module Control.Iterative.Internal (IterativeT(..), runIterativeT, MonadLog(..)) where

import Control.Monad.Reader (MonadReader(..))
import Control.Monad.State.Strict (MonadState(..), get, put)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Writer.Class (MonadWriter(..))
import Control.Monad.Catch (MonadThrow(..), throwM)


-- | Pure logging typeclass, using mtl's MonadWriter interface
class (Monad m, Monoid w) => MonadLog w m | m -> w where
  logMessage :: w -> m ()

-- | Iterative algorithms need configuration, state and logging; here we use a transformer stack of ReaderT + StateT.
-- Logging is done via any MonadWriter instance.
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
  IterativeT { unIterativeT :: ReaderT c (StateT s m) a } deriving (Functor, Applicative, Monad, MonadReader c, MonadState s)

-- Instance for MonadLog using the underlying MonadWriter
instance (Monad m, Monoid msg, MonadWriter msg m) => MonadLog msg (IterativeT c msg s m) where
  logMessage msg = IterativeT $ lift $ lift $ tell msg

instance MonadTrans (IterativeT c msg s) where
  lift = liftIterativeT

liftIterativeT :: Monad m => m a -> IterativeT c msg s m a
liftIterativeT = IterativeT . lift . lift
    
instance MonadThrow m => MonadThrow (IterativeT c msg s m) where
  throwM e = lift $ throwM e

instance MonadWriter w m => MonadWriter w (IterativeT c msg s m) where
  tell w = IterativeT $ lift $ lift $ tell w
  listen m = IterativeT $ ReaderT $ \r -> StateT $ \s -> do
    ((a, s'), w) <- listen $ runStateT (runReaderT (unIterativeT m) r) s
    return ((a, w), s')
  pass m = IterativeT $ ReaderT $ \r -> StateT $ \s -> pass $ do
    ((a, f), s') <- runStateT (runReaderT (unIterativeT m) r) s
    return ((a, s'), f)

-- | Run an 'IterativeT' computation in any monad with MonadWriter support
runIterativeT :: (Monad m, MonadWriter msg m)
              => r  -- ^ Configuration
              -> s  -- ^ Initial state
              -> IterativeT r msg s m a 
              -> m (a, s) -- ^ (result, final state) - logs are accumulated via MonadWriter
runIterativeT c x0 m =
  runStateT (runReaderT (unIterativeT m) c) x0

execIterativeT :: (Monad m, Functor m, MonadWriter msg m) =>
                  r
               -> s
               -> IterativeT r msg s m a
               -> m s
execIterativeT c x0 m = snd <$> runIterativeT c x0 m


      

