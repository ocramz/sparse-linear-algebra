{-# language DeriveFunctor, GeneralizedNewtypeDeriving, MultiParamTypeClasses, FlexibleInstances #-}
module Control.Iterative.Internal (IterativeT(..), runIterativeT, MonadLog(..)) where

import Control.Monad.Reader (MonadReader(..))
import Control.Monad.State.Strict (MonadState(..), get, put)
import Control.Monad.Trans.Class (MonadTrans(..), lift)
import Control.Monad.Trans.State.Strict (StateT(..), runStateT)
import Control.Monad.Trans.Reader (ReaderT(..), runReaderT)
import Control.Monad.Writer.Strict (MonadWriter(..), WriterT(..), runWriterT)
import Control.Monad.Catch (MonadThrow(..), throwM)


-- | Pure logging typeclass, replacing logging-effect's MonadLog
class Monad m => MonadLog msg m where
  logMessage :: msg -> m ()

-- | Iterative algorithms need configuration, state and logging; here we use a transformer stack of ReaderT + StateT + WriterT.
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
  IterativeT { unIterativeT :: ReaderT c (StateT s (WriterT [msg] m)) a } deriving (Functor, Applicative, Monad, MonadReader c, MonadState s)

instance Monad m => MonadLog msg (IterativeT c msg s m) where
  logMessage msg = IterativeT $ lift $ lift $ tell [msg]

instance MonadTrans (IterativeT c msg s) where
  lift = liftIterativeT

liftIterativeT :: Monad m => m a -> IterativeT c msg s m a
liftIterativeT = IterativeT . lift . lift . lift
    
instance MonadThrow m => MonadThrow (IterativeT c msg s m) where
  throwM e = lift $ throwM e

-- | Run an 'IterativeT' computation, return result, final state, and log messages
runIterativeT :: Monad m
              => r  -- ^ Configuration
              -> s  -- ^ Initial state
              -> IterativeT r message s m a 
              -> m ((a, s), [message]) -- ^ ((result, final state), log messages)
runIterativeT c x0 m =
  runWriterT (runStateT (runReaderT (unIterativeT m) c) x0)

execIterativeT :: (Monad m, Functor m) =>
                  r
               -> s
               -> IterativeT r message s m a
               -> m s
execIterativeT c x0 m = snd . fst <$> runIterativeT c x0 m


      

