module Control.Exception.Common where

-- import Control.Exception.Safe (Exception, MonadThrow, SomeException, throwM)
import Control.Exception
import Control.Monad.Catch (MonadThrow (..))
import Data.Typeable -- (TypeRep, Typeable, typeRep)

import Data.Sparse.Utils


data OutOfBoundsIndexError i = OutOfBoundsIndexError i deriving (Show, Eq, Typeable)
instance (Show i, Typeable i) => Exception (OutOfBoundsIndexError i)


checkIxBound :: MonadThrow m => Int -> UB -> m a -> m a
checkIxBound i n ff
  | not (inBounds0 n i) = throwM (OutOfBoundsIndexError i)
  | otherwise = ff

    
