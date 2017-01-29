module Control.Exception.Common where

-- import Control.Exception.Safe (Exception, MonadThrow, SomeException, throwM)
import Control.Exception
import Control.Monad.Catch (MonadThrow (..))
import Data.Typeable -- (TypeRep, Typeable, typeRep)

import Data.Sparse.Utils

-- | Out of bounds index error
data OutOfBoundsIndexError i = OOBIxError String i deriving (Eq, Typeable)
instance Show i => Show (OutOfBoundsIndexError i) where
  show (OOBIxError e i) = unwords [e, ": index", show i,"out of bounds"]
instance (Show i, Typeable i) => Exception (OutOfBoundsIndexError i)

checkIxBound :: MonadThrow m => String -> Int -> UB -> m a -> m a
checkIxBound e i n ff
  | not (inBounds0 n i) = throwM (OOBIxError e i)
  | otherwise = ff

    

-- | Inner product <.> between vectors of mismatched dimension
data DotSizeMismatch = DotSizeMismatch Int Int deriving (Eq, Typeable)
instance Show DotSizeMismatch where
  show (DotSizeMismatch na nb) = unwords ["<.> : Incompatible dimensions : ", show na, show nb]
instance Exception DotSizeMismatch



-- | Matrix-related errors

data MatrixException i = HugeConditionNumber i deriving (Eq, Typeable)
instance Show i => Show (MatrixException i) where
  show (HugeConditionNumber x) = unwords ["Rank-deficient system: condition number", show x]
instance (Show i, Typeable i) => Exception (MatrixException i)



-- | Numerical iteration errors
data IterationException = NotConverged  -- ^ residual norm is greater than tolerance
                        | Diverging     -- ^ residual norm increased
                        deriving (Show, Eq, Typeable)
instance Exception IterationException   
