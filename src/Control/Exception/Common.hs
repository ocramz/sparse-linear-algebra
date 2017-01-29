module Control.Exception.Common where

-- import Control.Exception.Safe (Exception, MonadThrow, SomeException, throwM)
import Control.Exception
import Control.Monad.Catch (MonadThrow (..))
import Data.Typeable -- (TypeRep, Typeable, typeRep)

import Data.Sparse.Utils


-- | Input error

data InputError = NonNegError String Int deriving (Eq, Typeable)
instance Show InputError where
  show (NonNegError s i) = unwords [s, ": parameter must be nonnegative, instead I got", show i]
instance Exception InputError



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


data MatrixShapeException = NonTriangularException String 
                          deriving (Eq, Typeable)
instance Show MatrixShapeException where
  show (NonTriangularException s )= unwords [s, ": matrix must be triangular"]
instance Exception MatrixShapeException



-- | Numerical iteration errors
data IterationException a = NotConverged String Int a 
                          | Diverging String a a    
                          deriving (Eq, Typeable)
instance Show a => Show (IterationException a) where
  show (NotConverged s niters x) = unwords [s, ": Could not converge within iteration budget", show niters, "; final state:", show x]
  show (Diverging s x0 x1) = unwords [s, ": Diverging iterations", show x0, show x1]
instance (Show a, Typeable a) => Exception (IterationException a)
