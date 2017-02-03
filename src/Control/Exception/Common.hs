module Control.Exception.Common where

-- import Control.Exception.Safe (Exception, MonadThrow, SomeException, throwM)
import Control.Exception
import Control.Monad.Catch (MonadThrow (..))
import Data.Typeable -- (TypeRep, Typeable, typeRep)

import Data.Sparse.Utils



data PartialFunctionError = EmptyList String deriving (Eq, Typeable)
instance Show PartialFunctionError where
  show (EmptyList s) = unwords [s, ": empty list"]
instance Exception PartialFunctionError


-- | Input error

data InputError = NonNegError String Int deriving (Eq, Typeable)
instance Show InputError where
  show (NonNegError s i) = unwords [s, ": parameter must be nonnegative, instead I got", show i]
instance Exception InputError



-- | Out of bounds index error
data OutOfBoundsIndexError i = OOBIxError String i
                             | OOBIxsError String [i]
                             | OOBNoCompatRows String (i,i) deriving (Eq, Typeable)
instance Show i => Show (OutOfBoundsIndexError i) where
  show (OOBIxError e i) = unwords [e, ": index", show i,"out of bounds"]
  show (OOBIxsError e ixs) = unwords [e, ":, indices", show ixs, "out of bounds"]
  show (OOBNoCompatRows e ij) = unwords [e, ": no compatible rows for indices", show ij]
instance (Show i, Typeable i) => Exception (OutOfBoundsIndexError i)

checkIxBound :: MonadThrow m => String -> Int -> UB -> m a -> m a
checkIxBound e i n ff
  | not (inBounds0 n i) = throwM (OOBIxError e i)
  | otherwise = ff

    

-- | Operand size mismatch errors
data OperandSizeMismatch = DotSizeMismatch Int Int
                         | NonTriangularException String
                         | MatVecSizeMismatchException String (Int, Int) Int deriving (Eq, Typeable)
instance Show OperandSizeMismatch where
  show (DotSizeMismatch na nb) = unwords ["<.> : Incompatible dimensions : ", show na, show nb]
  show (NonTriangularException s )= unwords [s, ": Matrix must be triangular"]
  show (MatVecSizeMismatchException s dm o) = unwords [s, ": Matrix-vector dimensions are incompatible: Matrix is",show dm,", whereas vector is",show o]
instance Exception OperandSizeMismatch




-- | Matrix exceptions
data MatrixException i = HugeConditionNumber i deriving (Eq, Typeable)
instance Show i => Show (MatrixException i) where
  show (HugeConditionNumber x) = unwords ["Rank-deficient system: condition number", show x]
instance (Show i, Typeable i) => Exception (MatrixException i)






-- | Numerical iteration errors
data IterationException a = NotConverged String Int a 
                          | Diverging String Int a a    
                          deriving (Eq, Typeable)
instance Show a => Show (IterationException a) where
  show (NotConverged s niters x) = unwords [s, ": Could not converge within iteration budget", show niters, "; final state:", show x]
  show (Diverging s niters x0 x1) = unwords [s, ": Diverging at iteration #", show niters, "; latest state:",show x0, "; previous state:",show x1]
instance (Show a, Typeable a) => Exception (IterationException a)
