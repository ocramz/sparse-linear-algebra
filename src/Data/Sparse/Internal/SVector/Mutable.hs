module Data.Sparse.Internal.SVector.Mutable where

import qualified Data.Vector as V 
import qualified Data.Vector.Mutable as VM

import Control.Monad.Primitive

data SMVector m a = SMV { smvDim :: {-# UNPACK #-} !Int,
                          smvIx :: V.Vector Int,
                          smvVal :: VM.MVector (PrimState m) a }


-- instance Show a => Show (SMVector m a) where
--   show (SMV n ix v) = unwords ["SMV (",show n,"),",show nz,"NZ:",show (V.zip ix v)]
--     where nz = V.length ix

-- instance Functor (SMVector m) where
--   fmap f (SMV n ix v) = SMV n ix (fmap f v)

-- instance Foldable (SMVector m) where
--   foldr f z (SMV _ _ v) = foldr f z v

-- instance Traversable (SMVector m) where
--   traverse f (SMV n ix v) = SMV n ix <$> traverse f v


