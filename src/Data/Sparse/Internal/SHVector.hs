module Data.Sparse.Internal.SHVector where

import qualified Data.Vector.Hybrid as VH
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V

import Control.Arrow (first, second)
import Data.Complex
-- import Control.DeepSeq

import Data.VectorSpace

import Numeric.LinearAlgebra.Class



-- | A sparse vector that uses Vector.Hybrid internally. This means that a vector of (Int, a) is represented by an _unboxed_ vector of Int's and a boxed vector of a's
data SHVector a = SHV {-# UNPACK #-} !Int !(VH.Vector VU.Vector V.Vector (Int, a)) deriving (Eq)

instance Show a => Show (SHVector a) where
  show (SHV n v) = unwords ["SHV (", show n,"),", show nz,"NZ:", show v]
    where
      nz = VH.length v      


-- -- instance Functor SHVector where
-- --   fmap f (SHV n v) = SHV n (fmap (second f) v)

-- fromList :: Int -> [(Int, a)] -> SHVector a
-- fromList n ll = SHV n (VH.fromList ll)
