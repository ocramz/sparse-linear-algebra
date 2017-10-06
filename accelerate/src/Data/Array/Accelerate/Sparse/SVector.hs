module Data.Array.Accelerate.Sparse.SVector where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)))

-- | sparse vectors
newtype SVector e = SVector (Array DIM1 (Int, e)) deriving (Eq, Show)
