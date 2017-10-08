{-# language FlexibleContexts #-}
{-# language ScopedTypeVariables #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.LinearAlgebra.Accelerate
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  BSD3 (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- `accelerate` instances for linear algebra
--
-----------------------------------------------------------------------------
module Numeric.LinearAlgebra.Accelerate where


import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate
          (Acc, Array, Vector, Segments, DIM1, DIM2, Exp, Any(Any), All(All), Z(Z), (:.)((:.)))
-- import Data.Array.Accelerate.IO (fromVectors, toVectors)
-- import Data.Vector
-- import Data.Vector.Algorithms.Merge
import Data.Array.Accelerate.Interpreter (run)

-- | Vector as newtype
-- newtype Vector e = Vector (Array DIM1 e) deriving (Eq, Show)
-- | segments : vector of segment lengths 
-- newtype Segments i = Segments (Vector i) deriving (Eq, Show)









