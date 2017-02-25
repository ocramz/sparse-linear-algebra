module Perf where


import Numeric.LinearAlgebra.Sparse
import Numeric.LinearAlgebra.Class
import Numeric.Eps
import qualified Data.Matrix.MatrixMarket as MM

import Data.VectorSpace hiding (magnitude)

import Data.Foldable
import Data.Scientific

fnameMat = "test/data/e05r0000.mtx"
fnameRhs = "test/data/e05r0000_rhs1.mtx"


-- main = return ()

main = do
  MM.RMatrix (nrow,ncol) _ _ dat <- MM.readMatrix fnameMat
  MM.RArray _ _ datr <- MM.readArray fnameRhs -- dim == ncol
  let b  = fromListSV ncol $ buildSparse datr
      aa = fromListSM (nrow,ncol) $ toD3 <$> dat
      -- nrow' = nrow
  -- return aa
  -- aa <\> b
  x <- linSolve0 BCG_ aa b (fromListSV nrow $ indexed $ replicate nrow 0.1)
  -- kappa <- conditionNumberSM aa
  let res = (aa #> x) ^-^ b
  return $ norm2 res
  
  
-- toD2 :: (a, Scientific) -> (a, Double)
-- toD2 (i, x) = (i, toRealFloat x)
toD3 :: Num i => (i, i, Scientific) -> (i, i, Double)
toD3 (i, j, x) = (i - 1, j - 1, toRealFloat x)

buildSparse :: (Num i, Epsilon t, RealFloat t) => [Scientific] -> [(i, t)]
buildSparse xxs = go 1 xxs where
  go i (x:xs) = let x' = toRealFloat x
                in if isNz x'
                   then (i - 1, x') : go (i+1) xs
                   else go (i+1) xs
  go _ [] = []
                                  





indexed :: Num i => [t] -> [(i, t)]
indexed xxs = go 0 xxs where
  go i (x:xs) = (i,x):go (i+1) xs
  go _ [] = []
