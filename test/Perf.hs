module Perf where


import Numeric.LinearAlgebra.Sparse
import qualified Data.Matrix.MatrixMarket as MM

import Data.Foldable
import Data.Scientific

fnameMat = "perf/data/e05r0000.mtx"
fnameRhs = "perf/data/e05r0000_rhs1.mtx"


-- main = return ()

main = do
  MM.RMatrix (nrow,ncol) _ _ dat <- MM.readMatrix fnameMat
  MM.RArray _ _ datr <- MM.readArray fnameRhs -- dim == ncol
  let b  = fromListSV ncol $ buildSparse datr
      aa = fromListSM (nrow,ncol) $ toD3 <$> dat
      -- nrow' = nrow
  -- return aa
  -- aa <\> b
  x <- linSolve0 CGS_ aa b (fromListSV nrow $ indexed $ replicate nrow 0.1)
  -- kappa <- conditionNumberSM aa
  return x
  
  
-- toD2 :: (a, Scientific) -> (a, Double)
-- toD2 (i, x) = (i, toRealFloat x)
toD3 :: Num a => (a, a, Scientific) -> (a, a, Double)
toD3 (i, j, x) = (i - 1, j - 1, toRealFloat x)

buildSparse :: Num a => [Scientific] -> [(a, Double)]
buildSparse xxs = go 1 xxs where
  go i (x:xs) | isNz x = (i - 1, toRealFloat x) : go (i+1) xs
              | otherwise = go (i+1) xs
  go _ [] = []

isNz :: (Ord a, Fractional a) => a -> Bool
isNz x = abs x < 1e-16

indexed :: Num t1 => [t] -> [(t1, t)]
indexed xxs = go 0 xxs where
  go i (x:xs) = (i,x):go (i+1) xs
  go _ [] = []
