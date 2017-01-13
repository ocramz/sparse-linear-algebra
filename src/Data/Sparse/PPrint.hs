-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module Data.Sparse.PPrint where

import Data.Complex

class PrintDense a where
  prd :: a -> IO ()

newline :: IO ()
newline = putStrLn ""  








newtype C a = C {unC :: Complex a} deriving Eq
instance (Num a, Ord a, Show a) => Show (C a) where
  show (C (r :+ i)) = unwords [show r, s, show i' ++ "j"] where
    s | signum i >= 0 = "+"
      | otherwise = "-"
    i' = abs i

-- c0 = C $ 1 :+ (-2)
-- c1 = C $ 3 :+ 2
