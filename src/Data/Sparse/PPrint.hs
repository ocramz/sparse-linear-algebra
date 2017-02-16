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
import Text.Printf

import Numeric.Eps



class PrintDense a where
  -- | Pretty-print with a descriptive header
  prd :: a -> IO ()
  -- | Pretty-print with no header
  prd0 :: a -> IO ()

newline :: IO ()
newline = putStrLn ""  




printfDouble :: (PrintfArg t, PrintfType t1) => PPrintOptions -> t -> t1
printfDouble opts x = printf pstr x where
  pstr = concat ["%" , show ni, ".", show nd, "f"]
  nd = pprintDecimals opts
  ni = pprintDigits opts - nd

data PPrintOptions =
  PPrintOptions {
     pprintDigits :: Int,
     pprintDecimals :: Int } deriving (Eq, Show)

pprintDefaults = PPrintOptions 5 2




-- | Cleaner way to display Complex values

newtype C a = C {unC :: Complex a} deriving Eq
instance (Num a, Epsilon a, Ord a, Show a) => Show (C a) where
  show (C (r :+ i)) = unwords [show r, oi] where
    oi | isNz i = unwords [s, show i' ++ "j"]
       | otherwise = []
    s | signum i >= 0 = "+"
      | otherwise = "-"
    i' = abs i

-- c0, c1 :: C Double
-- c0 = C $ 1 :+ (-2) 
-- c1 = C $ pi :+ 0

