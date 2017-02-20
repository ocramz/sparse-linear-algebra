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



data PPrintOptions =
  PPOpts {
     pprintDigits :: Int,
     pprintDecimals :: Int } deriving (Eq, Show)

pprintDefaults = PPOpts 1 3



prepD opts x = pstr -- printf pstr x
  where
  pstr | abs x > 10 || abs x < 0.5 = s ++ "e"
       | otherwise = s ++ "f"
    where
      s = concat ["%" , show ni, ".", show nd]
  nd = pprintDecimals opts
  ni = pprintDigits opts



-- printfComplex :: (PrintfArg t, Epsilon t, Ord t) =>
--      PPrintOptions -> Complex t -> String
prepC opts (r :+ i) = prepD opts r ++ oi where
    oi | isNz i = concat [s, "i ", prepD opts i']
       | otherwise = []
    s | signum i >= 0 = " + "
      | otherwise = " - "
    i' = abs i

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

