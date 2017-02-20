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

pprintDefaults :: PPrintOptions
pprintDefaults = PPOpts 1 3


-- | printf format string for a Fractional
prepD :: (Ord t, Fractional t) => PPrintOptions -> t -> String
prepD opts x = pstr -- printf pstr x
  where
  pstr | abs x > 10 || abs x < 0.5 = s ++ "e"
       | otherwise = s ++ "f"
    where
      s = concat ["%" , show ni, ".", show nd]
  nd = pprintDecimals opts
  ni = pprintDigits opts

-- | printf format string for a Complex Fractional
prepC :: (Epsilon t, Ord t) => PPrintOptions -> Complex t -> String
prepC opts (r :+ i) = prepD opts r ++ oi where
    oi | isNz i = concat [s, "i ", prepD opts i']
       | otherwise = []
    s | signum i >= 0 = " + "
      | otherwise = " - "
    i' = abs i




