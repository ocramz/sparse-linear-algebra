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
import Data.List
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
     pprintInt :: Int,
     pprintDec :: Int
     } deriving (Eq, Show)

pprintDefaults :: PPrintOptions
pprintDefaults = PPOpts 1 3




-- | printf format string for a Fractional
prepD :: (Ord t, Fractional t) => PPrintOptions -> t -> String
prepD opts x = pstr -- printf pstr x
  where
  pstr | abs x >= 10 || abs x < 0.5 = s ++ "e"
       | otherwise = s ++ "f"
    where
      s = concat ["%" , show ni, ".", show nd]
  nd = pprintDec opts
  ni = pprintInt opts

-- | printf format string for a Complex Fractional
prepC :: (Epsilon t, Ord t) => PPrintOptions -> Complex t -> String
prepC opts (r :+ i) = prepD opts r ++ oi where
    oi | isNz i = concat [s, "i ", prepD opts i']
       | otherwise = []
    s | signum i >= 0 = " + "
      | otherwise = " - "
    i' = abs i


-- | printf for a list of values
--
-- > printN prepD (PPOpts 1 3) [1,pi]
printN prepf opts xl
  | null xl = printf "\n"
  | null xs = printf (prepf opts x) x
  | n==1 = let [x1]=xs in printf s x x1
  | n==2 = let [x1,x2]=xs in printf s x x1 x2
  | n==3 = let [x1,x2,x3]=xs in printf s x x1 x2 x3
  | n==4 = let [x1,x2,x3,x4]=xs in printf s x x1 x2 x3 x4
  | otherwise = let [x1,x2,x3,x4,_]=xs
                    xfin=last xs in printf s' x x1 x2 x3 x4 xfin
  where
    (x:xs) = xl
    n = length xs
    s = unwords (replicate (n+1) (prepf opts x)) ++ "\n"
    s'= unwords [unwords (replicate n (prepf opts x)), "...", prepf opts x] ++ "\n"
  -- case length ll of 1 -> printf (prepf opts x) x where
  --                     x = head ll
  --                   2 -> printf (prepf opts )
  
