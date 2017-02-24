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

pdef :: PPrintOptions
pdef = PPOpts 1 3




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
    oi | isNz i = concat [s, prepD opts i', "i"]
       | otherwise = []
    s | signum i >= 0 = "+"
      | otherwise = "-"
    i' = abs i


-- | printf for a list of values
--
-- > printN prepD (PPOpts 1 3) [1,pi]
printDN opts xl
  | null xl = printf "\n"
  | null xs = printf (prepD opts x) x
  | n==1 = let [x1]=xs in printf (pr x) x
  -- | n==1 = let [x1]=xs in printf s x x1
  -- | n==2 = let [x1,x2]=xs in printf s x x1 x2
  -- | n==3 = let [x1,x2,x3]=xs in printf s x x1 x2 x3
  -- | n==4 = let [x1,x2,x3,x4]=xs in printf s x x1 x2 x3 x4
  -- | otherwise = let [x1,x2,x3,x4,_]=xs
  --                   xfin=last xs in printf s' x x1 x2 x3 x4 xfin
  where
    (x:xs) = xl
    n = length xs
    pr = prepD opts



-- > printCN pdef [(1:+pi), (3.5:+4.3), (pi:+(-3.4))]     
printCN opts xl
  | null xl = printf "\n"
  | null xs = printf (pr x) (re x) (aim x)
  | n==1 = let [x1]=xs
           in printf (commas (pr <$> xl)++"\n") (re x) (aim x) (re x1) (aim x1)
  | n==2 = let [x1,x2]=xs
           in printf (commas (pr <$> xl) ++"\n") (re x) (aim x) (re x1) (aim x1) (re x2) (aim x2)
  | n==3 = let [x1,x2,x3]=xs
           in printf (commas (pr <$> xl)++"\n") (re x) (aim x) (re x1) (aim x1) (re x2) (aim x2) (re x3) (aim x3)
  | n==4 = let [x1,x2,x3,x4]=xs
           in printf (commas (pr <$> xl)++"\n") (re x) (aim x) (re x1) (aim x1) (re x2) (aim x2) (re x3) (aim x3) (re x4) (aim x4)
  | otherwise = let xs@[x,x1,x2,x3]=take 4 xl
                    xfin = last xl
                in printf (commas (pr <$> xs) ++ ", ... , " ++ pr xfin++"\n") (re x) (aim x) (re x1) (aim x1) (re x2) (aim x2) (re x3) (aim x3) (re xfin) (aim xfin)
  where
    (x:xs) = xl
    pr = prepC opts
    n = length xs

commas :: [String] -> String    
commas = intercalate ", "

re :: Complex a -> a
re = realPart
aim :: Num a => Complex a -> a
aim = abs . imagPart
