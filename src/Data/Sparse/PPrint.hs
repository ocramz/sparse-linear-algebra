-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module Data.Sparse.PPrint (prd, prd0, PrintDense, newline,
                           PPrintOptions, prdef, prepD, prepC,
                           -- printDN,
                           printCN) where

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


-- | Pretty printing options: # of integer digits, # of decimal digits
data PPrintOptions =
  PPOpts {
     pprintInt :: Int,
     pprintDec :: Int
     } deriving (Eq, Show)

prdef :: PPrintOptions
prdef = PPOpts 1 2




-- | printf format string
prepD :: (Ord t, Epsilon t) => PPrintOptions -> t -> String
prepD opts x = pstr -- printf pstr x
  where
  pstr | nearZero x = "_"
       | abs x > 999 || abs x < 0.1 = s ++ "e"
       | otherwise = s ++ "f"
    where
      s = concat ["%" , show ni, ".", show nd]
  nd = pprintDec opts
  ni = pprintInt opts

-- | printf format string for a Complex
prepC :: (Epsilon t, Ord t) => PPrintOptions -> Complex t -> String
prepC opts (r :+ i) = prepD opts r ++ oi where
    oi | isNz i = concat [s, prepD opts i'," i"]
       | otherwise = []
    s | signum i >= 0 = " + "
      | otherwise = " - "
    i' = abs i


-- | printf for a list of values
--
-- > printDN (PPOpts 1 3) 2 [1,pi]
printDN opts xl0
  | n==0 = printf "\n"
  | n==1 = let [x1]=nums in printf strsc x1
  | n==2 = let [x1,x2]=nums in printf strsc x1 x2
  | n==3 = let [x1,x2,x3]=nums in printf strsc x1 x2 x3
  | n==4 = let [x1,x2,x3,x4]=nums in printf strsc x1 x2 x3 x4
  where
    nmax = 4
    (n, strs, nums) = printN opts prepD nmax xl0
    strsc = commas strs ++ "\n"
    -- pr = prepD opts

printN :: Epsilon a =>
     t -> (t -> a -> String) -> Int -> [a] -> (Int, [String], [a])
printN opts prepf nmax xl0 = go 0 xl0 [] [] where
  pr = prepf opts
  go i [x] ss ns | isNz x = go (i+1) [] (pr x : ss) (x : ns)
                 | otherwise = go i [] ("_" : ss) ns  
  go i (x:xs) ss ns | isNz x = if i < nmax - 2
                               then go (i+1) xs (pr x : ss) (x : ns)
                               else if i == nmax - 2
                                    then go (i+1) xs (" ... " : pr x : ss) (x : ns)
                                    else go i xs ss ns
                    | otherwise = if i < nmax - 2 
                                  then go i xs ("_" : ss) ns
                                  else if i == nmax - 2 
                                       then go (i+1) xs (" ... " : "_" : ss) ns
                                       else go i xs ss ns
  go nfin [] ss ns = (nfin, reverse ss , reverse ns)


-- | printf for list of complex values
-- 
-- > printCN pdef [(1:+pi), (3.5:+4.3), (pi:+(-3.4))]
printCN
  :: (PrintfArg t1, PrintfType t, Epsilon t1, Ord t1) =>
     PPrintOptions -> [Complex t1] -> t
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
