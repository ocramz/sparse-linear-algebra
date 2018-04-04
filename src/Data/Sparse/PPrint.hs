{-# language FlexibleContexts, RecordWildCards #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
--
-- Pretty printing helper functions
-----------------------------------------------------------------------------
module Data.Sparse.PPrint (prd, prd0, PrintDense, newline
                           , PPrintOptions(..)
                             -- , prdefR, prdefC
                           , printDN , printCN
                          ) where

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


-- | Pretty printing options: total length in # digits (including the decimal point), # of decimal digits
data PPrintOptions =
  PPOpts {
     pprintLen, pprintDec, pprintColWidth :: Int
     } deriving (Eq, Show)

-- | Some defaults
prdefR, prdefC :: PPrintOptions
prdefR = PPOpts 4 2 7   -- reals
prdefC = PPOpts 4 2 16  -- complex values


-- | Pretty print an array of real numbers
printDN :: (PrintfArg a, Epsilon a, Ord a, Floating a) =>
     Int -> Int -> PPrintOptions -> [a] -> String
printDN l n = printNpad l n f where
  f o x | isNz x = printf (prepD o x) x
        | otherwise = printf "_"


-- | Pretty print an array of complex numbers
printCN :: (PrintfArg a, Floating a, Epsilon a, Epsilon (Complex a), Ord a) =>
     Int -> Int -> PPrintOptions -> [Complex a] -> String
printCN l n = printNpad l n f where
  f o x | nearZero (re x) && isNz (imagPart x) =
               printf (prepD o (imagPart x)++ "i") (aim x)
        | nearZero (imagPart x) && isNz (re x) =
               printf (prepD o (realPart x)) (re x)
        | isNz x = printf (prepC o x) (re x) (aim x)
        | otherwise = printf "_"






-- | printf an array of items with padding space to render a fixed column width
printNpad ::
     Int     -- ^ Length of list to be provided
     -> Int  -- ^ 
     -> (PPrintOptions -> a -> String)
     -> PPrintOptions -> [a] -> String
printNpad llen nmax f o@PPOpts{..} xxl = commas [h,l] where
  h = commas $ take hlen ll
  l = last ll
  hlen = min (llen-1) (nmax-1)
  ll = unfoldr g (0, xxl) 
  g (i, x:xs) | i<nmax-2 || llen>=nmax-1 = Just (s', sxs)
              | i==nmax-2 = Just (dots', sxs)
              | null xs = Just (s', sxs)
              | otherwise = Just ("", sxs) where                  
                  s = f o x
                  sxs = (succ i, xs)
                  s' = s ++ spaces (n - length s) 
                  dots' = dots ++ spaces (n - length dots)
  g (_, []) = Nothing
  n = pprintColWidth
  dots = " ... "






-- | printf format string
prepD :: (Ord t, Floating t, Epsilon t) => PPrintOptions -> t -> String
prepD PPOpts{..} x = s where
  s | nearZero x  = "_"
    | abs x > magHi || abs x < magLo = s0 ++ "e"
    | otherwise = s0 ++ "f"
  s0 = concat ["%" , show nl, ".", show nd]
  -- s0 = "%1." ++ show nd
  nl = pprintLen  
  nd = pprintDec 
  nint = nl - nd  -- # of integer digits
  magLo = 10 ** (- fromIntegral nd)
  magHi = 10 ** fromIntegral nint
            


-- | printf format string for a Complex
prepC :: (Epsilon t, Floating t, Ord t) => PPrintOptions -> Complex t -> String
prepC opts (r :+ i) = prepD opts r ++ oi where
  oi = concat [s, prepD opts i', "i"]
  s | signum i >= 0 = " + "
    | otherwise = " - "
  i' = abs i




    

-- | Helpers

spaces :: Int -> String
spaces n = replicate n ' '

commas :: [String] -> String    
commas = intercalate ", "

re :: Complex a -> a
re = realPart

aim :: Num a => Complex a -> a
aim = abs . imagPart
