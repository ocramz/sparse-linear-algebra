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
                           , PPrintOptions(..), prdefR, prdefC
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
prdefR = PPOpts 6 2 7   -- reals
prdefC = PPOpts 6 2 16  -- complex values


-- | Pretty print an array of real numbers
printDN :: (PrintfArg a, Epsilon a, Ord a) =>
     Int -> PPrintOptions -> [a] -> String
printDN n = printNpad n f where
  f o x | isNz x = printf (prepD o x) x
        | otherwise = printf "_"


-- | Pretty print an array of complex numbers
printCN :: (PrintfArg a, Epsilon a, Epsilon (Complex a), Ord a) =>
     Int -> PPrintOptions -> [Complex a] -> String
printCN n = printNpad n f where
  f o x | nearZero (re x) = printf (prepD o (imagPart x)++ "i") (aim x)
        | nearZero (imagPart x) = printf (prepD o (realPart x)) (re x)
        | isNz x = printf (prepC o x) (re x) (aim x)
        | otherwise = printf "_"






-- | printf an array of items with padding space to render a fixed column width
printNpad ::
     Int -> (PPrintOptions -> a -> String) -> PPrintOptions -> [a] -> String
printNpad nmax f o@PPOpts{..} xxl = commas [h,l] where
  h = commas $ take (nmax-1) ll
  l = last ll
  ll = unfoldr g (0, xxl) 
  g (i, x:xs) | i<nmax-2 = Just (s', (succ i, xs))
              | i==nmax-2 = Just (dots', (succ i, xs))
              | null xs = Just (s', (succ i, []))
              | otherwise = Just (spaces n, (succ i, xs)) where                  
                  s = f o x
                  s' = s ++ spad
                  dots' = dots ++ spadd
                  spad = spaces (n - length s)
                  spadd = spaces (n - length dots)
  g (_, []) = Nothing
  n = pprintColWidth
  dots = " ... "



-- | printf format string
prepD :: (Ord t, Epsilon t) => PPrintOptions -> t -> String
prepD PPOpts{..} x = s where
  s | nearZero x  = "_"
    | abs x > magHi-1 || abs x < magLo = s0 ++ "e"
    | otherwise = s0 ++ "f"
  s0 = concat ["%" , show nl, ".", show nd]
  -- s0 = "%1." ++ show nd
  nl = pprintLen  
  nd = pprintDec 
  nint = nl - nd - 1 -- # of integer digits
  magLo = 10 ** (- fromIntegral nd)
  magHi = 10 ** fromIntegral nint
            


-- | printf format string for a Complex
prepC :: (Epsilon t, Ord t) => PPrintOptions -> Complex t -> String
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
