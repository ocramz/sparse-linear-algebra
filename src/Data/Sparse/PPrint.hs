{-# language FlexibleContexts, RecordWildCards #-}
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
                           PPrintOptions, prdef, prepD, prepC
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
prdef :: PPrintOptions
prdef = PPOpts 4 2 15



-- | Pretty print an array of real numbers
printDN :: (PrintfArg a, Epsilon a, Foldable f, Ord a) =>
     PPrintOptions -> f a -> String
printDN opts = printNpad f opts where
  f o x | isNz x = printf (prepD o x) x
        | otherwise = printf "_"


-- | Pretty print an array of complex numbers
printCN :: (PrintfArg a, Epsilon a, Epsilon (Complex a), Foldable f, Ord a) =>
     PPrintOptions -> f (Complex a) -> String
printCN opts = printNpad f opts where
  f o x | nearZero (re x) = printf (prepD o (imagPart x)++ "i") (aim x)
        | nearZero (imagPart x) = printf (prepD o (realPart x)) (re x)
        | isNz x = printf (prepC o x) (re x) (aim x)
        | otherwise = printf "_"





-- | printf an array of items with padding space to render a fixed column width
printNpad :: Foldable f =>
     (PPrintOptions -> a -> String) -> PPrintOptions -> f a -> String
printNpad f o@PPOpts{..} xl = concat $ foldr ins [] xl where
  ins x acc = (s ++ spad) : acc where
    n = pprintColWidth
    s = f o x
    spad = spaces (n - length s)

    



-- | printf format string
prepD :: (Ord t, Epsilon t) => PPrintOptions -> t -> String
prepD PPOpts{..} x = s where
  s | nearZero x  = "_"
    | abs x > 999 || abs x < 0.1 = s0 ++ "e"
    | otherwise = s0 ++ "f"
  s0 = concat ["%" , show ni, ".", show nd]
  nd = pprintDec 
  ni = pprintLen
            


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
