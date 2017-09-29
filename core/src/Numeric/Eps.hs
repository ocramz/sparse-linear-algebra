{-# language FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016-2017 Marco Zocca, 2012-2015 Edward Kmett
-- License     :  BSD3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-- Testing for values "near" zero
-----------------------------------------------------------------------------
module Numeric.Eps
  ( Epsilon(..), isNz, roundZero, roundOne, roundZeroOne,
    nearOne
  ) where
import Data.Complex
import Foreign.C.Types (CFloat, CDouble)

-- | Provides a test to see if a quantity is near zero.
--
-- >>> nearZero (1e-11 :: Double)
-- False
--
-- >>> nearZero (1e-17 :: Double)
-- True
--
-- >>> nearZero (1e-5 :: Float)
-- False
--
-- >>> nearZero (1e-7 :: Float)
-- True
class (Floating a, Num a) => Epsilon a where
  -- | Determine if a quantity is near zero.
  nearZero :: a -> Bool

-- | @'abs' a '<=' 1e-6@
instance Epsilon Float where
  nearZero a = abs a <= 1e-6

-- | @'abs' a '<=' 1e-12@
instance Epsilon Double where
  nearZero a = abs a <= 1e-12

-- | @'abs' a '<=' 1e-6@
instance Epsilon CFloat where
  nearZero a = abs a <= 1e-6

-- | @'abs' a '<=' 1e-12@
instance Epsilon CDouble where
  nearZero a = abs a <= 1e-12


-- | @'magnitude' a '<=' 1e-6@
instance Epsilon (Complex Float) where
  nearZero a = magnitude a <= 1e-6

-- | @'magnitude' a '<=' 1e-12@
instance Epsilon (Complex Double) where
  nearZero a = magnitude a <= 1e-12

-- | @'magnitude' a '<=' 1e-6@
instance Epsilon (Complex CFloat) where
  nearZero a = magnitude a <= 1e-6

-- | @'magnitude' a '<=' 1e-12@
instance Epsilon (Complex CDouble) where
  nearZero a = magnitude a <= 1e-12



-- * Rounding operations


-- | Is this quantity close to 1 ?
nearOne :: Epsilon a => a -> Bool
nearOne x = nearZero (1 - x)

-- | Is this quantity distinguishable from 0 ?
isNz :: Epsilon a => a -> Bool
isNz x = not (nearZero x)

withDefault :: (t -> Bool) -> t -> t -> t
withDefault q d x | q x = d
                  | otherwise = x

roundZero, roundOne, roundZeroOne :: Epsilon a => a -> a
roundZero = withDefault nearZero (fromIntegral (0 :: Int))
roundOne = withDefault nearOne (fromIntegral (1 :: Int))

with2Defaults :: (t -> Bool) -> (t -> Bool) -> t -> t -> t -> t
with2Defaults q1 q2 d1 d2 x | q1 x = d1
                            | q2 x = d2
                            | otherwise = x

-- | Round to respectively 0 or 1
roundZeroOne = with2Defaults nearZero nearOne (fromIntegral (0 :: Int)) (fromIntegral (1 :: Int))
