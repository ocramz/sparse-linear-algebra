module Numeric.Eps where

-- * Numerical tolerance for "near-0" tests
-- | eps = 1e-8 
eps :: Double
eps = 1e-8




-- * Rounding operations

-- | Rounding rule
almostZero, almostOne :: Double -> Bool
almostZero x = abs x <= eps
almostOne x = x >= (1-eps) && x < (1+eps)

isNz :: Double -> Bool
isNz = not . almostZero

withDefault :: (t -> Bool) -> t -> t -> t
withDefault q d x | q x = d
                  | otherwise = x

roundZero, roundOne :: Double -> Double
roundZero = withDefault almostZero 0
roundOne = withDefault almostOne 1

with2Defaults :: (t -> Bool) -> (t -> Bool) -> t -> t -> t -> t
with2Defaults q1 q2 d1 d2 x | q1 x = d1
                            | q2 x = d2
                            | otherwise = x

-- | Round to respectively 0 or 1 within some predefined numerical precision eps
roundZeroOne :: Double -> Double
roundZeroOne = with2Defaults almostZero almostOne 0 1
