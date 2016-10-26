module Numeric.Eps where

-- * Numerical tolerance for "near-0" tests
-- | eps = 1e-8 
-- eps :: Double
-- eps = fromRational $ toRational 1e-8

eps :: Fractional a => a       
eps = 0.00000001


-- almostZ x = abs (x - eps) == 0


-- * Rounding operations


-- | Rounding rule
almostZero, almostOne, isNz :: Real a => a -> Bool
almostZero x = abs (toRational x) <= eps
almostOne x = x' >= (1-eps) && x' < (1+eps) where x' = toRational x
isNz = not . almostZero

withDefault :: (t -> Bool) -> t -> t -> t
withDefault q d x | q x = d
                  | otherwise = x

roundZero, roundOne, roundZeroOne :: Real a => a -> a
roundZero = withDefault almostZero (fromIntegral 0)
roundOne = withDefault almostOne (fromIntegral 1)

with2Defaults :: (t -> Bool) -> (t -> Bool) -> t -> t -> t -> t
with2Defaults q1 q2 d1 d2 x | q1 x = d1
                            | q2 x = d2
                            | otherwise = x

-- | Round to respectively 0 or 1
roundZeroOne = with2Defaults almostZero almostOne (fromIntegral 0) (fromIntegral 1)
