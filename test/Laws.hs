{-# LANGUAGE PartialTypeSignatures #-}
module Laws where

import Test.Hspec
import Numeric.LinearAlgebra.Class
import Numeric.Eps

infix 1 `shouldBeNear`

shouldBeNear :: (Epsilon a, _) => a -> a -> _
shouldBeNear x y = x ^-^ y `shouldSatisfy` nearZero --x `shouldSatisfy` (near y)

nearAssociative :: (Epsilon a, AdditiveGroup a, _) => a -> a -> a -> _
nearAssociative a b c = (a ^+^ b) ^+^ c `shouldBeNear` a ^+^ (b ^+^ c)

nearCancellative :: (Epsilon a, AdditiveGroup a, _) => a -> _
nearCancellative a = a ^-^ a `shouldBeNear` zeroV

commutative :: (Eq a, AdditiveGroup a, _) => a -> a -> _
commutative a b = a ^+^ b `shouldBe` b ^+^ a

neutralZero :: (Eq a, AdditiveGroup a, _) => a -> _
neutralZero a = a ^+^ zeroV `shouldBe` a

nearAssociativeScalar :: (Epsilon v, VectorSpace v, _) => Scalar v -> Scalar v -> v -> _
nearAssociativeScalar a b v = (a * b) .* v `shouldBeNear` a .* (b .* v)

neutralScalar :: (Eq v, VectorSpace v, _) => v -> _
neutralScalar v = 1 .* v `shouldBe` v

nearDistributiveScalar :: (Epsilon v, VectorSpace v, _) => Scalar v -> Scalar v -> v -> _
nearDistributiveScalar a b v = (a + b) .* v `shouldBeNear` a .* v ^+^ b .* v

nearDistributiveScalar2 :: (Epsilon v, VectorSpace v, _) => Scalar v -> v -> v -> _
nearDistributiveScalar2 a v w = a .* (v ^+^ w) `shouldBeNear` a .* v ^+^ a .* w
