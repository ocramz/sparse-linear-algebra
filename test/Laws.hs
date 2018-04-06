{-# LANGUAGE PartialTypeSignatures, FlexibleContexts #-}
{-# GHC_OPTIONS -Wno-partial-type-signatures #-}
module Laws where

import Test.Hspec
import Numeric.LinearAlgebra.Class

-- All the axioms of an Abelian group
associative :: (Eq a, AdditiveGroup a, _) => a -> a -> a -> _
associative a b c = (a ^+^ b) ^+^ c `shouldBe` a ^+^ (b ^+^ c)

cancellative :: (Eq a, AdditiveGroup a, _) => a -> _
cancellative a = a ^-^ a `shouldBe` zeroV

commutative :: (Eq a, AdditiveGroup a, _) => a -> a -> _
commutative a b = a ^+^ b `shouldBe` b ^+^ a

neutralZero :: (Eq a, AdditiveGroup a, _) => a -> _
neutralZero a = a ^+^ zeroV `shouldBe` a

-- All the axioms of a module (action is associative and bilinear)
associativeScalar :: (Eq v, VectorSpace v, _) => Scalar v -> Scalar v -> v -> _
associativeScalar a b v = (a * b) .* v `shouldBe` a .* (b .* v)

neutralScalar :: (Eq v, VectorSpace v, _) => v -> _
neutralScalar v = 1 .* v `shouldBe` v

distributiveScalar :: (Eq v, VectorSpace v, _) => Scalar v -> Scalar v -> v -> _
distributiveScalar a b v = (a + b) .* v `shouldBe` a .* v ^+^ b .* v

distributiveScalar2 :: (Eq v, VectorSpace v, _) => Scalar v -> v -> v -> _
distributiveScalar2 a v w = a .* (v ^+^ w) `shouldBe` a .* v ^+^ a .* w

-- Inner product should be bilinear, commutative in the real case and positive
innerProductBilinear :: (Eq (Scalar v), InnerSpace v, _) => v -> v -> v -> _
innerProductBilinear a b c = (a ^+^ b) <.> c `shouldBe` (a <.> c) + (b <.> c)

innerProductBilinear2 :: (Eq (Scalar v), InnerSpace v, _) => Scalar v -> v -> v -> _
innerProductBilinear2 l b c = (l .* b) <.> c `shouldBe` l * (b <.> c)

innerProductCommutative :: (Eq (Scalar v), InnerSpace v, _) => v -> v -> _
innerProductCommutative v w = v <.> w `shouldBe` w <.> v

innerProductPositive :: (Ord (Scalar v), InnerSpace v, _) => v -> _
innerProductPositive v = (v <.> v) `shouldSatisfy` (>= 0)
