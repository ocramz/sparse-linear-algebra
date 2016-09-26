module LibSpec where

import qualified Data.IntMap as IM
  
import Test.Hspec
import Test.Hspec.QuickCheck

import Lib
import Math.Linear.Sparse



main :: IO ()
main = hspec spec

niter = 5

spec :: Spec
spec =
  describe "Lib" $ do
    -- it "works" $ do
    --   True `shouldBe` True
    -- prop "ourAdd is commutative" $ \x y ->
    --   ourAdd x y `shouldBe` ourAdd y x
    it "matVec : matrix-vector product" $
      normSq ((aa0 #> x0true) ^-^ b0 ) <= eps `shouldBe` True
    it "BiCGSTAB" $ 
      normSq (_xBicgstab (bicgstabN aa0 b0 x0 x0 niter) ^-^ x0true) <= eps `shouldBe` True
    it "CGS" $ 
      normSq (_x (cgsN aa0 b0 x0 x0 niter) ^-^ x0true) <= eps `shouldBe` True


{-

example 0 : 2x2 linear system

[1 2] [2] = [8]
[3 4] [3]   [18]

-}

aa0 :: SpMatrix Double
aa0 = SM (2,2) im where
  im = IM.fromList [(0, aa0r0), (1, aa0r1)]

aa0r0, aa0r1 :: IM.IntMap Double
aa0r0 = IM.fromList [(0,1),(1,2)]
aa0r1 = IM.fromList [(0,3),(1,4)]


-- b0, x0 : r.h.s and initial solution resp.
b0, x0, x0true :: SpVector Double
b0 = mkSpVectorD 2 [8,18]
x0 = mkSpVectorD 2 [0.3,1.4]


-- x0true : true solution
x0true = mkSpVectorD 2 [2,3]

