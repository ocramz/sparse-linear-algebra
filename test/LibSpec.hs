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
    it "BiCGSTAB" $ 
      normSq (_xBicgstab (bicgstabN aa0 b0 x0 x0 niter) ^-^ x0true) <= eps `shouldBe` True
    it "CGS" $ 
      normSq (_x (cgsN aa0 b0 x0 x0 niter) ^-^ x0true) <= eps `shouldBe` True



(m,n) = (2,2)

aa0 :: SpMatrix Double
aa0 = SM (m,n) im where
  row0 = mkSpVectorD n [1,2]
  row1 = mkSpVectorD n [3,4]
  im = IM.fromList [(0, row0), (1, row1)]

b0, x0 :: SpVector Double
b0 = mkSpVectorD m [8,18]

x0 = mkSpVectorD m [0.3,1.4]
-- r0hat = mkSpVectorD m [1.1, 0.9]

{-
[1 2] [2] = [8]
[3 4] [3]   [18]

-}

x0true = mkSpVectorD n [2,3]

