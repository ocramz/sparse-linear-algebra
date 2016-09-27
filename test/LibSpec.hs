module LibSpec where

import qualified Data.IntMap as IM

import Control.Monad (replicateM)

import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWC
       
import Test.Hspec
-- import Test.Hspec.QuickCheck

import Lib
import Math.Linear.Sparse



main :: IO ()
main = hspec spec

-- niter = 5

spec :: Spec
spec = do
  describe "Math.Linear.Sparse : library" $ do
    -- it "works" $ do
    --   True `shouldBe` True
    -- prop "ourAdd is commutative" $ \x y ->
    --   ourAdd x y `shouldBe` ourAdd y x
    it "matVec : matrix-vector product" $
      normSq ((aa0 #> x0true) ^-^ b0 ) <= eps `shouldBe` True
    it "eye : identity matrix" $
      infoSM (eye 10) `shouldBe` SMInfo 10 0.1
  describe "Math.Linear.Sparse : Linear solvers" $ do    
    it "BiCGSTAB" $ 
      normSq (_xBicgstab (bicgstab aa0 b0 x0 x0) ^-^ x0true) <= eps `shouldBe` True
    it "CGS" $ 
      normSq (_x (cgs aa0 b0 x0 x0) ^-^ x0true) <= eps `shouldBe` True


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




-- --

{-
example 1 : random linear system

-}


randMat1 n = do
  g <- MWC.create
  aav <- replicateM (n^2) (MWC.normal 0 1 g)
  let ii_ = [0 .. n-1]
      (ix_,iy_) = unzip $ concatMap (zip ii_ . replicate n) ii_
  return $ fromListSM (n,n) $ zip3 ix_ iy_ aav
  

randSol1 n = do
  g <- MWC.create
  bv <- replicateM n (MWC.normal 0 1 g)
  let ii_ = [0..n-1]
  return $ fromListSV n $ zip ii_ bv

randRhs n = do
  aa <- randMat1 n
  x <- randSol1 n
  return $ aa #> x
