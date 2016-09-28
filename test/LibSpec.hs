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
    it "matMat : matrix-matrix product" $
      (m1 `matMat` m2) `shouldBe` m1m2
    it "eye : identity matrix" $
      infoSM (eye 10) `shouldBe` SMInfo 10 0.1
  describe "Math.Linear.Sparse : Linear solvers" $ do    
    it "BiCGSTAB" $ 
      normSq (_xBicgstab (bicgstab aa0 b0 x0 x0) ^-^ x0true) <= eps `shouldBe` True
    it "CGS" $ 
      normSq (_x (cgs aa0 b0 x0 x0) ^-^ x0true) <= eps `shouldBe` True
  let n = 10
      nsp = 3
  describe ("random sparse linear system of size " ++ show n ++ " and sparsity " ++ show (fromIntegral nsp/fromIntegral n)) $ it "<\\>" $ do
    aa <- randSpMat n nsp
    xtrue <- randSpVec n nsp
    b <- randSpVec n nsp    
    let b = aa #> xtrue
    printDenseSM aa
    normSq (aa <\> b ^-^ xtrue) <= eps `shouldBe` True
  --     normSq (_xBicgstab (bicgstab aa b x0 x0) ^-^ x) <= eps `shouldBe` True



-- run N iterations 

-- runNBiC :: Int -> SpMatrix Double -> SpVector Double -> BICGSTAB
runNBiC n aa b = map _xBicgstab $ runAppendN' (bicgstabStep aa x0) n bicgsInit where
   x0 = mkSpVectorD nd $ replicate nd 0.9
   nd = dimSV r0
   r0 = b ^-^ (aa #> x0)    
   p0 = r0
   bicgsInit = BICGSTAB x0 r0 p0

-- runNCGS :: Int -> SpMatrix Double -> SpVector Double -> CGS
runNCGS n aa b = map _x $ runAppendN' (cgsStep aa x0) n cgsInit where
  x0 = mkSpVectorD nd $ replicate nd 0.1
  nd = dimSV r0
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = r0
  u0 = r0
  cgsInit = CGS x0 r0 p0 u0  


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



-- dense
solveRandom n = do
  aa0 <- randMat n
  let aa = aa0 ^+^ eye n
  xtrue <- randVec n
  -- x0 <- randVec n
  let b = aa #> xtrue
      dx = aa <\> b ^-^ xtrue
  return $ normSq dx
  -- let xhatB = _xBicgstab (bicgstab aa b x0 x0)
  --     xhatC = _x (cgs aa b x0 x0)
  -- return (aa, x, x0, b, xhatB, xhatC)

-- sparse
solveSpRandom :: Int -> Int -> IO Double
solveSpRandom n nsp = do
  aa0 <- randSpMat n nsp
  let aa = aa0 ^+^ eye n
  xtrue <- randSpVec n nsp
  let b = (aa ^+^ eye n) #> xtrue
      dx = aa <\> b ^-^ xtrue
  return $ normSq dx



-- `ndim` iterations

solveRandomN ndim nsp niter = do
  aa0 <- randSpMat ndim (nsp ^ 2)
  let aa = aa0 ^+^ eye ndim
  xtrue <- randSpVec ndim nsp
  let b = aa #> xtrue
      xhatB = head $ runNBiC niter aa b
      xhatC = head $ runNCGS niter aa b
  printDenseSM aa    
  return (normSq (xhatB ^-^ xtrue), normSq (xhatC ^-^ xtrue))

--

{-
matMat

[1, 2] [5, 6] = [19, 22]
[3, 4] [7, 8]   [43, 50]
-}

m1 = fromListDenseSM 2 2 [1,3,2,4]
m2 = fromListDenseSM 2 2 [5, 7, 6, 8]
m1m2 = fromListDenseSM 2 2 [19, 43, 22, 50]
