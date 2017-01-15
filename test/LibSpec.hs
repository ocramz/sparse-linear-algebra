{-# LANGUAGE FlexibleContexts, TypeFamilies #-}
{-# language ScopedTypeVariables, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Copyright   :  (C) 2016 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------
module LibSpec where

import Data.Sparse.Common
import Numeric.LinearAlgebra.Sparse
-- import Numeric.LinearAlgebra.Class

-- import Control.Applicative (liftA2)
-- -- import Control.Monad (liftM, liftM2, replicateM)
-- import Control.Monad.Primitive
-- import Data.Foldable (foldrM)

import Data.Complex
import Data.Either (either, isRight)

import Data.VectorSpace hiding (magnitude)

import Control.Monad.State.Strict (execState)

-- import qualified System.Random.MWC as MWC
-- import qualified System.Random.MWC.Distributions as MWC
       
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck




main :: IO ()
main = hspec spec

-- niter = 5

spec :: Spec
spec = do
  describe "Numeric.LinearAlgebra.Sparse : Library" $ do
    prop "Subtraction is cancellative" $ \(x :: SpVector Double) ->
      norm2Sq (x ^-^ x) `shouldBe` zero
    it "<.> : inner product (Real)" $
      tv0 <.> tv0 `shouldBe` 61
    it "<.> : inner product (Complex)" $
      tvc2 <.> tvc3 `shouldBe` 2 :+ (-2)  
    it "transposeSM : sparse matrix transpose" $
      transposeSM m1 `shouldBe` m1t
    it "(#>) : matrix-vector product (Real)" $
      nearZero ( norm2Sq ((aa0 #> x0true) ^-^ b0 )) `shouldBe` True
    it "(<#) : vector-matrix product (Real)" $
      nearZero ( norm2Sq ((x0true <# aa0) ^-^ aa0tx0 ))`shouldBe` True  
    it "(##) : matrix-matrix product (Real, square)" $ 
      (m1 ## m2) `shouldBe` m1m2
    it "(##) : matrix-matrix product (Real, rectangular)" $ do
      (m1' ## m2') `shouldBe` m1m2'
      (m2' ## m1') `shouldBe` m2m1'
    it "(##) : matrix-matrix product (Complex)" $
      (aa3c ## aa3c) `shouldBe` aa3cx 
    it "eye : identity matrix" $
      infoSM (eye 10) `shouldBe` SMInfo 10 0.1
    it "insertCol : insert a column in a SpMatrix" $
      insertCol (eye 3) (fromListDenseSV 3 [2,2,2]) 0 `shouldBe` fromListSM (3,3) [(0,0,2),(1,0,2),(1,1,1),(2,0,2),(2,2,1)]
    it "insertRow : insert a row in a SpMatrix" $
      insertRow (eye 3) (fromListDenseSV 3 [2,2,2]) 1 `shouldBe` fromListSM (3,3) [(0,0,1), (1,0,2), (1,1,2), (1,2,2), (2,2,1)]
    it "extractCol -> insertCol : identity" $
      insertCol (eye 3) (extractCol (eye 3) 1) 1 `shouldBe` eye 3
    it "extractRow -> insertRow : identity" $
      insertRow (eye 3) (extractRow (eye 3) 1) 1 `shouldBe` eye 3      
    it "countSubdiagonalNZ : # of nonzero elements below the diagonal" $
      countSubdiagonalNZSM m3 `shouldBe` 1
    it "permutPairsSM : permutation matrices are orthogonal" $ do
      let pm0 = permutPairsSM 3 [(0,2), (1,2)] :: SpMatrix Double
      pm0 #~#^ pm0 `shouldBe` eye 3
      pm0 #~^# pm0 `shouldBe` eye 3
    it "isLowerTriSM : checks whether matrix is lower triangular" $
      isLowerTriSM tm8' && isUpperTriSM tm8 `shouldBe` True
    it "modifyInspectN : early termination by iteration count" $
      execState (modifyInspectN 2 (nearZero . diffSqL) (/2)) (1 :: Double) `shouldBe` 1/8
    it "modifyInspectN : termination by value convergence" $
      nearZero (execState (modifyInspectN (2^16) (nearZero . head) (/2)) (1 :: Double)) `shouldBe` True
    -- prop "aa2 is positive semidefinite" $ \(v :: SpVector Double) ->
    --   prop_psd aa2 v
  describe "QuickCheck properties:" $ do
    prop "prop_spd : (m #^# m) is symmetric positive definite" $
      \(PropSPD (m :: SpMatrix Double) v) -> prop_spd m v
    prop "prop_dot : (v <.> v) ~= 1 if ||v|| == 1" $
      \(v :: SpVector Double) -> prop_dot v
  describe "Numeric.LinearAlgebra.Sparse : Iterative linear solvers (Real)" $ do
    -- it "TFQMR (2 x 2 dense)" $
    it "GMRES (2 x 2 dense)" $
      checkLinSolveR GMRES_ aa0 b0 x0true `shouldBe` True
    it "GMRES (3 x 3 sparse, symmetric pos.def.)" $
      checkLinSolveR GMRES_ aa2 b2 x2 `shouldBe` True
    it "GMRES (4 x 4 sparse)" $
      checkLinSolveR GMRES_ aa1 b1 x1 `shouldBe` True
    it "BCG (2 x 2 dense)" $
      checkLinSolveR BCG_ aa0 b0 x0true `shouldBe` True
    it "BCG (3 x 3 sparse, symmetric pos.def.)" $
      checkLinSolveR BCG_ aa2 b2 x2 `shouldBe` True
    -- it "BiCGSTAB (2 x 2 dense)" $ 
    --   nearZero (normSq (linSolve BICGSTAB_ aa0 b0 ^-^ x0true)) `shouldBe` True
    it "BiCGSTAB (3 x 3 sparse, symmetric pos.def.)" $ 
      checkLinSolveR BICGSTAB_ aa2 b2 x2 `shouldBe` True
    it "CGS (2 x 2 dense)" $ 
      checkLinSolveR CGS_ aa0 b0 x0true `shouldBe` True
    it "CGS (3 x 3 sparse, SPD)" $ 
      checkLinSolveR CGS_ aa2 b2 x2 `shouldBe` True
  describe "Numeric.LinearAlgebra.Sparse : Direct linear solvers (Real)" $ 
    it "luSolve (4 x 4 sparse)" $ 
      checkLuSolve aa1 b1 `shouldBe` True         
  describe "Numeric.LinearAlgebra.Sparse : QR factorization (Real)" $ do    
    it "qr (4 x 4 sparse)" $
      checkQr tm4 `shouldBe` True
    it "qr (3 x 3 dense)" $ 
      checkQr tm2 `shouldBe` True
    it "qr (10 x 10 sparse)" $
      checkQr tm7 `shouldBe` True
  describe "Numeric.LinearAlgebra.Sparse : QR factorization (Complex)" $ do
    it "qr (2 x 2 dense)" $
      checkQr aa3cx `shouldBe` True
  describe "Numeric.LinearAlgebra.Sparse : LU factorization (Real)" $ do
    it "lu (4 x 4 dense)" $
      checkLu tm6 `shouldBe` True
    it "lu (10 x 10 sparse)" $
      checkLu tm7 `shouldBe` True
  describe "Numeric.LinearAlgebra.Sparse : Cholesky factorization (Real, symmetric pos.def.)" $ 
    it "chol (5 x 5 sparse)" $
      checkChol tm7 `shouldBe` True
  describe "Numeric.LinearAlgebra.Sparse : Arnoldi iteration, early breakdown detection (Real)" $ do      
    it "arnoldi (4 x 4 dense)" $
      checkArnoldi tm6 4 `shouldBe` True
    it "arnoldi (5 x 5 sparse)" $
      checkArnoldi tm7 5 `shouldBe` True    


{- linear systems -}

checkLinSolve method aa b x x0r =
  either
    (error . show)
    (\xhat -> nearZero (norm2Sq (x ^-^ xhat)))
    (linSolve0 method aa b x0r)

checkLinSolveR
  :: LinSolveMethod ->
     SpMatrix Double ->       -- ^ operator
     SpVector Double ->       -- ^ r.h.s
     SpVector Double -> Bool  -- ^ candidate solution
checkLinSolveR method aa b x = checkLinSolve method aa b x x0r where
  x0r = mkSpVR n $ replicate n 0.1
  n = ncols aa

checkLinSolveC
  :: LinSolveMethod
     -> SpMatrix (Complex Double)
     -> SpVector (Complex Double)
     -> SpVector (Complex Double)
     -> Bool
checkLinSolveC method aa b x = checkLinSolve method aa b x x0r where
  x0r = mkSpVC n $ replicate n (0.1 :+ 0.1)
  n = ncols aa



{- Givens rotation-}
checkGivens1 :: (Elt a, MatrixRing (SpMatrix a), Epsilon a, Floating a) =>
     SpMatrix a -> IxRow -> IxCol -> (a, Bool)
checkGivens1 tm i j = (rij, nearZero rij) where
  g = givens tm i j
  r = g ## tm
  rij = r @@ (i, j)


{- QR-}

checkQr :: (Elt a, MatrixRing (SpMatrix a), Epsilon (MatrixNorm (SpMatrix a)), Epsilon a, Floating a) =>
     SpMatrix a -> Bool
checkQr a = c1 && c2 && c3 where
  (q, r) = qr a
  c1 = nearZero $ normFrobenius ((q #~# r) ^-^ a)
  c2 = isOrthogonalSM q
  c3 = isUpperTriSM r



{- LU -}
checkLu :: (Elt a, MatrixRing (SpMatrix a), VectorSpace (SpVector a), Epsilon (MatrixNorm (SpMatrix a)), Epsilon a) =>
     SpMatrix a -> Bool
checkLu a = c1 && c2 where
  (l, u) = lu a
  c1 = nearZero (normFrobenius ((l #~# u) ^-^ a))
  c2 = isUpperTriSM u && isLowerTriSM l



{- Cholesky -}

checkChol :: (Elt a, MatrixRing (SpMatrix a), Epsilon (MatrixNorm (SpMatrix a)), Epsilon a, Floating a) =>
     SpMatrix a -> Bool
checkChol a = c1 && c2 where
  l = chol a
  c1 = nearZero $ normFrobenius ((l ##^ l) ^-^ a)
  c2 = isLowerTriSM l


{- direct linear solver -}

checkLuSolve :: (Scalar (SpVector t) ~ t, MatrixType (SpVector t) ~ SpMatrix t,
      Elt t, Normed (SpVector t), LinearVectorSpace (SpVector t),
      Epsilon (Magnitude (SpVector t)), Epsilon t) =>
     SpMatrix t -> SpVector t -> Bool
checkLuSolve amat rhs = nearZero (norm2Sq ( (lmat #> (umat #> xlu)) ^-^ rhs ))
  where
     (lmat, umat) = lu amat
     xlu = luSolve lmat umat rhs
      
  
{- Arnoldi iteration -}
-- checkArnoldi :: (Epsilon a, Floating a, Eq a) => SpMatrix a -> Int -> Bool
checkArnoldi aa kn = nearZero (normFrobenius $ lhs ^-^ rhs) where
  b = onesSV (nrows aa)
  (q, h) = arnoldi aa b kn
  (m, n) = dim q
  q' = extractSubmatrix q (0, m - 1) (0, n - 2) -- q' = all but one column of q
  rhs = q #~# h
  lhs = aa #~# q'







-- * Arbitrary newtypes and instances for QuickCheck



-- | An Arbitrary SpVector such that at least one entry is nonzero
instance Arbitrary (SpVector Double) where
  arbitrary = sized genf `suchThat` any isNz where
    genf n = do
      v_ <- vector n
      return $ fromListDenseSV n v_-- SV n (IM.fromList (zip i_ v_))




-- instance QC.Arbitrary (SpMatrix Double) where
--   arbitrary = QC.sized genf `QC.suchThat` any isNz where
--     genf n = do
--       let i_ = [0 .. n - 1]
--       v_ <- QC.vector n
--       return $ SV n (IM.fromList (zip i_ v_))


sized2 f = sized $ \i -> sized $ \j -> f i j




-- | A square matrix and vector of compatible size
data PropMatVec a = PropMatVec (SpMatrix a) (SpVector a) deriving (Eq, Show)

instance Arbitrary (PropMatVec Double) where
  arbitrary = sized genf where
    genf n = do
      let d = 4 * n
      i_ <- vectorOf d (choose (0, n-1))
      j_ <- vectorOf d (choose (0, n-1))      
      x <- vector d
      let m = fromListSM (n,n) $ zip3 i_ j_ x
      iv <- vectorOf d (choose (0, n-1))
      xv <- vector d           
      let v = fromListSV n $ zip iv xv 
      return $ PropMatVec m v



-- | A symmetric positive definite matrix and vector of compatible size
data PropSPD a = PropSPD (SpMatrix a) (SpVector a) deriving (Eq, Show)

instance Arbitrary (PropSPD Double) where
  arbitrary = do
    PropMatVec m v <- arbitrary :: Gen (PropMatVec Double)
    return $ PropSPD (m #^# m) v


    








-- | QuickCheck properties

-- | Dot product of a normalized vector with itself is ~= 1
prop_dot :: (Normed v, Epsilon (Scalar v)) => v -> Bool
prop_dot v = let v' = normalize2 v in nearOne (v' <.> v')

-- | Positive semidefinite. 
prop_spd :: (LinearVectorSpace v, InnerSpace v, Ord (Scalar v), Num (Scalar v)) =>
     MatrixType v -> v -> Bool
prop_spd mm v = (v <.> (mm #> v)) >= 0

prop_spd' :: PropSPD Double -> Bool
prop_spd' (PropSPD m v) = prop_spd m v


-- -- test data




{-

example 0 : 2x2 linear system

[1 2] [2] = [8]
[3 4] [3]   [18]

[1 3] [2] = [11]
[2 4] [3]   [16]


-}


aa0 :: SpMatrix Double
aa0 = fromListDenseSM 2 [1,3,2,4]

-- b0, x0 : r.h.s and initial solution resp.
b0, x0, x0true, aa0tx0 :: SpVector Double
b0 = mkSpVR 2 [8,18]
x0 = mkSpVR 2 [0.3,1.4]


-- x0true : true solution
x0true = mkSpVR 2 [2,3]

aa0tx0 = mkSpVR 2 [11,16]







{- 4x4 system -}

aa1 :: SpMatrix Double
aa1 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]

x1, b1 :: SpVector Double
x1 = mkSpVR 4 [1,2,3,4]

b1 = mkSpVR 4 [30,56,60,101]



{- 3x3 system -}
aa2 :: SpMatrix Double
aa2 = sparsifySM $ fromListDenseSM 3 [2, -1, 0, -1, 2, -1, 0, -1, 2]
x2, b2 :: SpVector Double
x2 = mkSpVR 3 [3,2,3]

b2 = mkSpVR 3 [4,-2,4]


aa22 = fromListDenseSM 2 [2,1,1,2] :: SpMatrix Double





{- 2x2 Complex system -}

aa0c :: SpMatrix (Complex Double)
aa0c = fromListDenseSM 2 [ 3 :+ 1, (-3) :+ 2, (-2) :+ (-1), 1 :+ (-2)]

b0c = mkSpVC 2 [3 :+ (-4), (-1) :+ 0.5]

x1c = mkSpVC 2 [2 :+ 2, 2 :+ 3]
b1c = mkSpVC 2 [4 :+ (-2), (-10) :+ 1]

aa2c :: SpMatrix (Complex Double)
aa2c = fromListDenseSM 2 [3, -3, -2, 1]








-- matlab : aa = [1, 2-j; 2+j, 1-j]
aa3c, aa3cx :: SpMatrix (Complex Double)
aa3c = fromListDenseSM 2 [1, 2 :+ 1, 2 :+ (-1), 1 :+ (-1)]

-- matlab : aaxaa = aa * aa
aa3cx = fromListDenseSM 2 [6, 5, 3 :+ (-4), 5:+ (-2)]







{-
matMat

[1, 2] [5, 6] = [19, 22]
[3, 4] [7, 8]   [43, 50]
-}

m1, m2, m1m2, m1', m2', m1m2', m2m1' :: SpMatrix Double
m1 = fromListDenseSM 2 [1,3,2,4]
m2 = fromListDenseSM 2 [5, 7, 6, 8]     
m1m2 = fromListDenseSM 2 [19, 43, 22, 50]

m1' = fromListSM (2,3) [(0,0,2), (1,0,3), (1,2,4), (1,2,1)]
m2' = fromListSM (3,2) [(0,0,5), (0,1,3), (2,1,4)]
m1m2' = fromListDenseSM 2 [10,15,6,13] 
m2m1' = fromListSM (3,3) [(0,0,19),(2,0,12),(0,2,3),(2,2,4)]

-- transposeSM

m1t = fromListDenseSM 2 [1,2,3,4]


--

{-
countSubdiagonalNZ
-}

m3 = fromListSM (3,3) [(0,2,3),(2,0,4),(1,1,3)] 






{- eigenvalues -}


aa3 = fromListDenseSM 3 [1,1,3,2,2,2,3,1,1] :: SpMatrix Double

b3 = mkSpVR 3 [1,1,1] :: SpVector Double



-- aa4 : eigenvalues 1 (mult.=2) and -1
aa4 = fromListDenseSM 3 [3,2,-2,2,2,-1,6,5,-4] :: SpMatrix Double

aa4c = toC <$> aa4

b4 = fromListDenseSV 3 [-3,-3,-3] :: SpVector Double






tm0, tm1, tm2, tm3, tm4 :: SpMatrix Double
tm0 = fromListSM (2,2) [(0,0,pi), (1,0,sqrt 2), (0,1, exp 1), (1,1,sqrt 5)]

tv0, tv1 :: SpVector Double
tv0 = mkSpVR 2 [5, 6]

tv1 = fromListSV 2 [(0,1)] 

-- wikipedia test matrix for Givens rotation

tm1 = sparsifySM $ fromListDenseSM 3 [6,5,0,5,1,4,0,4,3]

-- tm1g1 = givens tm1 1 0
-- tm1a2 = tm1g1 ## tm1

-- tm1g2 = givens tm1a2 2 1
-- tm1a3 = tm1g2 ## tm1a2

-- tm1q = transposeSM (tm1g2 ## tm1g1)


-- wp test matrix for QR decomposition via Givens rotation

tm2 = fromListDenseSM 3 [12, 6, -4, -51, 167, 24, 4, -68, -41]


tm3 = transposeSM $ fromListDenseSM 3 [1 .. 9]

tm3g1 = fromListDenseSM 3 [1, 0,0, 0,c,-s, 0, s, c]
  where c= 0.4961
        s = 0.8682


--

tm4 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]



tm5 = fromListDenseSM 3 [2, -4, -4, -1, 6, -2, -2, 3, 8] :: SpMatrix Double


tm6 = fromListDenseSM 4 [1,3,4,2,2,5,2,10,3,6,8,11,4,7,9,12] :: SpMatrix Double

tm7 :: SpMatrix Double
tm7 = a ^+^ b ^+^ c where
  n = 5
  a = mkSubDiagonal n 1 $ replicate n (-1)
  b = mkSubDiagonal n 0 $ replicate n 2
  c = mkSubDiagonal n (-1) $ replicate n (-1)




tm8 :: SpMatrix Double
tm8 = fromListSM (2,2) [(0,0,1), (0,1,1), (1,1,1)]

tm8' :: SpMatrix Double
tm8' = fromListSM (2,2) [(0,0,1), (1,0,1), (1,1,1)]



tm9 :: SpMatrix Double
tm9 = fromListSM (4, 3) [(0,0,pi), (1,1, 3), (2,2,4), (3,2, 1), (3,1, 5)]





-- tvc0 <.> tvc1 = 5 
tvc0, tvc1, tvc2, tvc3 :: SpVector (Complex Double)
tvc0 = fromListSV 2 [(0,0), (1,2 :+ 1)]
tvc1 = fromListSV 2 [(0,0), (1, 2 :+ (-1))] 


-- dot([1+i, 2-i], [3-2i, 1+i]) = 2 - 2i
tvc2 = fromListDenseSV 2 [1 :+ 1,  2 :+ (-1)]
tvc3 = fromListDenseSV 2 [3 :+ (-2), 1 :+ 1 ]






-- l0 = [1,2,4,5,8]
-- l1 = [2,3,6]
-- l2 = [7]

-- v0,v1 :: V.Vector Int
-- v0 = V.fromList [0,1,2,5,6]
-- v1 = V.fromList [0,3,4,6]

-- -- e1, e2 :: V.Vector (Int, Double)
-- -- e1 = V.indexed $ V.fromList [1,0,0]
-- -- e2 = V.indexed $ V.fromList [0,1,0]

-- e1, e2:: CsrVector Double
-- e1 = fromListCV 4 [(0, 1)] 
-- e2 = fromListCV 4 [(1, 1)]
-- e3 = fromListCV 4 [(0, 1 :+ 2)] :: CsrVector (Complex Double)

-- e1c = V.indexed $ V.fromList [1,0,0] :: V.Vector (Int, Complex Double)

-- m0,m1,m2,m3 :: CsrMatrix Double
-- m0 = toCSR 2 2 $ V.fromList [(0,0, pi), (1,0,3), (1,1,2)]
-- m1 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (1,0,2), (1,1,3), (2,0,4), (2,3,1), (3,2,2)]
-- m2 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (2,0,4), (2,3,1), (3,2,2)]
-- m3 = toCSR 4 4 $ V.fromList [(1,0,5), (1,1,8), (2,2,3), (3,1,6)]







-- --


-- -- run N iterations 

-- -- runNBiC :: Int -> SpMatrix Double -> SpVector Double -> BICGSTAB
-- runNBiC n aa b = map _xBicgstab $ runAppendN' (bicgstabStep aa x0) n bicgsInit where
--    x0 = mkSpVectorD nd $ replicate nd 0.9
--    nd = dim r0
--    r0 = b ^-^ (aa #> x0)    
--    p0 = r0
--    bicgsInit = BICGSTAB x0 r0 p0

-- -- runNCGS :: Int -> SpMatrix Double -> SpVector Double -> CGS
-- runNCGS n aa b = map _x $ runAppendN' (cgsStep aa x0) n cgsInit where
--   x0 = mkSpVectorD nd $ replicate nd 0.1
--   nd = dim r0
--   r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
--   p0 = r0
--   u0 = r0
--   cgsInit = CGS x0 r0 p0 u0



-- solveRandomN ndim nsp niter = do
--   aa0 <- randSpMat ndim (nsp ^ 2)
--   let aa = aa0 ^+^ eye ndim
--   xtrue <- randSpVec ndim nsp
--   let b = aa #> xtrue
--       xhatB = head $ runNBiC niter aa b
--       xhatC = head $ runNCGS niter aa b
--   -- printDenseSM aa    
--   return (normSq (xhatB ^-^ xtrue), normSq (xhatC ^-^ xtrue))



{-
random linear system

-}



-- -- dense
-- solveRandom n = do
--   aa0 <- randMat n
--   let aa = aa0 ^+^ eye n
--   xtrue <- randVec n
--   -- x0 <- randVec n
--   let b = aa #> xtrue
--       dx = aa <\> b ^-^ xtrue
--   return $ normSq dx
--   -- let xhatB = _xBicgstab (bicgstab aa b x0 x0)
--   --     xhatC = _x (cgs aa b x0 x0)
--   -- return (aa, x, x0, b, xhatB, xhatC)

-- -- sparse
-- solveSpRandom :: Int -> Int -> IO Double
-- solveSpRandom n nsp = do
--   aa0 <- randSpMat n nsp
--   let aa = aa0 ^+^ eye n
--   xtrue <- randSpVec n nsp
--   let b = (aa ^+^ eye n) #> xtrue
--       dx = aa <\> b ^-^ xtrue
--   return $ normSq dx




-- solveRandomBanded n bw mu sig = do
--   let ndiags = 2*bw
--   bands <- replicateM (ndiags + 1) (randArray n mu sig)
--   xtrue <- randVec n
--   b <- randVec n
--   let
--     diags = [-bw .. bw - 1]

-- randDiagMat :: PrimMonad m =>
--      Rows -> Double -> Double -> Int -> m (SpMatrix Double)
-- randDiagMat n mu sig i = do
--   x <- randArray n mu sig
--   return $ mkSubDiagonal n i x
