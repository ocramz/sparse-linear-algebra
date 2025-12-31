{-# LANGUAGE FlexibleContexts, TypeFamilies #-}
{-# language ScopedTypeVariables, FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -Wno-missing-signatures #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  MatrixFactorizationsSpec
-- Copyright   :  (C) 2016-2025 Marco Zocca
-- License     :  GPL-3 (see LICENSE)
-- Maintainer  :  zocca.marco gmail
-- Stability   :  provisional
-- Portability :  portable
--
-- Test suite for matrix factorization algorithms (QR, Cholesky)
-----------------------------------------------------------------------------
module MatrixFactorizationsSpec where

import Control.Exception.Common
import Data.Sparse.Common
import Numeric.LinearAlgebra.Sparse

import Control.Monad.Catch
import Control.Monad.IO.Class
import qualified Debug.Trace as DT

import Data.Complex
       
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Test.QuickCheck.Monadic (monadicIO, run, assert)


-- * Main Spec

spec :: Spec
spec = do
  specQR
  specChol


-- * QR Factorization Tests

specQR :: Spec
specQR = do       
  describe "Numeric.LinearAlgebra.Sparse : QR factorization (Real)" $ do
    it "qr (3 x 3 dense)" $ 
      checkQr0 qr tm2 >>= (`shouldBe` True)
    it "qr (4 x 4 sparse)" $
      checkQr0 qr tm4 >>= (`shouldBe` True)
    it "qr (4 x 4 dense)" $
      checkQr0 qr tm6 >>= (`shouldBe` True)
    it "qr (issue test case: 4x4 matrix)" $
      checkQr0 qr issueMatrix >>= (`shouldBe` True)
  describe "Numeric.LinearAlgebra.Sparse : QR factorization (Complex)" $ do
    it "qr (2 x 2 dense)" $
      checkQr0 qr aa3cx >>= (`shouldBe` True)
    it "qr (3 x 3 dense)" $
      checkQr0 qr tmc4 >>= (`shouldBe` True)

checkQr0 :: (Elt a, MatrixRing (SpMatrix a), Epsilon a, MonadThrow m) =>
     (SpMatrix a -> m (SpMatrix a, SpMatrix a))
     -> SpMatrix a
     -> m Bool
checkQr0 mfqr a = do
  (q, r) <- mfqr a
  let c1 = nearZero $ normFrobenius $ sparsifySM ((q ## r) ^-^ a)
      c2 = isOrthogonalSM q
      c3 = isUpperTriSM r
      upperPart = ifilterSM (\i j _ -> i <= j) r
      lowerPart = ifilterSM (\i j _ -> i > j) r
      result = c1 && c2 && c3
  return $ DT.trace ("\n=== checkQr0: reconstruction=" ++ show c1 ++ ", orthogonal=" ++ show c2 ++ 
                     ", upperTri=" ++ show c3 ++ ", R==upper=" ++ show (r == upperPart) ++
                     ", nnz(lower)=" ++ show (nnz lowerPart)) result

prop_QR :: (MonadThrow m) => PropMatI Double -> m Bool
prop_QR (PropMatI m) = checkQr0 qr m

prop_QR_complex :: (MonadThrow m) => PropMatIC (Complex Double) -> m Bool
prop_QR_complex (PropMatIC m) = checkQr0 qr m


-- * Cholesky Factorization Tests

specChol :: Spec
specChol = do
  describe "Numeric.LinearAlgebra.Sparse : Cholesky factorization properties (Real)" $ do
    it "chol: specific test case (3x3 SPD matrix)" $
      checkChol testSPD3x3 >>= (`shouldBe` True)
    it "chol: specific test case (5x5 sparse tridiagonal)" $
      checkChol tm7spd >>= (`shouldBe` True)
    prop "chol: L L^T = A for symmetric positive definite matrices (Real)" $
      \(PropMatSPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_reconstruction m
        assert result
    prop "chol: diagonal elements are positive (Real)" $
      \(PropMatSPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_positive_diagonal m
        assert result
  describe "Numeric.LinearAlgebra.Sparse : Cholesky factorization properties (Complex)" $ do
    it "chol: specific test case (2x2 HPD matrix)" $
      checkChol testHPD2x2 >>= (`shouldBe` True)
    it "chol: specific test case (3x3 HPD matrix)" $
      checkChol testHPD3x3 >>= (`shouldBe` True)
    prop "chol: L L^H = A for Hermitian positive definite matrices (Complex)" $
      \(PropMatHPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_reconstruction_complex m
        assert result
    prop "chol: diagonal elements have positive real parts (Complex)" $
      \(PropMatHPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_positive_diagonal_complex m
        assert result
  describe "Numeric.LinearAlgebra.Sparse : Cholesky arrowhead matrices (Real)" $ do
    it "chol: Rails example from gist (8x8 arrowhead)" $
      checkChol railsMatrix >>= (`shouldBe` True)
    prop "chol: L L^T = A for arrowhead SPD matrices (Real)" $
      \(PropMatArrowheadSPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_reconstruction m
        assert result
    prop "chol: L is lower triangular for arrowhead matrices (Real)" $
      \(PropMatArrowheadSPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_lower_triangular m
        assert result
    prop "chol: diagonal elements are positive for arrowhead matrices (Real)" $
      \(PropMatArrowheadSPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_positive_diagonal m
        assert result
  describe "Numeric.LinearAlgebra.Sparse : Cholesky arrowhead matrices (Complex)" $ do
    prop "chol: L L^H = A for arrowhead HPD matrices (Complex)" $
      \(PropMatArrowheadHPD m) -> monadicIO $ do
        result <- run $ prop_Cholesky_reconstruction_complex m
        assert result

checkChol :: (Elt a, MatrixRing (SpMatrix a), Epsilon a,
              InnerSpace a, Scalar a ~ a,
              MonadThrow m) =>
     SpMatrix a -> m Bool
checkChol a = do
  l <- chol a
  let c1 = nearZero $ normFrobenius ((l ##^ l) ^-^ a)
      c2 = isLowerTriSM l
  return $ c1 && c2

prop_Cholesky_reconstruction :: (MonadThrow m) =>
     SpMatrix Double -> m Bool
prop_Cholesky_reconstruction a = do
  l <- chol a
  let reconstruction = l ##^ l
      residual = normFrobenius (reconstruction ^-^ a)
  return $ nearZero residual

prop_Cholesky_lower_triangular :: (MonadThrow m) =>
     SpMatrix Double -> m Bool
prop_Cholesky_lower_triangular a = do
  l <- chol a
  return $ isLowerTriSM l

prop_Cholesky_positive_diagonal :: (MonadThrow m) =>
     SpMatrix Double -> m Bool
prop_Cholesky_positive_diagonal a = do
  l <- chol a
  let diag = extractDiagDense l
      allPositive = all (\(_, x) -> x > 0) (toListSV diag)
  return allPositive

prop_Cholesky_reconstruction_complex :: (MonadThrow m) =>
     SpMatrix (Complex Double) -> m Bool
prop_Cholesky_reconstruction_complex a = do
  l <- chol a
  let reconstruction = l ##^ l
      residual = normFrobenius (reconstruction ^-^ a)
  return $ nearZero residual

prop_Cholesky_lower_triangular_complex :: (MonadThrow m) =>
     SpMatrix (Complex Double) -> m Bool
prop_Cholesky_lower_triangular_complex a = do
  l <- chol a
  return $ isLowerTriSM l

prop_Cholesky_positive_diagonal_complex :: (MonadThrow m) =>
     SpMatrix (Complex Double) -> m Bool
prop_Cholesky_positive_diagonal_complex a = do
  l <- chol a
  let diag = extractDiagDense l
      allPositiveReal = all (\(_, x) -> realPart x > 0) (toListSV diag)
  return allPositiveReal


-- * QuickCheck Generators

-- | An arbitrary SpMatrix with identity diagonal 
newtype PropMatI a = PropMatI {unPropMatI :: SpMatrix a} deriving (Eq)
instance Show a => Show (PropMatI a) where show = show . unPropMatI
instance Arbitrary (PropMatI Double) where
  arbitrary = sized (\m -> PropMatI <$> genSpMI m) `suchThat` ((> 2) . nrows . unPropMatI)

-- | An arbitrary complex SpMatrix with identity diagonal 
newtype PropMatIC a = PropMatIC {unPropMatIC :: SpMatrix a} deriving (Eq)
instance Show a => Show (PropMatIC a) where show = show . unPropMatIC
instance Arbitrary (PropMatIC (Complex Double)) where
  arbitrary = sized (\m -> PropMatIC <$> genSpMI m) `suchThat` ((> 2) . nrows . unPropMatIC)

genSpMI :: (AdditiveGroup a, Num a, Arbitrary a) => Int -> Gen (SpMatrix a)
genSpMI m = do
  mm <- genSpM m m
  return $ mm ^+^ eye m

-- | A symmetric, positive-definite matrix (Real)
newtype PropMatSPD a = PropMatSPD {unPropMatSPD :: SpMatrix a} deriving (Eq, Show)

instance Arbitrary (PropMatSPD Double) where
  arbitrary = sized (\m -> PropMatSPD <$> genSpM_SPD m) `suchThat` ((>= 2) . nrows . unPropMatSPD)

genSpM_SPD :: Int -> Gen (SpMatrix Double)
genSpM_SPD n = do
  m <- genSpMDense n n
  let spd = (m #^# m) ^+^ (1.0 .* eye n)
  return spd

-- | A Hermitian positive-definite matrix (Complex)
newtype PropMatHPD a = PropMatHPD {unPropMatHPD :: SpMatrix a} deriving (Eq, Show)

instance Arbitrary (PropMatHPD (Complex Double)) where
  arbitrary = sized (\m -> PropMatHPD <$> genSpM_HPD m) `suchThat` ((>= 2) . nrows . unPropMatHPD)

genSpM_HPD :: Int -> Gen (SpMatrix (Complex Double))
genSpM_HPD n = do
  m <- genSpMDense n n
  let hpd = (m ##^ m) ^+^ ((1.0 :+ 0.0) .* eye n)
  return hpd

-- | Arrowhead SPD matrix (Real)
newtype PropMatArrowheadSPD a = PropMatArrowheadSPD {unPropMatArrowheadSPD :: SpMatrix a} deriving (Eq, Show)

instance Arbitrary (PropMatArrowheadSPD Double) where
  arbitrary = sized (\m -> PropMatArrowheadSPD <$> genSpM_ArrowheadSPD m) `suchThat` ((>= 3) . nrows . unPropMatArrowheadSPD)

genSpM_ArrowheadSPD :: Int -> Gen (SpMatrix Double)
genSpM_ArrowheadSPD n = do
  diagElems <- vectorOf n $ choose (1.0, 10.0)
  lastRowElems <- vectorOf (n-1) $ choose (-5.0, 5.0)
  let diagSum = sum (map abs lastRowElems) + 1.0
      adjustedDiag = take (n-1) diagElems ++ [max (last diagElems) diagSum]
      diagMat = mkDiagonal n adjustedDiag
      lastRowMat = fromListSM (n, n) $ zip3 (replicate (n-1) (n-1)) [0..n-2] lastRowElems
      lastColMat = fromListSM (n, n) $ zip3 [0..n-2] (replicate (n-1) (n-1)) lastRowElems
  return $ diagMat ^+^ lastRowMat ^+^ lastColMat

-- | Arrowhead HPD matrix (Complex)
newtype PropMatArrowheadHPD a = PropMatArrowheadHPD {unPropMatArrowheadHPD :: SpMatrix a} deriving (Eq, Show)

instance Arbitrary (PropMatArrowheadHPD (Complex Double)) where
  arbitrary = sized (\m -> PropMatArrowheadHPD <$> genSpM_ArrowheadHPD_complex m) `suchThat` ((>= 3) . nrows . unPropMatArrowheadHPD)

genSpM_ArrowheadHPD_complex :: Int -> Gen (SpMatrix (Complex Double))
genSpM_ArrowheadHPD_complex n = do
  diagElems <- vectorOf n $ choose (1.0, 10.0)
  lastRowElemsReal <- vectorOf (n-1) $ choose (-5.0, 5.0)
  lastRowElemsImag <- vectorOf (n-1) $ choose (-5.0, 5.0)
  let lastRowElems = zipWith (:+) lastRowElemsReal lastRowElemsImag
      diagSum = sum (map magnitude lastRowElems) + 1.0
      adjustedDiag = take (n-1) diagElems ++ [max (last diagElems) diagSum]
      diagComplex = map (:+ 0.0) adjustedDiag
      diagMat = mkDiagonal n diagComplex
      lastRowMat = fromListSM (n, n) $ zip3 (replicate (n-1) (n-1)) [0..n-2] lastRowElems
      lastColMat = fromListSM (n, n) $ zip3 [0..n-2] (replicate (n-1) (n-1)) (map conjugate lastRowElems)
  return $ diagMat ^+^ lastRowMat ^+^ lastColMat

genSpMDense :: (Num a, Arbitrary a) => Int -> Int -> Gen (SpMatrix a)
genSpMDense m n = do
  xss <- vectorOf (m * n) arbitrary
  return $ fromListDenseSM m xss

genSpM :: (Num a, Arbitrary a) => Int -> Int -> Gen (SpMatrix a)
genSpM m n = do
  nnz <- choose (0, m * n `div` 2)
  entries <- vectorOf nnz genEntry
  return $ fromListSM (m, n) entries
  where
    genEntry = do
      i <- choose (0, m - 1)
      j <- choose (0, n - 1)
      x <- arbitrary
      return (i, j, x)


-- * Test Matrices

-- QR test matrices
tm2 :: SpMatrix Double
tm2 = fromListDenseSM 3 [12, 6, -4, -51, 167, 24, 4, -68, -41]

tm4 :: SpMatrix Double
tm4 = fromListSM (4, 4) [(0, 0, 10.0), (0, 1, 2.0), (1, 1, 5.0), (2, 2, 8.0), (3, 3, 3.0), (3, 1, 1.0)]

tm6 :: SpMatrix Double
tm6 = fromListDenseSM 4 [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

issueMatrix :: SpMatrix Double
issueMatrix = fromListDenseSM 4 [4,0,0,1, 0,4,1,0, 0,1,4,0, 1,0,0,4]

aa3cx :: SpMatrix (Complex Double)
aa3cx = fromListDenseSM 2 [4 :+ 0, 1 :+ (-1), 1 :+ 1, 4 :+ 0]

tmc4 :: SpMatrix (Complex Double)
tmc4 = fromListDenseSM 3 [1 :+ 0, 2 :+ 1, 3 :+ 2, 0 :+ 1, 2 :+ 0, 1 :+ 3, 3 :+ 2, 1 :+ 3, 1 :+ 0]

-- Cholesky test matrices
testSPD3x3 :: SpMatrix Double
testSPD3x3 = fromListDenseSM 3 [4, 12, -16, 12, 37, -43, -16, -43, 98]

tm7spd :: SpMatrix Double
tm7spd = a ^+^ b ^+^ c where
  n = 5
  a = mkSubDiagonal n 1 $ replicate (n-1) (-1)
  b = mkSubDiagonal n 0 $ replicate n 2
  c = mkSubDiagonal n (-1) $ replicate (n-1) (-1)

testHPD2x2 :: SpMatrix (Complex Double)
testHPD2x2 = fromListDenseSM 2 [4 :+ 0, 1 :+ (-1), 1 :+ 1, 2 :+ 0]

testHPD3x3 :: SpMatrix (Complex Double)
testHPD3x3 = fromListDenseSM 3 
  [4 :+ 0,     1 :+ (-1),  2 :+ 0,
   1 :+ 1,     5 :+ 0,     0 :+ (-2),
   2 :+ 0,     0 :+ 2,     6 :+ 0]

railsMatrix :: SpMatrix Double
railsMatrix = fromListSM (8, 8)
  [ (0, 0, 2.0), (0, 7, 1.0)
  , (1, 1, 2.0), (1, 7, 1.0)
  , (2, 2, 2.0), (2, 7, 1.0)
  , (3, 3, 2.0), (3, 7, 1.0)
  , (4, 4, 2.0), (4, 7, 1.0)
  , (5, 5, 2.0), (5, 7, 1.0)
  , (6, 6, 2.0), (6, 7, 1.0)
  , (7, 0, 1.0), (7, 1, 1.0), (7, 2, 1.0), (7, 3, 1.0)
  , (7, 4, 1.0), (7, 5, 1.0), (7, 6, 1.0), (7, 7, 8.0)
  ]
