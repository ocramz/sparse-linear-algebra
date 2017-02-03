{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts, TypeFamilies, MultiParamTypeClasses, FlexibleInstances  #-}
-- {-# OPTIONS_GHC -O2 -rtsopts -with-rtsopts=-K32m -prof#-}
{-|

This module exposes the high-level functionality of the library.

-}

module Numeric.LinearAlgebra.Sparse
       (
         -- * Matrix factorizations
         qr,
         lu,
         chol,
         -- * Condition number
         conditionNumberSM,
         -- * Householder reflection
         hhMat, hhRefl,
         -- * Householder bidiagonalization

         -- * Givens' rotation
         givens,
         -- * Arnoldi iteration
         arnoldi, 
         -- * Eigensolvers
         -- eigsQR,
         -- eigRayleigh,
         -- * Linear solvers
         -- ** Iterative methods
         -- linSolve0, LinSolveMethod(..), (<\>),
         -- pinv,
         -- ** Direct methods
         luSolve, triLowerSolve, triUpperSolve,
         -- * Preconditioners
         ilu0, mSsor,
         -- * Matrix partitioning
         diagPartitions,
         -- -- * Random arrays
         -- randArray,
         -- -- * Random matrices and vectors
         -- randMat, randVec, 
         -- -- ** Sparse "
         -- randSpMat, randSpVec,
         -- * Iteration combinators
         modifyInspectN, modifyInspectGuarded,
         runAppendN', untilConverged,
         diffSqL
       )
       where


import Control.Exception.Common
import Data.Sparse.Common

import Control.Exception
import Control.Monad.Catch
import Data.Typeable

-- import Control.Applicative ((<|>))

-- import Control.Monad (replicateM)
import Control.Monad.State.Strict
-- import Control.Monad.Writer
-- import Control.Monad.Trans
import qualified Control.Monad.Trans.State  as MTS -- (runStateT)
-- import Control.Monad.Trans.Writer (runWriterT)
import Data.Complex

import Data.VectorSpace hiding (magnitude)

import qualified Data.Sparse.Internal.IntM as I
-- import Data.Utils.StrictFold (foldlStrict) -- hidden in `containers`

-- import qualified System.Random.MWC as MWC
-- import qualified System.Random.MWC.Distributions as MWC

-- import Data.Monoid
-- import qualified Data.Foldable as F
-- import qualified Data.Traversable as T

-- import qualified Data.List as L
import Data.Maybe

import qualified Data.Vector as V



-- | A lumped constraint for the numerical types  
type Data x = (Epsilon x, Elt x, Show x, Ord x)





-- * Matrix condition number

-- |uses the R matrix from the QR factorization
conditionNumberSM :: (MonadThrow m, MatrixRing (SpMatrix a), Data a, Typeable a) =>
     SpMatrix a -> m a
conditionNumberSM m | nearZero lmin = throwM (HugeConditionNumber kappa)
                    | otherwise = return kappa where
  kappa = lmax / lmin
  (_, r) = qr m
  u = extractDiagDense r  -- FIXME : need to extract with default element 0 
  lmax = abs (maximum u)
  lmin = abs (minimum u)







-- * Householder transformation

hhMat :: Num a => a -> SpVector a -> SpMatrix a
hhMat beta x = eye n ^-^ beta `scale` (x >< x) where
  n = dim x


{-| a vector `x` uniquely defines an orthogonal plane; the Householder operator reflects any point `v` with respect to this plane:
 v' = (I - 2 x >< x) v
-}
hhRefl :: Num a => SpVector a -> SpMatrix a
hhRefl = hhMat (fromInteger 2)











-- * Givens rotation matrix




-- class Signum x where
--   type SigFloat x :: *
--   signum' :: x -> SigFloat x

-- instance Signum Double where {type SigFloat Double = Double; signum' = signum}
-- instance Signum (Complex Double) where {type SigFloat (Complex Double) = Double; signum' = signum . magnitude}


-- -- -- | Givens coefficients (using stable algorithm shown in  Anderson, Edward (4 December 2000). "Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem". LAPACK Working Note)
-- -- -- givensCoef0 :: (Ord b, Floating a, Eq a) => (a -> b) -> a -> a -> (a, a, a)
-- -- givensCoef0
-- --   :: (Signum t, Ord a, Floating t, Eq t) =>
-- --      (t -> a) -> t -> t -> (t, t, t)
-- givensCoef0 ff a b  -- returns (c, s, r) where r = norm (a, b)
--   | b==0 = (signum' a, 0, abs a)
--   | a==0 = (0, signum' b, abs b)
--   | ff a > ff b = let t = b/a
--                       u = signum' a * abs ( sqrt (1 + t**2))
--                   in (1/u, - t/u, a*u)
--   | otherwise = let t = a/b
--                     u = signum' b * abs ( sqrt (1 + t**2))
--                 in (t/u, - 1/u, b*u)

-- -- | Givens coefficients, real-valued
-- givensCoef :: (Floating a, Ord a) => a -> a -> (a, a, a)
-- givensCoef = givensCoef0 abs

-- -- | Givens coefficients, complex-valued
-- givensCoefC :: RealFloat a => Complex a -> Complex a -> (Complex a, Complex a, Complex a)
-- givensCoefC = givensCoef0 magnitude

-- | Givens coefficients
givensCoef :: (Elt t, Floating t) => t -> t -> (t, t, t)
givensCoef a b = (c0/d, s0/d, d) where
  c0 = conj a
  s0 = - conj b
  d = hypot c0 s0 -- sqrt $ (norm1 c0)**2 + (norm1 s0)**2

-- | Stable hypotenuse calculation
hypot :: Floating a => a -> a -> a
hypot x y = abs x * (sqrt (1 + (y/x)**2))

                  
{- |
Givens method, row version: choose other row index i' s.t. i' is :
* below the diagonal
* corresponding element is nonzero

QR.C1 ) To zero out entry A(i, j) we must find row k such that A(k, j) is
non-zero but A has zeros in row k for all columns less than j.

NB: the current version is quite inefficient in that:
1. the Givens' matrix `G_i` is different from Identity only in 4 entries
2. at each iteration `i` we multiply `G_i` by the previous partial result `M`. Since this corresponds to a rotation, and the `givensCoef` function already computes the value of the resulting non-zero component (output `r`), `G_i ## M` can be simplified by just changing two entries of `M` (i.e. zeroing one out and changing the other into `r`).
-}
{-# inline givens #-}
givens :: (Elt a, Floating a) => SpMatrix a -> IxRow -> IxCol -> SpMatrix a
givens mm i j 
  | isValidIxSM mm (i,j) && nrows mm >= ncols mm =
       -- fromListSM' [(i,i,c),(j,j,conj c),(j,i,- (conj s)),(i,j,s)] (eye (nrows mm))
       fromListSM'
         [(i,i, c), (j,j, conj c), (j,i, - conj s), (i,j, s)] (eye (nrows mm))       
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    i' = head $ fromMaybe (error $ "givens: no compatible rows for entry " ++ show (i,j)) (candidateRows (immSM mm) i j)
    a = mm @@ (i', j)
    b = mm @@ (i, j)   -- element to zero out

-- |Returns a set of rows {k} that satisfy QR.C1
-- candidateRows :: IM.IntMap (IM.IntMap a) -> IxRow -> IxCol -> Maybe [IM.Key]
candidateRows mm i j | null u = Nothing
                     | otherwise = Just (I.keys u) where
  u = I.filterWithKey (\irow row -> irow /= i &&
                                    firstNonZeroColumn row j) mm

-- |Is the `k`th the first nonzero column in the row?
{-# inline firstNonZeroColumn #-}
-- firstNonZeroColumn :: IM.IntMap a -> IxRow -> Bool
firstNonZeroColumn mm k = isJust (I.lookup k mm) &&
                          isNothing (I.lookupLT k mm)






{-# inline givens' #-}
givens' :: (Elt a, MonadThrow m) => SpMatrix a -> Int -> Int -> m (SpMatrix a)
givens' mm i j 
  | isValidIxSM mm (i,j) && nrows mm >= ncols mm = do
      i' <- candidateRows' (immSM mm) i j
      return $ givensMat mm i i' j
  | otherwise = throwM (OOBIxsError "Givens" [i, j])
  where
  givensMat mm i i' j = fromListSM'
           [(i,i, c), (j,j, conj c), (j,i, - conj s), (i,j, s)] $ eye (nrows mm)
           where
             (c, s, _) = givensCoef a b
             a = mm @@ (i', j)
             b = mm @@ (i, j)   -- element to zero out
  candidateRows' mm i j | null u = throwM (OOBNoCompatRows "Givens" (i,j))
                        | otherwise = return $ head (I.keys u) where
    u = I.filterWithKey (\irow row -> irow /= i &&
                                      firstNZColumn row j) mm
    firstNZColumn m k = isJust (I.lookup k m) &&
                        isNothing (I.lookupLT k m)
    




                          




-- * QR decomposition


-- | Given a matrix A, returns a pair of matrices (Q, R) such that Q R = A, where Q is orthogonal and R is upper triangular. Applies Givens rotation iteratively to zero out sub-diagonal elements.
-- qr :: (Elt a, MatrixRing (SpMatrix a), Epsilon a, Floating a) =>
--      SpMatrix a -> (SpMatrix a, SpMatrix a)
qr mm = (transpose qt, r) where
  (qt, r, _) = execState (modifyUntil haltf stepf) gminit
  haltf (_, _, iis) = null iis
  stepf (qmatt, m, iis) = (qmatt', m', tail iis) where
    (i, j) = head iis
    g = givens m i j
    qmatt' = g #~# qmatt  -- update Q'
    m' = g #~# m          -- update R
  gminit = (eye (nrows mm), mm, subdiagIndicesSM mm)


qr' :: (Elt a, MatrixRing (SpMatrix a), Epsilon a, MonadThrow m) =>
     SpMatrix a -> m (SpMatrix a, SpMatrix a)
qr' mm = do 
    (qt, r, _) <- MTS.execStateT (modifyUntilM haltf qrstepf) gminit
    return (transpose qt, r) 
  where
    gminit = (eye (nrows mm), mm, subdiagIndicesSM mm)
    haltf (_, _, iis) = null iis
    qrstepf (qmatt, m, iis) = do
        let (i, j) = head iis
        g <- givens' m i j
        let
          qmatt' = g #~# qmatt  -- update Q'
          m' = g #~# m          -- update R
        return (qmatt', m', tail iis)

  








-- * Eigenvalue algorithms

-- ** QR algorithm

-- -- | `eigsQR n mm` performs `n` iterations of the QR algorithm on matrix `mm`, and returns a SpVector containing all eigenvalues
-- eigsQR :: (Elt a, Normed (SpVector a), MatrixRing (SpMatrix a), Epsilon (Scalar (SpVector a)), Epsilon a, Floating (Scalar (SpVector a)), Floating a) =>
--      Int -> SpMatrix a -> SpVector a
-- eigsQR nitermax m = extractDiagDense $ execState (convergtest eigsStep) m where
--   eigsStep mm = q #~^# (m ## q) -- r #~# q
--    where
--     (q, _) = qr mm


eigsQRconvergtest nitermax g = modifyInspectN nitermax f g where
      f [m1, m2] = let dm1 = extractDiagDense m1
                       dm2 = extractDiagDense m2
                   in nearZero $ norm2' (dm1 ^-^ dm2)
      f _ = False


eigsQR' nitermax m = untilConvergedGM "eigsQR" nitermax pf (const True) -- stepf
  where
    pf = extractDiagDense
    -- stepf mm = do
    --   (q, _) <- qr mm
    --   return $ q #~^# (m ## q) -- r #~# q
    -- summf [m1, m2] = let dm1 = extractDiagDense m1
    --                      dm2 = extractDiagDense m2
    --                  in (dm1 ^-^ dm2)
    -- -- summf _ = False    
      




-- eigsQR' nitermax m = do
--   vm <- MTS.execStateT (convergtest stepf) m
--   return $ extractDiagDense vm
--   where
--     stepf mm = do
--       (q, _) <- qr mm
--       return $ q #~^# (m ## q) -- r #~# q
--     convergtest g = modifyInspectN nitermax f g where
--       f [m1, m2] = let dm1 = extractDiagDense m1
--                        dm2 = extractDiagDense m2
--                    in nearZero $ norm2' (dm1 ^-^ dm2)
--       f _ = False



-- ** Rayleigh iteration

-- | `eigsRayleigh n mm` performs `n` iterations of the Rayleigh algorithm on matrix `mm` and returns the eigenpair closest to the initialization. It displays cubic-order convergence, but it also requires an educated guess on the initial eigenpair
-- eigRayleigh :: (Epsilon a, Floating a) => Int -- max # iterations
--      -> SpMatrix a           -- matrix
--      -> (SpVector a, a) -- initial guess of (eigenvector, eigenvalue)
--      -> (SpVector a, a) -- final estimate of (eigenvector, eigenvalue)
-- eigRayleigh :: Int -> SpMatrix Double -> (SpVector Double, Double) -> (SpVector Double, Double)
-- eigRayleigh :: Int -> SpMatrix (Complex Double) -> (SpVector (Complex Double), (Complex Double)) -> (SpVector (Complex Double), (Complex Double))
-- eigRayleigh nitermax m = execState (convergtest (rayleighStep m)) where
--   convergtest g = modifyInspectN nitermax f g where
--     f [(b1, _), (b2, _)] = nearZero $ norm2' (b2 ^-^ b1)
--   rayleighStep aa (b, mu) = (b', mu') where
--       ii = eye (nrows aa)
--       nom = (aa ^-^ (mu `matScale` ii)) <\> b
--       b' = normalize2' nom
--       mu' = (b' <.> (aa #> b')) / (b' <.> b')




-- * Householder vector 

-- (Golub & Van Loan, Alg. 5.1.1, function `house`)
hhV x = (v, beta) where
  tx = tailSV x
  sigma = tx <.> tx 
  vtemp = singletonSV 1 `concatSV` tx
  (v, beta) | nearZero sigma = (vtemp, 0)
            | otherwise = let mu = sqrt (headSV x**2 + sigma)
                              xh = headSV x
                              vh | mag xh <= 1 = xh - mu
                                 | otherwise = - sigma / (xh + mu)
                              vnew = (1 / vh) `scale` insertSpVector 0 vh vtemp     
                          in (vnew, 2 * xh**2 / (sigma + vh**2))

                         



-- * Bidiagonalization

-- | Golub-Kahan-Lanczos bidiagonalization (see "Restarted Lanczos Bidiagonalization for the SVD", SLEPc STR-8, http://slepc.upv.es/documentation/manual.htm )
gklBidiag aa q1nn | dim q1nn == n = (pp, bb, qq)
                  | otherwise = error "hhBidiag : dimension mismatch. Provide q1 compatible with aa #> q1"
  where
  (m,n) = dim aa
  aat = transpose aa
  bb = fromListSM (n,n) bl
  (ql, _, pl, _, pp, bl, qq) = execState (modifyUntil tf bidiagStep) bidiagInit
  tf (_, _, _, i, _, _, _) = i == n 
  bidiagInit = (q2n, beta1, p1n, 1 :: Int, pp, bb', qq)
   where
    q1 = normalize2' q1nn
    
    p1 = aa #> q1
    alpha1 = norm2' p1
    p1n = p1 ./ alpha1
    
    q2 = (aat #> p1) ^-^ (alpha1 .* q1)
    beta1 = norm2' q2
    q2n = q2 ./ beta1
    
    pp = insertCol (zeroSM m n) p1n 0
    qq = insertCol (zeroSM n n) q2n 0
    bb' = [(0, 0, alpha1)]
  bidiagStep (qj , betajm, pjm , j   , pp,  bb, qq ) =
             (qjp, betaj , pj, succ j, pp', bb', qq') where

    u = (aa #> qj) ^-^ (betajm .* pjm)
    alphaj = norm2' u
    pj = u ./ alphaj
    
    v = (aat #> pj) ^-^ (alphaj .* qj)
    betaj = norm2' v
    qjp = v ./ betaj
  
    pp' = insertCol pp pj j
    bb' = [(j-1, j, betaj),
           (j ,j, alphaj)] ++ bb
    qq' = insertCol qq qjp j


-- fromColsL :: [SpVector a] -> SpMatrix a
-- fromColsL = fromCols . V.fromList

toCols :: SpMatrix a -> [SpVector a]
toCols aa = map (extractCol aa) [0 .. n-1] where
  (m,n) = dim aa



-- | Example 5.4.2 from G & VL
aa1 :: SpMatrix Double
aa1 = transpose $ fromListDenseSM 3 [1..12]



-- aa1 :: SpMatrix Double
-- aa1 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]




-- * SVD

{- Golub & Van Loan, sec 8.6.2 (p 452 segg.)

SVD of A, Golub-Kahan method

* reduce A to upper bidiagonal form B (Alg. 5.4.2, Householder bidiagonalization)
* compute SVD of B (implicit-shift QR step applied to B^T B, Alg. 8.3.2)

-}








-- * Cholesky factorization

-- | Given a positive semidefinite matrix A, returns a lower-triangular matrix L such that L L^T = A . This is an implementation of the Cholesky–Banachiewicz algorithm, i.e. proceeding row by row from the upper-left corner.
chol :: (Epsilon a, Floating a, Elt a) => SpMatrix a -> SpMatrix a
chol aa = lfin where
  (_, lfin) = execState (modifyUntil q cholUpd) cholInit
  q (i, _) = i == nrows aa              -- stopping criterion
  cholInit = cholUpd (0, zeroSM n n)    -- initialization
  n = nrows aa
  cholUpd (i, ll) = (i + 1, ll') where
    ll' = cholDiagUpd (cholSDRowUpd ll) -- first upd subdiagonal entries in the row
    cholSDRowUpd ll_ = insertRow ll_ lrs i where
       lrs = fromListSV (i + 1) $ onRangeSparse (cholSubDiag ll i) [0 .. i-1]
    cholDiagUpd ll_ = insertSpMatrix i i (cholDiag ll_ i) ll_ 
  cholSubDiag ll i j = 1/ljj*(aij - inn) where
    ljj = ll@@(j, j)
    aij = aa@@(i, j)
    inn = contractSub ll ll i j (j - 1)
  cholDiag ll i | i == 0 = sqrt aai
                | otherwise = sqrt $ aai - sum (fmap (**2) lrow) -- i > 0
    where
        aai = aa@@(i,i)
        lrow = ifilterSV (\j _ -> j < i) (extractRow ll i) -- sub-diagonal elems of L

















-- * LU factorization
-- ** Doolittle algorithm
{- Doolittle algorithm for factoring A' = P A, where P is a permutation matrix such that A' has a nonzero as its (0, 0) entry -}

-- | Given a matrix A, returns a pair of matrices (L, U) where L is lower triangular and U is upper triangular such that L U = A
lu :: (VectorSpace (SpVector t), Epsilon t, Fractional t, Elt t) =>
     SpMatrix t -> (SpMatrix t, SpMatrix t)
lu aa = (lf, ufin) where
  (ixf, lf, uf) = execState (modifyUntil q luUpd) luInit
  ufin = uUpdSparse (ixf, lf, uf) -- final U update
  q (i, _, _) = i == (nrows aa - 1)
  n = nrows aa
  luInit = (1, l0, u0) where
--  l0 = insertCol (eye n) ((1/u00) .* extractSubCol aa 0 (1,n - 1)) 0  -- initial L
    l0 = insertCol (eye n) (scale (1/u00) (extractSubCol aa 0 (1, n-1))) 0
    u0 = insertRow (zeroSM n n) (extractRow aa 0) 0                     -- initial U
    u00 = u0 @@ (0,0)  -- make sure this is non-zero by applying permutation
  luUpd (i, l, u) = (i + 1, l', u') where
    u' = uUpdSparse (i, l, u)  -- update U
    l' = lUpdSparse (i, l, u') -- update L
  uUpdSparse (ix, lmat, umat) = insertRow umat (fromListSV n us) ix where
    us = onRangeSparse (solveForUij ix) [ix .. n - 1]
    solveForUij i j = a - p where
      a = aa @@! (i, j)
      p = contractSub lmat umat i j (i - 1)
  lUpdSparse (ix, lmat, umat) = insertCol lmat (fromListSV n ls) ix where
    ls = onRangeSparse (`solveForLij` ix) [ix + 1 .. n - 1]
    solveForLij i j
     | isNz ujj = (a - p)/ujj
     | otherwise =
        error $ unwords ["solveForLij : U",
                       show (j ,j),
                       "is close to 0. Permute rows in order to have a nonzero diagonal of U"]
      where
       a = aa @@! (i, j)
       ujj = umat @@! (j , j)   -- NB this must be /= 0
       p = contractSub lmat umat i j (i - 1)



-- | Apply a function over a range of integer indices, zip the result with it and filter out the almost-zero entries
onRangeSparse :: Epsilon b => (Int -> b) -> [Int] -> [(Int, b)]
onRangeSparse f ixs = filter (isNz . snd) $ zip ixs $ map f ixs




















-- -- Produces the permutation matrix necessary to have a nonzero in position (iref, jref). This is used in the LU factorization
-- permutAA :: Num b => IxRow -> IxCol -> SpMatrix a -> Maybe (SpMatrix b)
-- permutAA iref jref (SM (nro,_) mm) 
--   | isJust (lookupIM2 iref jref mm) = Nothing -- eye nro
--   | otherwise = Just $ permutationSM nro [head u] where
--       u = IM.keys (ifilterIM2 ff mm)
--       ff i j _ = i /= iref &&
--                  j == jref
















-- * Arnoldi iteration

-- | Given a matrix A, a vector b and a positive integer `n`, this procedure finds the basis of an order `n` Krylov subspace (as the columns of matrix Q), along with an upper Hessenberg matrix H, such that A = Q^T H Q.
-- At the i`th iteration, it finds (i + 1) coefficients (the i`th column of the Hessenberg matrix H) and the (i + 1)`th Krylov vector.

arnoldi :: (MatrixType (SpVector a) ~ SpMatrix a, V (SpVector a) , Scalar (SpVector a) ~ a, Floating a, Epsilon a) =>
     SpMatrix a
     -> SpVector a
     -> Int
     -> (SpMatrix a, SpMatrix a)
arnoldi aa b kn = (fromCols qvfin, fromListSM (nmax + 1, nmax) hhfin)
  where
  (qvfin, hhfin, nmax, _) = execState (modifyUntil tf arnoldiStep) arnInit 
  tf (_, _, ii, fbreak) = ii == kn || fbreak -- termination criterion
  (m, n) = (nrows aa, ncols aa)
  arnInit = (qv1, hh1, 1, False) where
      q0 = normalize2 b   -- starting basis vector
      aq0 = aa #> q0       -- A q0
      h11 = q0 `dot` aq0          
      q1nn = (aq0 ^-^ (h11 .* q0))
      hh1 = V.fromList [(0, 0, h11), (1, 0, h21)] where        
        h21 = norm2' q1nn
      q1 = normalize2 q1nn       -- q1 `dot` q0 ~ 0
      qv1 = V.fromList [q0, q1]
  arnoldiStep (qv, hh, i, _) = (qv', hh', i + 1, breakf) where
    qi = V.last qv
    aqi = aa #> qi
    hhcoli = fmap (`dot` aqi) qv -- H_{1, i}, H_{2, i}, .. , H_{m + 1, i}
    zv = zeroSV m
    qipnn =
      aqi ^-^ (V.foldl' (^+^) zv (V.zipWith (.*) hhcoli qv)) -- unnormalized q_{i+1}
    qipnorm = norm2' qipnn      -- normalization factor H_{i+1, i}
    qip = normalize2 qipnn              -- q_{i + 1}
    hh' = (V.++) hh (indexed2 $ V.snoc hhcoli qipnorm) where -- update H
      indexed2 v = V.zip3 ii jj v
      ii = V.fromList [0 .. n]    -- nth col of upper Hessenberg has `n+1` nz
      jj = V.replicate (n + 1) i  -- `n+1` replicas of `i`
    qv' = V.snoc qv qip        -- append q_{i+1} to Krylov basis Q_i
    breakf | nearZero qipnorm = True  -- breakdown condition
           | otherwise = False







-- * Preconditioning

-- | Partition a matrix into strictly subdiagonal, diagonal and strictly superdiagonal parts
diagPartitions :: SpMatrix a -> (SpMatrix a, SpMatrix a, SpMatrix a)
diagPartitions aa = (e,d,f) where
  e = extractSubDiag aa
  d = extractDiag aa
  f = extractSuperDiag aa


-- -- ** Jacobi preconditioner

-- -- | Returns the reciprocal of the diagonal 
-- jacobiPreconditioner :: SpMatrix Double -> SpMatrix Double
-- jacobiPreconditioner = reciprocal . extractDiag


-- ** Incomplete LU

-- | Used for Incomplete LU : remove entries in `m` corresponding to zero entries in `m2`
ilu0 :: (Elt a, VectorSpace (SpVector a), Epsilon a) =>
   SpMatrix a -> (SpMatrix a, SpMatrix a)
ilu0 aa = (lh, uh) where
  (l, u) = lu aa
  lh = sparsifyLU l aa
  uh = sparsifyLU u aa
  sparsifyLU m m2 = ifilterSM f m where
    f i j _ = isJust (lookupSM m2 i j)


-- ** SSOR

-- | `mSsor aa omega` : if `omega = 1` it returns the symmetric Gauss-Seidel preconditioner. When ω = 1, the SOR reduces to Gauss-Seidel; when ω > 1 and ω < 1, it corresponds to over-relaxation and under-relaxation, respectively.
mSsor :: (MatrixRing (SpMatrix b), Fractional b) =>
     SpMatrix b -> b -> (SpMatrix b, SpMatrix b)
mSsor aa omega = (l, r) where
  (e, d, f) = diagPartitions aa
  n = nrows e
  l = (eye n ^-^ scale omega e) ## reciprocal d
  r = d ^-^ scale omega f 








-- * Linear solver, LU-based

-- | Direct solver based on a triangular factorization of the system matrix.
luSolve :: (Scalar (SpVector t) ~ t, MonadThrow m, Elt t, InnerSpace (SpVector t), Epsilon t) =>
     SpMatrix t -> SpMatrix t -> SpVector t -> m (SpVector t)
luSolve ll uu b
  | isLowerTriSM ll && isUpperTriSM uu = return $ triUpperSolve uu (triLowerSolve ll b)
  | otherwise = throwM (NonTriangularException "luSolve")-- error "luSolve : factors must be triangular matrices" 

-- triLowerSolve :: (Epsilon a, RealFrac a, Elt a, AdditiveGroup a) => SpMatrix a -> SpVector a -> SpVector a
triLowerSolve :: (Scalar (SpVector t) ~ t, InnerSpace (SpVector t), Epsilon t, Elt t) => SpMatrix t -> SpVector t -> SpVector t
triLowerSolve ll b = sparsifySV v where
  (v, _) = execState (modifyUntil q lStep) lInit where
  q (_, i) = i == dim b
  lStep (ww, i) = (ww', i + 1) where
    lii = ll @@ (i, i)
    bi = b @@ i
    wi = (bi - r)/lii where
      r = extractSubRow ll i (0, i-1) `dot` takeSV i ww
    ww' = insertSpVector i wi ww
  lInit = (ww0, 1) where
    l00 = ll @@ (0, 0)
    b0 = b @@ 0
    w0 = b0 / l00
    ww0 = insertSpVector 0 w0 $ zeroSV (dim b)  

-- | NB in the computation of `xi` we must rebalance the subrow indices because `dropSV` does that too, in order to take the inner product with consistent index pairs

triUpperSolve ::
  (Scalar (SpVector t) ~ t, InnerSpace (SpVector t), Epsilon t, Elt t) =>
  SpMatrix t -> SpVector t -> SpVector t
triUpperSolve uu w = sparsifySV x where
  (x, _) = execState (modifyUntil q uStep) uInit
  q (_, i) = i == (- 1)
  uStep (xx, i) = (xx', i - 1) where
    uii = uu @@ (i, i)
    wi = w @@ i
    xi = (wi - r) / uii where
        r = extractSubRow_RK uu i (i + 1, dim w - 1) `dot` dropSV (i + 1) xx
    xx' = insertSpVector i xi xx
  uInit = (xx0, i - 1) where
    i = dim w - 1
    u00 = uu @@ (i, i)
    w0 = w @@ i
    x0 = w0 / u00
    xx0 = insertSpVector i x0 $ zeroSV (dim w)






-- * Iterative linear solvers


-- ** GMRES

-- | Given a linear system `A x = b` where `A` is an (m x m) real-valued matrix, the GMRES method finds an approximate solution `xhat` such that the Euclidean norm of the residual `A xhat - b` is minimized. `xhat` is spanned by the order-`n` Krylov subspace of (A, b).
-- In this implementation:
-- 1) the Arnoldi iteration is carried out until numerical breakdown (therefore yielding _at_most_ `m+1` Krylov basis vectors)
-- 2) the resulting Hessenberg matrix H is factorized in QR form (H = Q R)
-- 3) the Krylov-subspace solution `yhat` is found by backsubstitution (since R is upper-triangular)
-- 4) the approximate solution in the original space `xhat` is computed using the Krylov basis, `xhat = Q_n yhat`
--
-- Many optimizations are possible, for example interleaving the QR factorization (and the subsequent triangular solve) with the Arnoldi process (and employing an updating QR factorization which only requires one Givens' rotation at every update). 

gmres :: (Scalar (SpVector t) ~ t,
          MatrixType (SpVector t) ~ SpMatrix t,
      Elt t, V (SpVector t), MatrixRing (SpMatrix t), Epsilon t) =>
     SpMatrix t -> SpVector t -> SpVector t
gmres aa b = qa' #> yhat where
  m = ncols aa
  (qa, ha) = arnoldi aa b m   -- at most m steps of Arnoldi (aa, b)
  -- b' = (transposeSe qa) #> b
  b' = norm2' b .* (ei mp1 1)  -- b rotated back to canonical basis by Q^T
     where mp1 = nrows ha     -- = 1 + (# Arnoldi iterations)
  (qh, rh) = qr ha            -- QR factors of H
  yhat = triUpperSolve rh' rhs' where
    rhs' = takeSV (dim b' - 1) (transpose qh #> b')
    rh' = takeRows (nrows rh - 1) rh -- last row of `rh` is 0
  qa' = takeCols (ncols qa - 1) qa   -- we don't use last column of Krylov basis















-- ** CGNE

data CGNE a =
  CGNE {_xCgne , _rCgne, _pCgne :: SpVector a} deriving Eq
instance Show a => Show (CGNE a) where
    show (CGNE x r p) = "x = " ++ show x ++ "\n" ++
                        "r = " ++ show r ++ "\n" ++
                        "p = " ++ show p ++ "\n"

-- cgne :: (RealFrac a, Elt a, AdditiveGroup a, Epsilon a) => SpMatrix a -> SpVector a -> SpVector a -> CGNE a
cgne aa b x0 = execState (untilConverged _xCgne (cgneStep aa)) cgneInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = transposeSM aa #> r0
  cgneInit = CGNE x0 r0 p0
  cgneStep aa (CGNE x r p) = CGNE x1 r1 p1 where
    alphai = (r `dot` r) / (p `dot` p)
    x1 = x ^+^ (alphai .* p)
    r1 = r ^-^ (alphai .* (aa #> p))
    beta = (r1 `dot` r1) / (r `dot` r)
    p1 = transposeSM aa #> r ^+^ (beta .* p)


-- ** TFQMR

-- tfqmr :: (Epsilon a, Floating a) => SpMatrix a -> SpVector a -> SpVector a -> TFQMR a
tfqmr aa b x0 = execState (untilConverged _xTfq (tfqmrStep aa r0)) tfqmrInit where
  n = dim b
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  w0 = r0
  u0 = r0
  v0 = aa #> u0
  d0 = zeroSV n
  r0hat = r0
  rho0 = r0hat `dot` r0
  alpha0 = rho0 / (v0 `dot` r0hat)
  m = 0
  tau0 = norm2' r0
  theta0 = 0
  eta0 = 0
  tfqmrInit = TFQMR x0 w0 u0 v0 d0 m tau0 theta0 eta0 rho0 alpha0
  tfqmrStep aa r0hat (TFQMR x w u v d m tau theta eta rho alpha) =
    TFQMR x1 w1 u1 v1 d1 (m+1) tau1 theta1 eta1 rho1 alpha1
    where
    w1 = w ^-^ (alpha .* (aa #> u))
    d1 = u ^+^ ((theta**2/alpha*eta) .* d)
    theta1 = norm2' w1 / tau
    c = recip $ sqrt (1 + theta1**2)
    tau1 = tau * theta1 * c
    eta1 = c**2 * alpha
    x1 = x^+^ (eta1 .* d1)
    (alpha1, u1, rho1, v1)
      | even m = let
                     alpha' = rho / (v `dot` r0hat)
                     u' = u ^-^ (alpha' .* v)
                 in
                     (alpha', u', rho, v)
      | otherwise = let
                     rho' = w1 `dot` r0hat
                     beta = rho'/rho
                     u' = w1 ^+^ (beta .* u)
                     v' = (aa #> u') ^+^ (beta .* (aa #> u ^+^ (beta .* v)) )
                    in (alpha, u', rho', v')

data TFQMR a =
  TFQMR { _xTfq, _wTfq, _uTfq, _vTfq, _dTfq :: SpVector a,
          _mTfq :: Int,
          _tauTfq, _thetaTfq, _etaTfq, _rhoTfq, _alphaTfq :: a}
  deriving Eq
instance Show a => Show (TFQMR a) where
    show (TFQMR x _ _ _ _ _ _ _ _ _ _) = "x = " ++ show x ++ "\n"




-- ** BCG

data BCG a =
  BCG { _xBcg, _rBcg, _rHatBcg, _pBcg, _pHatBcg :: SpVector a } deriving Eq

-- bcg :: (Epsilon a, Fractional a) => SpMatrix a -> SpVector a -> SpVector a -> BCG a
bcg :: (V (SpVector a), MatrixType (SpVector a) ~ SpMatrix a, 
      Fractional (Scalar (SpVector a))) =>
     SpMatrix a -> SpVector a -> SpVector a -> BCG a
bcg aa b x0 = execState (untilConverged _xBcg (bcgStep aa)) bcgInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  r0hat = r0
  p0 = r0
  p0hat = r0
  bcgInit = BCG x0 r0 r0hat p0 p0hat

bcgStep aa (BCG x r rhat p phat) = BCG x1 r1 rhat1 p1 phat1 where
    aap = aa #> p
    alpha = (r `dot` rhat) / (aap `dot` phat)
    x1 = x ^+^ (alpha .* p)
    r1 = r ^-^ (alpha .* aap)
    rhat1 = rhat ^-^ (alpha .* (transposeSM aa #> phat))
    beta = (r1 `dot` rhat1) / (r `dot` rhat)
    p1 = r1 ^+^ (beta .* p)
    phat1 = rhat1 ^+^ (beta .* phat)

instance Show a => Show (BCG a) where
  show (BCG x r rhat p phat) = "x = " ++ show x ++ "\n" ++
                       "r = " ++ show r ++ "\n" ++
                       "r_hat = " ++ show rhat ++ "\n" ++
                       "p = " ++ show p ++ "\n" ++
                       "p_hat = " ++ show phat ++ "\n"


-- ** CGS

data CGS a = CGS { _x, _r, _p, _u :: SpVector a} deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
-- cgs :: (Epsilon a, Fractional a) =>
--   SpMatrix a -> SpVector a -> SpVector a -> SpVector a -> CGS a
cgs aa b x0 rhat =
  execState (untilConverged _x (cgsStep aa rhat)) (cgsInit aa b x0)
  
cgsInit aa b x0 = CGS x0 r0 r0 r0 where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution

cgsStep :: (V (SpVector a), Fractional (Scalar (SpVector a))) =>
     MatrixType (SpVector a) -> SpVector a -> CGS a -> CGS a
cgsStep aa rhat (CGS x r p u) = CGS xj1 rj1 pj1 uj1
    where
    aap = aa #> p
    alphaj = (r `dot` rhat) / (aap `dot` rhat)
    q = u ^-^ (alphaj .* aap)
    xj1 = x ^+^ (alphaj .* (u ^+^ q))         -- updated solution
    rj1 = r ^-^ (alphaj .* (aa #> (u ^+^ q))) -- updated residual
    betaj = (rj1 `dot` rhat) / (r `dot` rhat)
    uj1 = rj1 ^+^ (betaj .* q)
    pj1 = uj1 ^+^ (betaj .* (q ^+^ (betaj .* p)))


instance (Show a) => Show (CGS a) where
  show (CGS x r p u) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n" ++
                                "u = " ++ show u ++ "\n"





-- ** BiCGSTAB

-- _aa :: SpMatrix Double,    -- matrix
-- _b :: SpVector Double,     -- rhs
-- _r0 :: SpVector Double,    -- initial residual
-- _r0hat :: SpVector Double, -- candidate solution: r0hat `dot` r0 >= 0

-- | one step of BiCGSTAB
-- bicgstabStep :: Fractional a => SpMatrix a -> SpVector a -> BICGSTAB a -> BICGSTAB a


data BICGSTAB a =
  BICGSTAB { _xBicgstab, _rBicgstab, _pBicgstab :: SpVector a} deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
-- bicgstab ::
--   (Epsilon a, Fractional a) =>
--      SpMatrix a -> SpVector a -> SpVector a -> SpVector a -> BICGSTAB a
bicgstab aa b x0 r0hat =
  execState (untilConverged _xBicgstab (bicgstabStep aa r0hat)) (bicgsInit aa b x0)

bicgstabStep :: (V (SpVector a), Fractional (Scalar (SpVector a))) =>
     MatrixType (SpVector a) -> SpVector a -> BICGSTAB a -> BICGSTAB a
bicgstabStep aa r0hat (BICGSTAB x r p) = BICGSTAB xj1 rj1 pj1 where
     aap = aa #> p
     alphaj = (r <.> r0hat) / (aap <.> r0hat)
     sj = r ^-^ (alphaj .* aap)
     aasj = aa #> sj
     omegaj = (aasj <.> sj) / (aasj <.> aasj)
     xj1 = x ^+^ (alphaj .* p) ^+^ (omegaj .* sj)    -- updated solution
     rj1 = sj ^-^ (omegaj .* aasj)
     betaj = (rj1 <.> r0hat)/(r <.> r0hat) * alphaj / omegaj
     pj1 = rj1 ^+^ (betaj .* (p ^-^ (omegaj .* aap)))

instance Show a => Show (BICGSTAB a) where
  show (BICGSTAB x r p) = "x = " ++ show x ++ "\n" ++
                          "r = " ++ show r ++ "\n" ++
                          "p = " ++ show p ++ "\n"

bicgsInit aa b x0 = BICGSTAB x0 r0 r0 where
  r0 = b ^-^ (aa #> x0)   -- residual of initial guess solution











-- * Moore-Penrose pseudoinverse
-- | Least-squares approximation of a rectangular system of equaitons. Uses <\\> for the linear solve
-- pinv :: (Epsilon a, RealFloat a, Elt a, AdditiveGroup a) => SpMatrix a -> SpVector a -> SpVector a
-- pinv :: Epsilon a => SpMatrix a -> SpVector a -> Either LinSysError (SpVector a)
pinv aa b = aa #~^# aa <\> atb where
  atb = transpose aa #> b




-- linSolve0
--   :: (MonadThrow m,
--       Typeable (Magnitude (SpVector a)), Typeable a,
--       Show a, Fractional (Scalar (SpVector a))) =>
--      LinSolveMethod
--      -> MatrixType (SpVector a)
--      -> SpVector a
--      -> SpVector a
--      -> SpVector a
--      -> m (SpVector a)
linSolve0 method aa b x0 r0hat
  | n /= nb = throwM (MatVecSizeMismatchException "linSolve0" dm nb)
  | otherwise = solve aa b where
     solve aa' b' | isDiagonalSM aa' = return $ reciprocal aa' #> b' -- diagonal solve
                  | otherwise = xHat
     xHat = case method of
       BICGSTAB_ -> solver "BICGSTAB" nitermax _xBicgstab (bicgstabStep aa r0hat) (bicgsInit aa b x0)
       CGS_ -> solver "CGS" nitermax _x  (cgsStep aa r0hat) (cgsInit aa b x0)
       -- GMRES_ -> return $ gmres aa b
     nitermax = 200
     dm@(m,n) = dim aa
     nb = dim b


solver :: (Normed b, MonadThrow m, Typeable (Magnitude b), Typeable a,
      Show (Magnitude b), Show a, Ord (Magnitude b)) =>
     String -> Int -> (a -> b) -> (a -> a) -> a -> m b
solver fname nitermax fproj stepf initf = do
  xf <- untilConvergedG fname nitermax fproj (const True) stepf initf
  return $ fproj xf

  

-- * Linear solver interface

data LinSolveMethod = GMRES_ | CGNE_ | TFQMR_ | BCG_ | CGS_ | BICGSTAB_ deriving (Eq, Show) 


-- linSolve0 method aa b x0
--   | n /= nb = error "linSolve : operand dimensions mismatch"
--   | otherwise = solve aa b where
--       solve aa' b' | isDiagonalSM aa' = Right $ reciprocal aa' #> b' -- diagonal solve
--                    | otherwise = solnE
--       solnE | nearZero (norm2 ((aa #> xHat) ^-^ b)) = Right xHat
--             | otherwise = Left (NotConverged "linSolve0" nits xHat)
--       xHat = case method of
--         GMRES_ -> gmres aa b
--         CGNE_ -> _xCgne (cgne aa b x0)
--         TFQMR_ -> _xTfq (tfqmr aa b x0)
--         BCG_ -> _xBcg (bcg aa b x0)
--         BICGSTAB_ ->  _xBicgstab (bicgstab aa b x0 x0)
--         CGS_ -> _x (cgs aa b x0 x0)
--       (m, n) = dim aa
--       nb     = dim b

-- -- -- | linSolve using the GMRES method as default
-- instance LinearSystem (SpVector Double) where
--   aa <\> b = linSolve0 GMRES_ aa b (mkSpVR n $ replicate n 0.1)
--     where n = ncols aa

-- instance LinearSystem (SpVector (Complex Double)) where
--   aa <\> b = linSolve0 GMRES_ aa b (mkSpVC n $ replicate n 0.1)
--     where n = ncols aa












-- | TODO : if system is poorly conditioned, is it better to warn the user or just switch solvers (e.g. via the pseudoinverse) ?

-- linSolveQR aa b init f1 stepf
--   | isInfinite k = do
--        tell "linSolveQR : rank-deficient system"
--   | otherwise = do
--        solv aa b init
--     where
--      (q, r) = qr aa
--      k = conditionNumberSM r
--      solv aa b init = execState (untilConverged f1 stepf) init




















-- * Control primitives for bounded iteration with convergence check

-- | transform state until a condition is met
modifyUntil :: MonadState s m => (s -> Bool) -> (s -> s) -> m s
modifyUntil q f = do
  x <- get
  let y = f x
  put y
  if q y then return y
         else modifyUntil q f

modifyUntilT q f = do
  x <- MTS.get
  let y = f x
  MTS.put y
  if q y then return y
         else modifyUntilT q f     

modifyUntilM q f = do
  x <- get
  y <- f x
  put y
  if q y then return y
         else modifyUntilM q f   

-- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached  
loopUntilAcc :: Int -> ([t] -> Bool) -> (t -> t)  -> t -> t
loopUntilAcc nitermax q f x = go 0 [] x where
  go i ll xx | length ll < 2 = go (i + 1) (y : ll) y 
             | otherwise = if q ll || i == nitermax
                           then xx
                           else go (i + 1) (take 2 $ y : ll) y
                where y = f xx


-- | iterate until convergence is verified or we run out of a fixed iteration budget
-- untilConverged :: (MonadState s m, Epsilon a) => (s -> SpVector a) -> (s -> s) -> m s
untilConverged :: (MonadState s m, Normed v, Epsilon (Magnitude v)) =>
     (s -> v) -> (s -> s) -> m s
untilConverged fproj = modifyInspectN 200 (normDiffConverged fproj)


-- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached (i.e. same thing as `loopUntilAcc` but this one runs in the State monad)
modifyInspectN ::
  MonadState s m =>
    Int ->           -- iteration budget
    ([s] -> Bool) -> -- convergence criterion
    (s -> s) ->      -- state stepping function
    m s
modifyInspectN nitermax q f 
  | nitermax > 0 = go 0 []
  | otherwise = error "modifyInspectN : n must be > 0" where
      go i ll = do
        x <- get
        let y = f x
        if length ll < 2
          then do put y
                  go (i + 1) (y : ll)
          else if q ll || i == nitermax
               then do put y
                       return y
               else do put y
                       go (i + 1) (take 2 $ y : ll)


-- | This function makes some default choices on the `modifyInspectGuarded` machinery: convergence is assessed using the squared L2 distance between consecutive states, and divergence is detected when this function is increasing between pairs of measurements.
untilConvergedG :: (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable s,
      Show (Magnitude v), Show s, Ord (Magnitude v)) =>
        String
     -> Int 
     -> (s -> v) 
     -> (s -> Bool) 
     -> (s -> s)               -- ^ state update _function_
     -> s 
     -> m s
untilConvergedG fname nitermax fproj qfinal =
  modifyInspectGuarded fname nitermax (convergf fproj) nearZero qdiverg qfinal

untilConvergedGM ::
  (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable t, Show t) =>
     String
     -> Int
     -> (t -> v)
     -> (t -> Bool)
     -> (t -> MTS.StateT t m t)  -- ^ state update _arrow_
     -> t
     -> m t
untilConvergedGM fname nitermax fproj qfinal =
  modifyInspectGuardedM fname nitermax (convergf fproj) nearZero qdiverg qfinal     

qdiverg :: Ord a => a -> a -> Bool
qdiverg = (>)

convergf :: Normed v => (t -> v) -> [t] -> Magnitude v
convergf fp [s1, s0] = norm2 (fp s1 ^-^ fp s0)
convergf _ _ = 1/0





-- | `untilConvergedG0` is a special case of `untilConvergedG` that assesses convergence based on the L2 distance to a known solution `xKnown`
untilConvergedG0 ::
  (Normed v, MonadThrow m, Typeable (Magnitude v), Typeable s, Show (Magnitude v), Show s, Ord (Magnitude v)) =>
        String 
     -> Int
     -> (s -> v)
     -> v        -- ^ Known solution
     -> (s -> s)
     -> s
     -> m s
untilConvergedG0 fname nitermax fproj xKnown =
  untilConvergedG fname nitermax fproj (\s -> nearZero (norm2 $ fproj s ^-^ xKnown))



-- | `modifyInspectGuarded` is a high-order abstraction of a numerical iterative process. It accumulates a rolling window of 3 states and compares a summary `q` of the latest 2 with that of the previous two in order to assess divergence (e.g. if `q latest2 > q prev2` then it). The process ends when either we hit an iteration budget or relative convergence is verified. The function then assesses the final state with a predicate `qfinal` (e.g. against a known solution; if this is not known, the user can just supply `const True`)
modifyInspectGuarded ::
  (MonadThrow m, Typeable s, Typeable a, Show s, Show a) =>
        String             -- ^ Calling function name
     -> Int                -- ^ Iteration budget
     -> ([s] -> a)         -- ^ State array projection
     -> (a -> Bool)        -- ^ Convergence criterion
     -> (a -> a -> Bool)   -- ^ Divergence criterion
     -> (s -> Bool)        -- ^ Final state evaluation
     -> (s -> s)           -- ^ State evolution
     -> s                  -- ^ Initial state
     -> m s                -- ^ Final state
modifyInspectGuarded fname nitermax q qconverg qdiverg qfinal f x0
  | nitermax > 0 = checkFinal 
  | otherwise = throwM (NonNegError fname nitermax)
  where
    checkFinal = do
      xfinal <- MTS.execStateT (go 0 []) x0
      if qfinal xfinal
        then return xfinal
        else throwM (NotConverged fname nitermax xfinal)
    go i ll = do
      x <- MTS.get
      let y = f x
      if length ll < 3
      then do MTS.put y
              go (i + 1) (y : ll) -- accumulate a l=3 rolling state window to observe
      else do
         let qi = q (init ll)     -- summary of latest 2 states
             qt = q (tail ll)     -- "       "  previous 2 states
         if qconverg qi           -- relative convergence  
         then do MTS.put y
                 return ()
         else if i == nitermax    -- end of iterations w/o convergence
              then do
                MTS.put y
                throwM (NotConverged fname nitermax y)
              else do
                if qdiverg qi qt  -- diverging
                then throwM (Diverging fname i qi qt)
                else do MTS.put y -- not diverging, keep iterating
                        go (i + 1) (take 3 $ y : ll)


-- | ", monadic version
modifyInspectGuardedM
  :: (MonadThrow m, Typeable s, Typeable a, Show s, Show a) =>
     String
     -> Int
     -> ([s] -> a)
     -> (a -> Bool)
     -> (a -> a -> Bool)
     -> (s -> Bool)
     -> (s -> MTS.StateT s m s)
     -> s
     -> m s
modifyInspectGuardedM fname nitermax q qconverg qdiverg qfinal f x0
  | nitermax > 0 = checkFinal 
  | otherwise = throwM (NonNegError fname nitermax)
  where
    checkFinal = do
      xfinal <- MTS.execStateT (go 0 []) x0
      if qfinal xfinal
        then return xfinal
        else throwM (NotConverged fname nitermax xfinal)
    go i ll = do
      x <- MTS.get
      y <- f x
      if length ll < 3
      then do MTS.put y
              go (i + 1) (y : ll) -- accumulate a l=3 rolling state window to observe
      else do
         let qi = q (init ll)     -- summary of latest 2 states
             qt = q (tail ll)     -- "       "  previous 2 states
         if qconverg qi           -- relative convergence  
         then do MTS.put y
                 return ()
         else if i == nitermax    -- end of iterations w/o convergence
              then do
                MTS.put y
                throwM (NotConverged fname nitermax y)
              else do
                if qdiverg qi qt  -- diverging
                then throwM (Diverging fname i qi qt)
                else do MTS.put y -- not diverging, keep iterating
                        go (i + 1) (take 3 $ y : ll)






-- helper functions for estimating convergence
-- meanl :: (Foldable t, Fractional a) => t a -> a
-- meanl xx = 1/fromIntegral (length xx) * sum xx

-- norm2l :: (Foldable t, Functor t, Floating a) => t a -> a
-- norm2l xx = sqrt $ sum (fmap (**2) xx)

-- | Squared difference of a 2-element list.
-- | NB: unsafe !
diffSqL :: Floating a => [a] -> a
diffSqL xx = (x1 - x2)**2 where [x1, x2] = [head xx, xx!!1]


-- | Relative tolerance :
-- relTol a b := ||a - b|| / (1 + min (||norm2 a||, ||norm2 b||))
relTol :: (Normed v, Ord (Magnitude v)) => v -> v -> Magnitude v
relTol a b = norm2 (a ^-^ b) / m where
  m = 1 + min (norm2 a) (norm2 b)


-- | convergence test
normDiff :: Normed v => (t -> v) -> [t] -> Magnitude v
normDiff fp [x1, x0] = norm2Sq (fp x0 ^-^ fp x1)

normDiffConverged :: (Normed v, Epsilon (Magnitude v)) => (s -> v) -> [s] -> Bool
normDiffConverged fp ll = nearZero (normDiff fp ll)






-- | run `niter` iterations and append the state `x` to a list `xs`, stop when either the `xs` satisfies a predicate `q` or when the counter reaches 0
runAppendN :: ([t] -> Bool) -> (t -> t) -> Int -> t -> [t]
runAppendN qq ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                          | otherwise = go qq ff niter x0 [] where
  go q f n z xs = 
    let x = f z in
    if n <= 0 || q xs then xs
                      else go q f (n-1) x (x : xs)

-- | ", NO convergence check 
runAppendN' :: (t -> t) -> Int -> t -> [t]
runAppendN' ff niter x0 | niter<0 = error "runAppendN : niter must be > 0"
                        | otherwise = go ff niter x0 [] where
  go f n z xs = 
    let x = f z in
    if n <= 0 then xs
              else go f (n-1) x (x : xs)

  












-- -- * Random arrays

-- randArray :: PrimMonad m => Int -> Double -> Double -> m [Double]
-- randArray n mu sig = do
--   g <- MWC.create
--   replicateM n (MWC.normal mu sig g)
  



-- -- * Random matrices and vectors

-- -- |Dense SpMatrix
-- randMat :: PrimMonad m => Int -> m (SpMatrix Double)
-- randMat n = do
--   g <- MWC.create
--   aav <- replicateM (n^2) (MWC.normal 0 1 g)
--   let ii_ = [0 .. n-1]
--       (ix_,iy_) = unzip $ concatMap (zip ii_ . replicate n) ii_
--   return $ fromListSM (n,n) $ zip3 ix_ iy_ aav

-- -- | Dense SpVector  
-- randVec :: PrimMonad m => Int -> m (SpVector Double)
-- randVec n = do
--   g <- MWC.create
--   bv <- replicateM n (MWC.normal 0 1 g)
--   let ii_ = [0..n-1]
--   return $ fromListSV n $ zip ii_ bv



-- -- | Sparse SpMatrix
-- randSpMat :: Int -> Int -> IO (SpMatrix Double)
-- randSpMat n nsp | nsp > n^2 = error "randSpMat : nsp must be < n^2 "
--                 | otherwise = do
--   g <- MWC.create
--   aav <- replicateM nsp (MWC.normal 0 1 g)
--   ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
--   jj <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
--   return $ fromListSM (n,n) $ zip3 ii jj aav

-- -- | Sparse SpVector
-- randSpVec :: Int -> Int -> IO (SpVector Double)
-- randSpVec n nsp | nsp > n = error "randSpVec : nsp must be < n"
--                 | otherwise = do
--   g <- MWC.create
--   aav <- replicateM nsp (MWC.normal 0 1 g)
--   ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
--   return $ fromListSV n $ zip ii aav


