{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts, TypeFamilies, MultiParamTypeClasses, FlexibleInstances  #-}
{-# language RankNTypes #-}
-- {-# language ApplicativeDo #-}
-- {-# OPTIONS_GHC -O2 -rtsopts -with-rtsopts=-K32m -prof#-}
{-|

This module exposes the user interface to the library.

-}

module Numeric.LinearAlgebra.Sparse
       (
         -- * Linear solvers
         -- ** Iterative methods
         (<\>),
         -- -- ** Preconditioners
         -- jacobiPre, ilu0Pre, mSsorPre,   
         -- *** Moore-Penrose pseudoinverse
         pinv,
         -- ** Direct methods
         -- luSolve,
         -- *** Forward substitution
         -- triLowerSolve,
         -- *** Backward substitution
         -- triUpperSolve,
         -- * Eigensolvers
         -- eigsQR,
         -- eigRayleigh,
         -- eigsArnoldi,
         -- * Matrix factorization algorithms
         -- ** QR
         -- qr, 
         -- ** LU
         lu,
         -- ** Cholesky
         -- chol, 
         -- ** Arnoldi iteration
         arnoldi,    
         -- * Utilities
         -- ** Givens' rotation
         givens,
         -- ** Condition number
         -- conditionNumberSM,
         -- ** Householder reflection
         hhRefl,
         -- -- * Householder bidiagonalization
         -- ** Matrix partitioning
         diagPartitions,
         -- -- * Random arrays
         -- randArray,
         -- -- * Random matrices and vectors
         -- randMat, randVec, 
         -- -- ** Sparse "
         -- randSpMat, randSpVec,
         -- * Creation and conversion of sparse data
         -- ** SpVector
         -- *** Sparse
         fromListSV, toListSV,
         -- *** Dense
         -- **** " from a list of entries
         vr, vc,
         -- **** " from/to a Vector of entries
         fromVector, toVectorDense,
         -- **** " having constant elements
         constv,
         -- ** SpMatrix
         fromListSM, toListSM,

         -- ** Packing and unpacking, rows or columns of a sparse matrix
         -- *** ", using lists as container
         fromRowsL, toRowsL,
         fromColsL, toColsL,
         -- *** ", using Vector as container
         fromRowsV, fromColsV,
         -- ** Block operations
         (-=-), (-||-), fromBlocksDiag,
         -- ** Special matrices
         eye, mkDiagonal, mkSubDiagonal,
         permutationSM, permutPairsSM,
         -- * Predicates
         isOrthogonalSM, isDiagonalSM,
         -- * Manipulation of sparse data
         filterSV, ifilterSV,
         -- * Sparsity-related predicates
         nearZero, nearOne, isNz,
         -- * Operators
         -- ** Scaling
         (.*), (./), 
         -- ** Inner product
         (<.>),
         -- ** Matrix-vector product
         (#>), (<#),
         -- ** Matrix-matrix product
         (##), (#^#), (##^),
         -- *** Sparsifying matrix-matrix product
         (#~#), (#~^#), (#~#^),
         -- ** Vector outer product
         (><),
         -- * Common operations
         dim, nnz, spy,
         -- ** Vector spaces
         cvx,
         -- *** Norms and normalization
         norm, norm2, norm2', normalize, normalize2, normalize2',
         norm1, hilbertDistSq,
         -- ** Matrix-related
         transpose, trace, normFrobenius,
         -- * Pretty-printing
         prd, prd0,
         -- * Iteration combinators
         -- untilConvergedG0,
         -- -- untilConvergedG, untilConvergedGM,
         -- modifyInspectGuarded, modifyInspectGuardedM,
         modifyInspectGuardedM, 
         -- IterationConfig (..),
         modifyUntil, modifyUntilM,
         -- * Internal
         -- linSolve0,
         LinSolveMethod(..),
         -- * Exceptions
         PartialFunctionError,InputError, OutOfBoundsIndexError,
         OperandSizeMismatch, MatrixException, IterationException
             
       )
       where


import Control.Exception.Common
import Control.Iterative
import Data.Sparse.Common

import Control.Monad.Catch 
import Data.Typeable

-- import Control.Applicative ((<|>))

import Control.Monad.State.Strict
import qualified Control.Monad.Trans.State  as MTS 
import Data.Complex

import qualified Data.Sparse.Internal.IntM as I

import Data.Maybe

import qualified Data.Vector as V



-- | A lumped constraint for numerical types
type Num' x = (Epsilon x, Elt x, Show x, Ord x, Typeable x)





-- * Matrix condition number

-- | Matrix condition number: computes the QR factorization and extracts the extremal eigenvalues from the R factor

-- conditionNumberSM m = do
--   (_, r) <- qr m
--   let
--    u = extractDiagDense r 
--    lmax = abs (maximum u)
--    lmin = abs (minimum u)
--    kappa = lmax / lmin                     
--   if nearZero lmin
--   then throwM (HugeConditionNumber "conditionNumberSM" kappa)
--   else return kappa 
                          





-- * Householder transformation

hhMat :: (Num a, AdditiveGroup a) => a -> SpVector a -> SpMatrix a
hhMat beta x = eye n ^-^ beta `scale` (x >< x) where
  n = dim x


-- | Householder reflection: a vector `x` uniquely defines an orthogonal (hyper)plane, i.e. an orthogonal subspace; the Householder operator reflects any point `v` through this subspace: v' = (I - 2 x >< x) v
hhRefl :: (Num a, AdditiveGroup a) => SpVector a -> SpMatrix a
hhRefl = hhMat 2











-- * Givens rotation matrix

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



                  

-- | Givens matrix : a planar rotation embedded in R^n. 
--
-- >>> aa = fromListSM 
-- >>> g <- givens aa 1 0
-- 
-- Row version of the method: given a matrix element below the diagonal, indexed (i,j), choose a row index i' that is below the diagonal as well and distinct from i such that the corresponding element is nonzero.
--
-- To zero out entry A(i, j) we must find row i' such that A(i', j) is non-zero but A has zeros in row i' for all column indices < j.
-- 
-- NB: The Givens' matrix differs from Identity in 4 entries
-- 
-- NB2: The form of a Complex rotation matrix in R^2 is as follows (@*@ indicates complex conjugation):
--
-- @
--     ( c    s )
--  G =(        )
--     ( -s*  c*)
-- @
{-# inline givens #-}
givens :: (Elt a, MonadThrow m) => SpMatrix a -> IxRow -> IxCol -> m (SpMatrix a)
givens aa i j 
  | isValidIxSM aa (i,j) && nrows aa >= ncols aa = do
      i' <- candidateRows' (immSM aa) i j
      return $ givensMat aa i i' j
  | otherwise = throwM (OOBIxsError "Givens" [i, j])
  where
  givensMat mm i i' j =
    fromListSM'
      [(i,i, conj c), (i,j, - conj s),
       (j,i, s),      (j,j, c)]
      (eye (nrows mm))
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
    

-- | Givens coefficients and norm of associated vector
givensCoef :: Elt t => t -> t -> (t, t, t)
givensCoef u v = (c0/r, s0/r, r) where
  c0 = conj u
  s0 = conj v
  r = hypot u v

hypot :: Elt a => a -> a -> a
hypot x y = sqrt (mag2 x + mag2 y) where
    mag2 i = i * conj i





  


-- * QR decomposition

{-|
Given a matrix A, returns a pair of matrices (Q, R) such that Q R = A, where Q is orthogonal and R is upper triangular. Applies Givens rotation iteratively to zero out sub-diagonal elements.

NB: at each iteration `i` we multiply the Givens matrix `G_i` by the previous partial result `M`. Since this corresponds to a rotation, and the `givensCoef` function already computes the value of the resulting non-zero component (output `r`), `G_i ## M` can be simplified by just updating two entries of `M` (i.e. zeroing one out and changing the other into `r`).

However, we must also accumulate the `G_i` in order to build `Q`, and the present formulation follows this definition closely.
-}

-- {-# inline qr #-}
-- qr :: (Elt a, MatrixRing (SpMatrix a), PrintDense (SpMatrix a),
--       Epsilon a, MonadThrow m, MonadLog String m) =>
--      SpMatrix a
--      -> m (SpMatrix a, SpMatrix a) -- ^ Q, R
-- qr mm = do 
--      (qt, r, _) <- modifyUntilM' config haltf qrstepf gminit
--      return (transpose qt, r) 
--   where
--     gminit = (eye (nrows mm), mm, subdiagIndicesSM mm)
--     haltf (_, _, iis) = null iis
--     config = IterConf 0 False fst2 prd2 where
--       fst2 (x,y,_) = (x,y)
--       prd2 (x,y) = do
--         prd0 x
--         prd0 y
--     qrstepf (qmatt, m, iis) = do
--         let (i, j) = head iis
--         g <- givens m i j
--         let
--           qmatt' = g #~# qmatt  -- update Q'
--           m' = g #~# m          -- update R
--         return (qmatt', m', tail iis)    







-- * Eigenvalue algorithms

-- ** QR algorithm

-- | @eigsQR n mm@ performs at most @n@ iterations of the QR algorithm on matrix @mm@, and returns a 'SpVector' containing all eigenvalues.

-- eigsQR nitermax debq m = pf <$> untilConvergedGM "eigsQR" c (const True) stepf m
--   where
--     pf = extractDiagDense
--     c = IterConf nitermax debq pf prd
--     stepf mm = do
--       (q, _) <- qr mm
--       return $ q #~^# (m ## q) -- r #~# q





-- ** Arnoldi-QR

-- | @eigsArnoldi n aa b@ computes at most n iterations of the Arnoldi algorithm to find a Krylov subspace of (A, b), denoted Q, along with a Hessenberg matrix of coefficients H.
--
-- After that, it computes the QR decomposition of H, denoted (O, R) and the eigenvalues {λ_i} of A are listed on the diagonal of the R factor.

-- eigsArnoldi :: (Scalar (SpVector t) ~ t, MatrixType (SpVector t) ~ SpMatrix t,
--       Elt t, V (SpVector t), Epsilon t, PrintDense (SpMatrix t),
--       MatrixRing (SpMatrix t), MonadThrow m, MonadLog String m) =>
--      Int
--      -> SpMatrix t
--      -> SpVector t
--      -> m (SpMatrix t, SpMatrix t, SpVector t) -- ^ Q, O, {λ_i}
-- eigsArnoldi nitermax aa b = do
--   (q, h) <- arnoldi aa b nitermax
--   (o, r) <- qr h
--   return (q, o, extractDiagDense r)
  



-- * Householder vector 

-- (Golub & Van Loan, Alg. 5.1.1, function `house`)
hhV :: (Scalar (SpVector t) ~ t, Elt t, InnerSpace (SpVector t), Epsilon t) =>
     SpVector t -> (SpVector t, t)
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

                         



-- -- * Bidiagonalization



-- -- * SVD

{- Golub & Van Loan, sec 8.6.2 (p 452 segg.)

SVD of A, Golub-Kahan method

* reduce A to upper bidiagonal form B (Alg. 5.4.2, Householder bidiagonalization)
* compute SVD of B (implicit-shift QR step applied to B^T B, Alg. 8.3.2)

-}








-- * Cholesky factorization

-- | Given a positive semidefinite matrix A, returns a lower-triangular matrix L such that L L^T = A . This is an implementation of the Cholesky–Banachiewicz algorithm, i.e. proceeding row by row from the upper-left corner.
-- | NB: The algorithm throws an exception if some diagonal element of A is zero.

-- chol :: (Elt a, Epsilon a, MonadThrow m, MonadLog String m, PrintDense (SpMatrix a)) =>
--         SpMatrix a
--      -> m (SpMatrix a)  -- ^ L
-- chol aa = do
--   let n = nrows aa
--       q (i, _) = i == n
--       config = IterConf 0 False snd prd0
--   l0 <- cholUpd aa (0, zeroSM n n)
--   (_, lfin) <- modifyUntilM' config q (cholUpd aa) l0
--   return lfin
--   where
--    oops i = throwM (NeedsPivoting "chol" (unwords ["L", show (i,i)]) :: MatrixException Double)
--    cholUpd aa (i, ll) = do 
--      sd <- cholSDRowUpd aa ll i  -- update subdiagonal entries
--      ll' <- cholDiagUpd aa sd i  -- update diagonal entries
--      return (i + 1, ll')
--    cholSDRowUpd aa ll i = do 
--      lrs <- fromListSV (i + 1) <$> onRangeSparseM cholSubDiag [0 .. i-1]
--      return $ insertRow ll lrs i where
--        cholSubDiag j | isNz ljj = return $ 1/ljj*(aij - inn)
--                      | otherwise = oops j
--         where
--           ljj = ll @@! (j, j)
--           aij = aa @@! (i, j)
--           inn = contractSub ll ll i j (j - 1)
--    cholDiagUpd aa ll i = do
--      cd <- cholDiag 
--      return $ insertSpMatrix i i cd ll where
--        cholDiag | i == 0 = sqrt <$> aai
--                 | otherwise = do
--                     a <- aai
--                     let l = sum (fmap (**2) lrow)
--                     return $ sqrt (a - l)
--          where
--           lrow = ifilterSV (\j _ -> j < i) (extractRow ll i) -- sub-diagonal elems of L
--           aai | isNz aaii = return aaii
--               | otherwise = oops i
--             where
--              aaii = aa @@! (i,i)
                 
















-- * LU factorization

-- | Given a matrix A, returns a pair of matrices (L, U) where L is lower triangular and U is upper triangular such that L U = A .
-- 
-- Implements the Doolittle algorithm, which sets the diagonal of the L matrix to ones and expects all diagonal entries of A to be nonzero. Apply pivoting (row or column permutation) to enforce a nonzero diagonal of the A matrix (the algorithm throws an appropriate exception otherwise).
lu :: (Scalar (SpVector t) ~ t, Elt t, VectorSpace (SpVector t), Epsilon t,
        MonadThrow m) =>
     SpMatrix t
     -> m (SpMatrix t, SpMatrix t) -- ^ L, U
lu aa = do
  let oops j = throwM (NeedsPivoting "solveForLij" ("U" ++ show (j, j)) :: MatrixException Double)
      n = nrows aa
      q (i, _, _) = i == n - 1
      luInit | isNz u00 = return (1, l0, u0)
             | otherwise = oops (0 :: Int)
        where
          l0 = insertCol (eye n) (extractSubCol aa 0 (1, n-1) ./ u00 ) 0
          u0 = insertRow (zeroSM n n) (extractRow aa 0) 0   -- initial U
          u00 = u0 @@! (0,0)  -- make sure this is non-zero by applying permutation
      luUpd (i, l, u) = do -- (i + 1, l', u') 
        u' <- uUpd aa n (i, l, u)  -- update U
        l' <- lUpd (i, l, u') -- update L
        return (i + 1, l', u')
      uUpd aa n (ix, lmat, umat) = do
        let us = onRangeSparse (solveForUij ix) [ix .. n - 1]
            solveForUij i j = a - p where
              a = aa @@! (i, j)
              p = contractSub lmat umat i j (i - 1)
        return $ insertRow umat (fromListSV n us) ix
      lUpd (ix, lmat, umat) = do -- insertCol lmat (fromListSV n ls) ix
        ls <- lsm
        return $ insertCol lmat (fromListSV n ls) ix
        where
          lsm = onRangeSparseM (`solveForLij` ix) [ix + 1 .. n - 1]
          solveForLij i j
            | isNz ujj = return $ (a - p)/ujj
            | otherwise = oops j
            where
             a = aa @@! (i, j)
             ujj = umat @@! (j , j)   -- NB this must be /= 0
             p = contractSub lmat umat i j (i - 1)    
  s0 <- luInit
  (ixf, lf, uf) <- MTS.execStateT (modifyUntilM q luUpd) s0
  ufin <- uUpd aa n (ixf, lf, uf)   -- final U update
  return (lf, ufin)
    





-- lu' aa = do
--   let oops j = throwM (NeedsPivoting "solveForLij" ("U" ++ show (j, j)) :: MatrixException Double)
--       n = nrows aa
--       q (i, _, _) = i == n - 1
--       luInit | isNz u00 = return (1, l0, u0)
--              | otherwise = oops (0 :: Int)
--         where
--           l0 = insertCol (eye n) ((extractSubCol aa 0 (1, n-1)) ./ u00 ) 0
--           u0 = insertRow (zeroSM n n) (extractRow aa 0) 0   -- initial U
--           u00 = u0 @@! (0,0)  -- make sure this is non-zero by applying permutation
--       luUpd (i, l, u) = do -- (i + 1, l', u') 
--         u' <- uUpd aa n (i, l, u)  -- update U
--         l' <- lUpd (i, l, u') -- update L
--         return (i + 1, l', u')
--       uUpd aa n (ix, lmat, umat) = do
--         let us = onRangeSparse (solveForUij ix) [ix .. n - 1]
--             solveForUij i j = a - p where
--               a = aa @@! (i, j)
--               p = contractSub lmat umat i j (i - 1)
--         return $ insertRow umat (fromListSV n us) ix
--       lUpd (ix, lmat, umat) = do -- insertCol lmat (fromListSV n ls) ix
--         ls <- lsm
--         return $ insertCol lmat (fromListSV n ls) ix
--         where
--           lsm = onRangeSparseM (`solveForLij` ix) [ix + 1 .. n - 1]
--           solveForLij i j
--             | isNz ujj = return $ (a - p)/ujj
--             | otherwise = oops j
--             where
--              a = aa @@! (i, j)
--              ujj = umat @@! (j , j)   -- NB this must be /= 0
--              p = contractSub lmat umat i j (i - 1)    
--   s0 <- luInit
--   let config = IterConf 0 True vf prd2 where
--         vf (_, l, u) = (l, u)
--         prd2 (x, y) = do
--           prd0 x
--           prd0 y
--   (ixf, lf, uf) <- modifyUntilM' config q luUpd s0
--   ufin <- uUpd aa n (ixf, lf, uf)   -- final U update
--   return (lf, ufin)
         


tmc4, tmc5, tmc6 :: SpMatrix (Complex Double)
tmc4 = fromListDenseSM 3 [3:+1, 4:+(-1), (-5):+3, 2:+2, 3:+(-2), 5:+0.2, 7:+(-2), 9:+(-1), 2:+3]

tvc4 = vc [1:+3,2:+2,1:+9]

tmc5 = fromListDenseSM 4 $ zipWith (:+) [16..31] [17,14..]


tmc6 = fromListDenseSM 2 $ zipWith (:+) [1,2,3,4] [4,3,2,1] 



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
-- 
-- At the i`th iteration, it finds (i + 1) coefficients (the i`th column of the Hessenberg matrix H) and the (i + 1)`th Krylov vector.
arnoldi :: (MatrixType (SpVector a) ~ SpMatrix a, V (SpVector a) ,
            Scalar (SpVector a) ~ a, Epsilon a, MonadThrow m) =>
     SpMatrix a                    -- ^ System matrix
     -> SpVector a                 -- ^ r.h.s.
     -> Int                        -- ^ Max. # of iterations
     -> m (SpMatrix a, SpMatrix a) -- ^ Q, H
arnoldi aa b kn | n == nb = return (fromColsV qvfin, fromListSM (nmax + 1, nmax) hhfin)
                | otherwise = throwM (MatVecSizeMismatchException "arnoldi" (m,n) nb)
  where
  (qvfin, hhfin, nmax, _) = execState (modifyUntil tf arnoldiStep) arnInit 
  tf (_, _, ii, fbreak) = ii == kn || fbreak -- termination criterion
  (m, n) = (nrows aa, ncols aa)
  nb = dim b
  arnInit = (qv1, hh1, 1, False) where
      q0 = normalize2 b   -- starting basis vector
      aq0 = aa #> q0       -- A q0
      h11 = q0 `dot` aq0          
      q1nn = aq0 ^-^ (h11 .* q0)
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
      aqi ^-^ V.foldl' (^+^) zv (V.zipWith (.*) hhcoli qv) -- unnormalized q_{i+1}
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
diagPartitions :: SpMatrix a
               -> (SpMatrix a, SpMatrix a, SpMatrix a) -- ^ Subdiagonal, diagonal, superdiagonal partitions
diagPartitions aa = (e,d,f) where
  e = extractSubDiag aa
  d = extractDiag aa
  f = extractSuperDiag aa


-- ** Jacobi preconditioner

-- | The Jacobi preconditioner is just the reciprocal of the diagonal 
jacobiPre :: Fractional a => SpMatrix a -> SpMatrix a
jacobiPre x = recip <$> extractDiag x


-- ** Incomplete LU

-- | Used for Incomplete LU : remove entries in the output matrix corresponding to zero entries in the input matrix (this is called ILU(0) in the preconditioner literature)
ilu0Pre :: (Scalar (SpVector t) ~ t, Elt t, VectorSpace (SpVector t),
      Epsilon t, MonadThrow m) =>
     SpMatrix t
     -> m (SpMatrix t, SpMatrix t) -- ^ L, U (with holes)
ilu0Pre aa = do
  (l, u) <- lu aa
  let lh = sparsifyLU l aa
      uh = sparsifyLU u aa
      sparsifyLU m m2 = ifilterSM f m where
        f i j _ = isJust (lookupSM m2 i j)
  return (lh, uh)


-- ** SSOR

-- | Symmetric Successive Over-Relaxation. `mSsor aa omega` : if `omega = 1` it returns the symmetric Gauss-Seidel preconditioner. When ω = 1, the SOR reduces to Gauss-Seidel; when ω > 1 and ω < 1, it corresponds to over-relaxation and under-relaxation, respectively.
mSsorPre :: (MatrixRing (SpMatrix b), Fractional b) =>
     SpMatrix b
     -> b  -- ^ relaxation factor
     -> (SpMatrix b, SpMatrix b) -- ^ Left, right factors
mSsorPre aa omega = (l, r) where
  (e, d, f) = diagPartitions aa
  n = nrows e
  l = (eye n ^-^ scale omega e) ## reciprocal d
  r = d ^-^ scale omega f 








-- * Linear solver, LU-based

-- luSolveConfig :: PrintDense (SpVector t) => IterationConfig (SpVector t, IxRow) (SpVector t)
-- luSolveConfig = IterConf 0 False fst prd0



-- -- | Direct solver based on a triangular factorization of the system matrix.
-- luSolve :: (Scalar (SpVector t) ~ t, MonadThrow m, Elt t, InnerSpace (SpVector t),
--             Epsilon t, PrintDense (SpVector t), MonadLog String m) =>
--      SpMatrix t    -- ^ Lower triangular
--      -> SpMatrix t -- ^ Upper triangular
--      -> SpVector t -- ^ r.h.s.
--      -> m (SpVector t)
-- luSolve ll uu b
--   | isLowerTriSM ll && isUpperTriSM uu = do
--       w <- triLowerSolve0 luSolveConfig ll b
--       triUpperSolve0 luSolveConfig uu w
--   | otherwise = throwM (NonTriangularException "luSolve")

-- | Forward substitution solver

-- triLowerSolve = triLowerSolve0 luSolveConfig




-- triLowerSolve0 config ll b = do
--   let q (_, i) = i == nb
--       nb = svDim b
--       oops i = throwM (NeedsPivoting "triLowerSolve" (unwords ["L", show (i, i)]) :: MatrixException Double)
--       lStep (ww, i) = do -- (ww', i + 1) where
--         let
--           lii = ll @@ (i, i)
--           bi = b @@ i
--           wi = (bi - r)/lii where
--             r = extractSubRow ll i (0, i-1) `dot` takeSV i ww
--         if isNz lii
--           then return (insertSpVector i wi ww, i + 1)
--           else oops i
--       lInit = do -- (ww0, 1) where
--         let
--           l00 = ll @@ (0, 0)
--           b0 = b @@ 0
--           w0 = b0 / l00
--         if isNz l00
--           then return (insertSpVector 0 w0 $ zeroSV (dim b), 1)
--           else oops (0 :: Int)
--   l0 <- lInit             
--   (v, _) <- modifyUntilM' config q lStep l0
--   return $ sparsifySV v



-- NB in the computation of `xi` we must rebalance the subrow indices (extractSubRow_RK) because `dropSV` does that too, in order to take the inner product with consistent index pairs
-- | Backward substitution solver

-- triUpperSolve = triUpperSolve0 luSolveConfig

-- triUpperSolve0 conf uu w = do 
--   let q (_, i) = i == (- 1)
--       nw = svDim w
--       oops i = throwM (NeedsPivoting "triUpperSolve" (unwords ["U", show (i, i)]) :: MatrixException Double)
--       uStep (xx, i) = do
--         let uii = uu @@ (i, i) 
--             wi = w @@ i
--             r = extractSubRow_RK uu i (i + 1, nw - 1) `dot` dropSV (i + 1) xx  
--             xi = (wi - r) / uii
--         if isNz uii
--           then return (insertSpVector i xi xx, i - 1)
--           else oops i 
--       uInit = do 
--         let i = nw - 1
--             u00 = uu @@! (i, i)
--             w0 = w @@ i
--             x0 = w0 / u00
--         if isNz u00
--           then return (insertSpVector i x0 (zeroSV nw), i - 1)
--           else oops (0 :: Int)
--   u0 <- uInit             
--   (x, _) <- modifyUntilM' conf q uStep u0
--   return $ sparsifySV x      





    





-- * Iterative linear solvers


-- ** GMRES

-- | Given a linear system `A x = b` where `A` is an (m x m) real-valued matrix, the GMRES method finds an approximate solution `xhat` such that the Euclidean norm of the residual `A xhat - b` is minimized. `xhat` is spanned by the order-`n` Krylov subspace of (A, b).
-- In this implementation:
-- 1) the Arnoldi iteration is carried out until numerical breakdown (therefore yielding _at_most_ `m+1` Krylov basis vectors)
-- 2) the resulting Hessenberg matrix H is factorized in QR form (H = Q R)
-- 3) the Krylov-subspace solution `yhat` is found by backsubstitution (since R is upper-triangular)
-- 4) the approximate solution in the original space `xhat` is computed using the Krylov basis, `xhat = Q_n yhat`
--
-- A common optimization involves interleaving the QR factorization (and the subsequent triangular solve) with the Arnoldi process (and employing an updating QR factorization which only requires one Givens' rotation at every update). 

-- gmres aa b = do
--   let m = ncols aa
--   (qa, ha) <- arnoldi aa b m   -- at most m steps of Arnoldi (aa, b)
--     -- b' = (transposeSe qa) #> b
--   let b' = norm2' b .* ei mp1 1  -- b rotated back to canonical basis by Q^T
--         where mp1 = nrows ha     -- = 1 + (# Arnoldi iterations)
--   (qh, rh) <- qr ha
--   let rhs' = takeSV (dim b' - 1) (transpose qh #> b')
--       rh' = takeRows (nrows rh - 1) rh -- last row of `rh` is 0
--   yhat <- triUpperSolve rh' rhs'
--   let qa' = takeCols (ncols qa - 1) qa  -- we don't use last column of Krylov basis   
--   return $ qa' #> yhat



  
-- ** CGNE

data CGNE a =
  CGNE {_xCgne , _rCgne, _pCgne :: SpVector a} deriving Eq
instance Show a => Show (CGNE a) where
    show (CGNE x r p) = "x = " ++ show x ++ "\n" ++
                        "r = " ++ show r ++ "\n" ++
                        "p = " ++ show p ++ "\n"

cgneInit :: (MatrixType (SpVector a) ~ SpMatrix a,
      LinearVectorSpace (SpVector a)) =>
     SpMatrix a -> SpVector a -> SpVector a -> CGNE a
cgneInit aa b x0 = CGNE x0 r0 p0 where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = transposeSM aa #> r0

cgneStep :: (MatrixType (SpVector a) ~ SpMatrix a,
      LinearVectorSpace (SpVector a), InnerSpace (SpVector a),
      MatrixRing (SpMatrix a), Fractional (Scalar (SpVector a))) =>
     SpMatrix a -> CGNE a -> CGNE a
cgneStep aa (CGNE x r p) = CGNE x1 r1 p1 where
    alphai = (r `dot` r) / (p `dot` p)
    x1 = x ^+^ (alphai .* p)
    r1 = r ^-^ (alphai .* (aa #> p))
    beta = (r1 `dot` r1) / (r `dot` r)
    p1 = transpose aa #> r ^+^ (beta .* p)





-- ** BCG

data BCG a =
  BCG { _xBcg, _rBcg, _rHatBcg, _pBcg, _pHatBcg :: SpVector a } deriving Eq

bcgInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> BCG a
bcgInit aa b x0 = BCG x0 r0 r0hat p0 p0hat where
  r0 = b ^-^ (aa #> x0)    
  r0hat = r0
  p0 = r0
  p0hat = r0

bcgStep :: (MatrixType (SpVector a) ~ SpMatrix a,
      LinearVectorSpace (SpVector a), InnerSpace (SpVector a),
      MatrixRing (SpMatrix a), Fractional (Scalar (SpVector a))) =>
     SpMatrix a -> BCG a -> BCG a 
bcgStep aa (BCG x r rhat p phat) = BCG x1 r1 rhat1 p1 phat1 where
    aap = aa #> p
    alpha = (r `dot` rhat) / (aap `dot` phat)
    x1 = x ^+^ (alpha .* p)
    r1 = r ^-^ (alpha .* aap)
    rhat1 = rhat ^-^ (alpha .* (transpose aa #> phat))
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

cgsInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> CGS a
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

data BICGSTAB a =
  BICGSTAB { _xBicgstab, _rBicgstab, _pBicgstab :: SpVector a} deriving Eq

bicgsInit :: LinearVectorSpace (SpVector a) =>
     MatrixType (SpVector a) -> SpVector a -> SpVector a -> BICGSTAB a
bicgsInit aa b x0 = BICGSTAB x0 r0 r0 where
  r0 = b ^-^ (aa #> x0)   -- residual of initial guess solution

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






-- * Moore-Penrose pseudoinverse
-- | Least-squares approximation of a rectangular system of equations.
-- pinv :: (LinearSystem v, MatrixRing (MatrixType v), MonadThrow m, MonadIO m) =>
--      MatrixType v -> v -> m v
pinv aa b = (aa #^# aa) <\> atb where
  atb = transpose aa #> b




-- * Linear solver interface

-- -- TFQMR is in LinearSolvers.Experimental for now
-- | Iterative methods for linear systems
data LinSolveMethod = GMRES_  -- ^ Generalized Minimal RESidual 
                    | CGNE_   -- ^ Conjugate Gradient on the Normal Equations
                    | BCG_    -- ^ BiConjugate Gradient
                    | CGS_    -- ^ Conjugate Gradient Squared
                    | BICGSTAB_ -- ^ BiConjugate Gradient Stabilized
                    deriving (Eq, Show)

-- -- | Interface method to individual linear solvers
-- linSolve0 fh flog method aa b x0
--   | m /= nb = throwM (MatVecSizeMismatchException "linSolve0" dm nb)
--   | otherwise = solve aa b where
--      solve aa' b' | isDiagonalSM aa' =
--                         return $ reciprocal aa' #> b' -- diagonal solve
--                   | otherwise = xHat
--      xHat = case method of
--        -- BICGSTAB_ -> solver "BICGSTAB" nits _xBicgstab (bicgstabStep aa r0hat) (bicgsInit aa b x0)
--        -- BCG_ -> solver "BCG" nits _xBcg (bcgStep aa) (bcgInit aa b x0)
--        CGS_ -> solver "CGS" nits _x  (cgsStep aa r0hat) (cgsInit aa b x0)
--        -- GMRES_ -> gmres aa b          
--        CGNE_ -> solver "CGNE" nits _xCgne (cgneStep aa) (cgneInit aa b x0)
--      r0hat = b ^-^ (aa #> x0)
--      nits = 200
--      dm@(m, _) = dim aa
--      nb = dim b
--      lwindow = 3
--      solver fname nitermax fproj stepf initf =
--        solver' fname fh flog nitermax lwindow fproj stepf initf



-- solver' name fh flog nitermax lwindow fproj stepf initf = do
--   xf <- untilConvergedG fh name nitermax lwindow fproj (flog . fproj) stepf initf
--   return $ fproj xf    

class IterativeSolver s where
  -- solver :: 
  



-- -- | <\> uses the GMRES method as default

-- instance LinearSystem (SpVector Double) where
--   aa <\> b = linSolve0 GMRES_ aa b (mkSpVR n $ replicate n 0.1)
--     where n = ncols aa

-- -- instance LinearSystem (SpVector (Complex Double)) where
-- --   aa <\> b = linSolve0 GMRES_ aa b (mkSpVC n $ replicate n 0.1)
-- --     where n = ncols aa











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



  





-- test data

-- -- aa4 : eigenvalues 1 (mult.=2) and -1
-- aa4 :: SpMatrix Double
-- aa4 = fromListDenseSM 3 [3,2,-2,2,2,-1,6,5,-4] 

-- aa4c :: SpMatrix (Complex Double)
-- aa4c = toC <$> aa4

