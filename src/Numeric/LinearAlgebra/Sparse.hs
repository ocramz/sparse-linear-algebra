{-# LANGUAGE FlexibleContexts, TypeFamilies, MultiParamTypeClasses, FlexibleInstances  #-}
-- {-# OPTIONS_GHC -O2 -rtsopts -with-rtsopts=-K32m -prof#-}
module Numeric.LinearAlgebra.Sparse
       (
         -- * Matrix factorizations
         qr, lu,
         chol,
         -- * Condition number
         conditionNumberSM,
         -- * Householder reflection
         hhMat, hhRefl,
         -- * Givens' rotation
         givens,
         -- * Arnoldi iteration
         arnoldi,
         -- * Eigensolvers
         eigsQR, eigRayleigh,
         -- * Linear solvers
         -- ** Iterative methods
         linSolve, LinSolveMethod(..), (<\>),
         pinv,
         -- ** Direct methods
         luSolve,
         -- * Preconditioners
         ilu0, mSsor,
         -- * Matrix partitioning
         diagPartitions,
         -- * Random arrays
         randArray,
         -- * Random matrices and vectors
         randMat, randVec, 
         -- ** Sparse "
         randSpMat, randSpVec,
         -- * Sparsify data
         sparsifySV,
         -- * Iteration combinators
         modifyInspectN, runAppendN', untilConverged,
         diffSqL
       )
       where


import Data.Sparse.Common


import Control.Monad.Primitive
import Control.Monad (mapM_, forM_, replicateM)
import Control.Monad.State.Strict
import Control.Monad.Writer
-- import Control.Monad.Trans
-- import Control.Monad.Trans.State (runStateT)
-- import Control.Monad.Trans.Writer (runWriterT)

import qualified Data.IntMap.Strict as IM
-- import Data.Utils.StrictFold (foldlStrict) -- hidden in `containers`

import qualified System.Random.MWC as MWC
import qualified System.Random.MWC.Distributions as MWC

import Data.Monoid
import qualified Data.Foldable as F
import qualified Data.Traversable as T

-- import qualified Data.List as L
import Data.Maybe

import qualified Data.Vector as V



  
-- * Sparsify : remove almost-0 elements (|x| < eps)
-- | Sparsify an SpVector
sparsifySV :: Epsilon a => SpVector a -> SpVector a
sparsifySV = filterSV isNz





-- * Matrix condition number

-- |uses the R matrix from the QR factorization
conditionNumberSM :: (Epsilon a, RealFloat a) => SpMatrix a -> a
conditionNumberSM m | isInfinite kappa = error "Infinite condition number : rank-deficient system"
                    | otherwise = kappa where
  kappa = lmax / lmin
  (_, r) = qr m
  u = extractDiagDense r  -- FIXME : need to extract with default element 0 
  lmax = abs (maximum u)
  lmin = abs (minimum u)







-- * Householder transformation

hhMat :: Num a => a -> SpVector a -> SpMatrix a
hhMat beta x = eye n ^-^ scale beta (x >< x) where
  n = dim x


{-| a vector `x` uniquely defines an orthogonal plane; the Householder operator reflects any point `v` with respect to this plane:
 v' = (I - 2 x >< x) v
-}
hhRefl :: Num a => SpVector a -> SpMatrix a
hhRefl = hhMat (fromInteger 2)











-- * Givens rotation matrix


hypot :: Floating a => a -> a -> a
hypot x y = abs x * (sqrt (1 + y/x)**2)

sign :: (Ord a, Num a) => a -> a
sign x
  | x > 0 = 1
  | x == 0 = 0
  | otherwise = -1 

-- | Givens coefficients (using stable algorithm shown in  Anderson, Edward (4 December 2000). "Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem". LAPACK Working Note)
givensCoef :: (Ord a, Floating a) => a -> a -> (a, a, a)
givensCoef a b  -- returns (c, s, r) where r = norm (a, b)
  | b==0 = (sign a, 0, abs a)
  | a==0 = (0, sign b, abs b)
  | abs a > abs b = let t = b/a
                        u = sign a * abs ( sqrt (1 + t**2))
                      in (1/u, - t/u, a*u)
  | otherwise = let t = a/b
                    u = sign b * abs ( sqrt (1 + t**2))
                in (t/u, - 1/u, b*u)


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
givens :: (Floating a, Epsilon a, Ord a) => SpMatrix a -> IxRow -> IxCol -> SpMatrix a
givens mm i j 
  | isValidIxSM mm (i,j) && isSquareSM mm =
       sparsifySM $ fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye (nrows mm))
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    i' = head $ fromMaybe (error $ "givens: no compatible rows for entry " ++ show (i,j)) (candidateRows (immSM mm) i j)
    a = mm @@ (i', j)
    b = mm @@ (i, j)   -- element to zero out

-- |Returns a set of rows {k} that satisfy QR.C1
candidateRows :: IM.IntMap (IM.IntMap a) -> IxRow -> IxCol -> Maybe [IM.Key]
candidateRows mm i j | IM.null u = Nothing
                     | otherwise = Just (IM.keys u) where
  u = IM.filterWithKey (\irow row -> irow /= i &&
                                     firstNonZeroColumn row j) mm

-- |Is the `k`th the first nonzero column in the row?
{-# inline firstNonZeroColumn #-}
firstNonZeroColumn :: IM.IntMap a -> IxRow -> Bool
firstNonZeroColumn mm k = isJust (IM.lookup k mm) &&
                          isNothing (IM.lookupLT k mm)





-- * QR decomposition


-- | Given a matrix A, returns a pair of matrices (Q, R) such that Q R = A, where Q is orthogonal and R is upper triangular. Applies Givens rotation iteratively to zero out sub-diagonal elements.
qr :: (Epsilon a, Floating a, Real a) => SpMatrix a -> (SpMatrix a, SpMatrix a)
qr mm = (transposeSM qmatt, rmat)  where
  qmatt = F.foldl' (#~#) ee $ gmats mm -- Q^T = (G_n * G_n-1 ... * G_1)
  rmat = qmatt #~# mm                  -- R = Q^T A
  ee = eye (nrows mm)
      
-- | Givens matrices in order [G1, G2, .. , G_N ]
gmats :: (Epsilon a, Real a, Floating a) => SpMatrix a -> [SpMatrix a]
gmats mm = gm mm (subdiagIndicesSM mm) where
 gm m ((i,j):is) = let g = givens m i j
                   in g : gm (g #~# m) is
 gm _ [] = []





-- -- | QR algorithm, state transformer version
-- gmatST0 (m, (i,j):is) = (m', is) where    -- WRONG, possible access to []
--   g = givens m i j                        
--   m' = g #~# m
-- gmatST0 (m, []) = (eye (nrows m), [])

-- gmatST m = gmatST0 (m, subdiagIndicesSM m)






-- * Eigenvalue algorithms

-- ** QR algorithm

-- | `eigsQR n mm` performs `n` iterations of the QR algorithm on matrix `mm`, and returns a SpVector containing all eigenvalues
eigsQR :: (Epsilon a, Real a, Floating a) => Int -> SpMatrix a -> SpVector a
eigsQR nitermax m = extractDiagDense $ execState (convergtest eigsStep) m where
  eigsStep m = r #~# q where (q, r) = qr m
  convergtest g = modifyInspectN nitermax f g where
    f [m1, m2] = let dm1 = extractDiagDense m1
                     dm2 = extractDiagDense m2
                 in nearZero $ norm2 (dm1 ^-^ dm2)






-- ** Rayleigh iteration

-- | `eigsRayleigh n mm` performs `n` iterations of the Rayleigh algorithm on matrix `mm` and returns the eigenpair closest to the initialization. It displays cubic-order convergence, but it also requires an educated guess on the initial eigenpair
-- eigRayleigh :: Int                -- max # iterations
--      -> SpMatrix Double           -- matrix
--      -> (SpVector Double, Double) -- initial guess of (eigenvector, eigenvalue)
--      -> (SpVector Double, Double) -- final estimate of (eigenvector, eigenvalue)
eigRayleigh nitermax m = execState (convergtest (rayleighStep m)) where
  convergtest g = modifyInspectN nitermax f g where
    f [(b1, _), (b2, _)] = nearZero $ norm2 (b2 ^-^ b1)
  rayleighStep aa (b, mu) = (b', mu') where
      ii = eye (nrows aa)
      nom = (aa ^-^ (mu `matScale` ii)) <\> b
      b' = normalize 2 nom
      mu' = b' `dot` (aa #> b') / (b' `dot` b')




-- * Householder vector 

-- (Golub & Van Loan, Alg. 5.1.1, function `house`)
hhV :: (Epsilon a, Real a, Floating a) => SpVector a -> (SpVector a, a)
hhV x = (v, beta) where
  n = dim x
  tx = tailSV x
  sigma = tx `dot` tx
  vtemp = singletonSV 1 `concatSV` tx
  (v, beta) | nearZero sigma = (vtemp, 0)
            | otherwise = let mu = sqrt (headSV x**2 + sigma)
                              xh = headSV x
                              vh | xh <= 1 = xh - mu
                                 | otherwise = - sigma / (xh + mu)
                              vnew = (1 / vh) .* insertSpVector 0 vh vtemp     
                          in (vnew, 2 * xh**2 / (sigma + vh**2))

                         



-- * Householder bidiagonalization

{- G & VL Alg. 5.4.2 -}






-- * SVD

{- Golub & Van Loan, sec 8.6.2 (p 452 segg.)

SVD of A, Golub-Kahan method

* reduce A to upper bidiagonal form B (Alg. 5.4.2, Householder bidiagonalization)
* compute SVD of B (implicit-shift QR step applied to B^T B, Alg. 8.3.2)

-}








-- * Cholesky factorization

-- ** Cholesky–Banachiewicz algorithm

-- | Given a positive semidefinite matrix A, returns a lower-triangular matrix L such that L L^T = A
chol :: (Epsilon a, Real a, Floating a) => SpMatrix a -> SpMatrix a
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
                | i > 0 = sqrt $ aai - sum (fmap (**2) lrow)
                | otherwise = error "cholDiag : index must be nonnegative" where
        aai = aa@@(i,i)
        lrow = ifilterSV (\j _ -> j < i) (extractRow ll i) -- sub-diagonal elems of L

















-- * LU factorization
-- ** Doolittle algorithm
{- Doolittle algorithm for factoring A' = P A, where P is a permutation matrix such that A' has a nonzero as its (0, 0) entry -}

-- | Given a matrix A, returns a pair of matrices (L, U) where L is lower triangular and U is upper triangular such that L U = A
lu :: (Epsilon a, Fractional a, Real a) => SpMatrix a -> (SpMatrix a, SpMatrix a)
lu aa = (lf, ufin) where
  (ixf, lf, uf) = execState (modifyUntil q luUpd) luInit
  ufin = uUpdSparse (ixf, lf, uf) -- final U update
  q (i, _, _) = i == (nrows aa - 1)
  n = nrows aa
  luInit = (1, l0, u0) where
    l0 = insertCol (eye n) ((1/u00) .* extractSubCol aa 0 (1,n - 1)) 0  -- initial L
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
                       show (j ,j ),
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











-- * Incomplete LU

-- | used for Incomplete LU : remove entries in `m` corresponding to zero entries in `m2`


ilu0 :: (Epsilon a, Real a, Fractional a) => SpMatrix a -> (SpMatrix a, SpMatrix a)
ilu0 aa = (lh, uh) where
  (l, u) = lu aa
  lh = sparsifyLU l aa
  uh = sparsifyLU u aa
  sparsifyLU m m2 = ifilterSM f m where
    f i j _ = isJust (lookupSM m2 i j)






-- * Arnoldi iteration

-- | Given a matrix A and a positive integer `n`, this procedure finds the basis of an order `n` Krylov subspace (as the columns of matrix Q), along with an upper Hessenberg matrix H, such that A = Q^T H Q.
-- At the i`th iteration, it finds (i + 1) coefficients (the i`th column of the Hessenberg matrix H) and the (i + 1)`th Krylov vector.
arnoldi ::
  (Floating a, Eq a) => SpMatrix a -> Int -> (SpMatrix a, SpMatrix a)
arnoldi aa kn = (fromCols qvfin, hhfin)
  where
  (qvfin, hhfin, _) = execState (modifyUntil tf arnoldiStep) arnInit 
  tf (_, _, ii) = ii == kn -- termination criterion
  (m, n) = dim aa
  arnInit = (qv1, hh1, 1) where      
      q0 = normalize 2 $ onesSV n -- starting basis vector
      aq0 = aa #> q0              -- A q0
      h11 = q0 `dot` aq0          
      q1nn = (aq0 ^-^ (h11 .* q0))
      hh1 = fromListSM (m + 1, n) [(0, 0, h11), (1, 0, h21)] where        
        h21 = norm 2 q1nn
      q1 = normalize 2 q1nn       -- q1 `dot` q0 ~ 0
      qv1 = V.fromList [q0, q1]
  arnoldiStep (qv, hh, i) = (qv', hh', i + 1) where
    qi = V.last qv
    aqi = aa #> qi
    hhcoli = fmap (`dot` aqi) qv -- H_{1, i}, H_{2, i}, .. , H_{m + 1, i}
    zv = zeroSV m
    qipnn =
      aqi ^-^ (V.foldl' (^+^) zv (V.zipWith (.*) hhcoli qv)) -- unnormalized q_{i+1}
    qipnorm = singletonSV $ norm 2 qipnn      -- normalization factor H_{i+1, i}
    qip = normalize 2 qipnn              -- q_{i + 1}
    hh' = insertCol hh (concatSV (fromVector hhcoli) qipnorm) i -- update H
    qv' = V.snoc qv qip        -- append q_{i+1} to Krylov basis Q_i


-- -- -- why does q2 lose orthogonality ? because the unnormalized q vectors have decreasing norm and the Arnoldi iteration breaks down when this norm is ~ 0
-- n = 3
-- q0 = normalize 2 $ onesSV n

-- aq0 = aa #> q0
-- h00 = q0 `dot` aq0
-- q1nn = aq0 ^-^ (h00 .* q0)
-- q1 = normalize 2 q1nn

-- aq1 = aa #> q1
-- h01 = q0 `dot` aq1
-- h11 = q1 `dot` aq1
-- q2nn = aq1 ^-^ ((h01 .* q0) ^+^ (h11 .* q1))
-- q2 = normalize 2 q2nn

-- -- --



-- -- -- test data
-- (qv, hh) = arnoldi aa 3

-- -- -- columns of qv should be orthonormal
-- -- q1 = extractCol qv 1
-- -- q2 = extractCol qv 2


-- aa :: SpMatrix Double
-- aa = sparsifySM $ fromListDenseSM 3 [2, -1, 0, -1, 2, -1, 0, -1, 2]

-- x2, b2, x2i :: SpVector Double
-- x2 = mkSpVectorD 3 [3,2,3]
-- b2 = mkSpVectorD 3 [4,-2,4]

-- x2i = mkSpVectorD 3 [1,1,1]



-- * Preconditioning

-- | Partition a matrix into strictly subdiagonal, diagonal and strictly superdiagonal parts
diagPartitions :: SpMatrix a -> (SpMatrix a, SpMatrix a, SpMatrix a)
diagPartitions aa = (e,d,f) where
  e = extractSubDiag aa
  d = extractDiag aa
  f = extractSuperDiag aa


-- ** SSOR

-- | `mSsor aa omega` : if `omega = 1` it returns the symmetric Gauss-Seidel preconditioner
mSsor :: Fractional a => SpMatrix a -> a -> (SpMatrix a, SpMatrix a)
mSsor aa omega = (l, r) where
  (e, d, f) = diagPartitions aa
  n = nrows e
  l = (eye n ^-^ scale omega e) ## reciprocal d
  r = d ^-^ scale omega f 








-- * Linear solver, LU-based

-- | Direct solver based on a triangular factorization of the system matrix.
luSolve ::
  (Fractional a, Eq a, Epsilon a) => SpMatrix a -> SpMatrix a -> SpVector a -> SpVector a
luSolve ll uu b
  | isLowerTriSM ll && isUpperTriSM uu = lubwSolve uu (lufwSolve ll b)
  | otherwise = error "luSolve : factors must be triangular matrices" 

lufwSolve ll b = sparsifySV v where
  (v, _) = execState (modifyUntil q lStep) lInit where
  q (_, i) = i == dim b
  lStep (ww, i) = (wwi, i + 1) where
    lii = ll @@ (i, i)
    bi = b @@ i
    wi = (bi - r)/lii where
      r = extractSubRow ll i (0, i-1) `dot` takeSV i ww
    wwi = insertSpVector i wi ww
  lInit = (ww0, 1) where
    l00 = ll @@ (0, 0)
    b0 = b @@ 0
    w0 = b0 / l00
    ww0 = insertSpVector 0 w0 $ zeroSV (dim b)  

-- | NB in the computation of `xi` we must rebalance the subrow indices because `dropSV` does that too, in order to take the inner product with consistent index pairs
lubwSolve uu w = sparsifySV x where
  (x, _) = execState (modifyUntil q uStep) uInit
  q (_, i) = i == (- 1)
  uStep (xx, i) = (xxi, i - 1) where
    uii = uu @@ (i, i)
    wi = w @@ i
    xi = (wi - r) / uii where
        r = extractSubRow_RK uu i (i + 1, dim w - 1) `dot` dropSV (i + 1) xx
    xxi = insertSpVector i xi xx
  uInit = (xx0, i - 1) where
    i = dim w - 1
    u00 = uu @@ (i, i)
    w0 = w @@ i
    x0 = w0 / u00
    xx0 = insertSpVector i x0 $ zeroSV (dim w)






-- * Iterative linear solvers


-- ** GMRES

-- *** Left-preconditioning





-- ** CGNE

cgneStep :: SpMatrix Double -> CGNE -> CGNE
cgneStep aa (CGNE x r p) = CGNE x1 r1 p1 where
  alphai = r `dot` r / (p `dot` p)
  x1 = x ^+^ (alphai .* p)
  r1 = r ^-^ (alphai .* (aa #> p))
  beta = r1 `dot` r1 / (r `dot` r)
  p1 = transposeSM aa #> r ^+^ (beta .* p)

data CGNE =
  CGNE {_xCgne , _rCgne, _pCgne :: SpVector Double} deriving Eq
instance Show CGNE where
    show (CGNE x r p) = "x = " ++ show x ++ "\n" ++
                       "r = " ++ show r ++ "\n" ++
                       "p = " ++ show p ++ "\n"

cgne :: SpMatrix Double -> SpVector Double -> SpVector Double -> CGNE
cgne aa b x0 = execState (untilConverged _xCgne (cgneStep aa)) cgneInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = transposeSM aa #> r0
  cgneInit = CGNE x0 r0 p0


-- ** TFQMR

-- | one step of TFQMR
tfqmrStep :: SpMatrix Double -> SpVector Double -> TFQMR -> TFQMR
tfqmrStep aa r0hat (TFQMR x w u v d m tau theta eta rho alpha) =
  TFQMR x1 w1 u1 v1 d1 (m+1) tau1 theta1 eta1 rho1 alpha1
  where
  -- alpham = alpha
  w1 = w ^-^ (alpha .* (aa #> u))
  d1 = u ^+^ ((theta**2/alpha*eta) .* d)
  theta1 = norm2 w1 / tau
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

tfqmr :: SpMatrix Double -> SpVector Double -> SpVector Double -> TFQMR  
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
  tau0 = norm2 r0
  theta0 = 0
  eta0 = 0
  tfqmrInit = TFQMR x0 w0 u0 v0 d0 m tau0 theta0 eta0 rho0 alpha0

data TFQMR =
  TFQMR { _xTfq, _wTfq, _uTfq, _vTfq, _dTfq :: SpVector Double,
          _mTfq :: Int,
          _tauTfq, _thetaTfq, _etaTfq, _rhoTfq, _alphaTfq :: Double}
  deriving Eq
instance Show TFQMR where
    show (TFQMR x _ _ _ _ _ _ _ _ _ _) = "x = " ++ show x ++ "\n"


-- ** BCG

-- | one step of BCG
bcgStep :: SpMatrix Double -> BCG -> BCG
bcgStep aa (BCG x r rhat p phat) = BCG x1 r1 rhat1 p1 phat1 where
  aap = aa #> p
  alpha = (r `dot` rhat) / (aap `dot` phat)
  x1 = x ^+^ (alpha .* p)
  r1 = r ^-^ (alpha .* aap)
  rhat1 = rhat ^-^ (alpha .* (transposeSM aa #> phat))
  beta = (r1 `dot` rhat1) / (r `dot` rhat)
  p1 = r1 ^+^ (beta .* p)
  phat1 = rhat1 ^+^ (beta .* phat)

data BCG =
  BCG { _xBcg, _rBcg, _rHatBcg, _pBcg, _pHatBcg :: SpVector Double } deriving Eq

bcg :: SpMatrix Double -> SpVector Double -> SpVector Double -> BCG
bcg aa b x0 = execState (untilConverged _xBcg (bcgStep aa)) bcgInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  r0hat = r0
  p0 = r0
  p0hat = r0
  bcgInit = BCG x0 r0 r0hat p0 p0hat

instance Show BCG where
  show (BCG x r rhat p phat) = "x = " ++ show x ++ "\n" ++
                       "r = " ++ show r ++ "\n" ++
                       "r_hat = " ++ show rhat ++ "\n" ++
                       "p = " ++ show p ++ "\n" ++
                       "p_hat = " ++ show phat ++ "\n"


-- ** CGS

-- | one step of CGS
cgsStep :: SpMatrix Double -> SpVector Double -> CGS -> CGS
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

data CGS = CGS { _x, _r, _p, _u :: SpVector Double} deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
cgs ::
  SpMatrix Double ->
  SpVector Double ->
  SpVector Double ->
  SpVector Double ->
  CGS
cgs aa b x0 rhat =
  execState (untilConverged _x (cgsStep aa rhat)) cgsInit where
  r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
  p0 = r0
  u0 = r0
  cgsInit = CGS x0 r0 p0 u0


instance Show CGS where
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
bicgstabStep :: SpMatrix Double -> SpVector Double -> BICGSTAB -> BICGSTAB
bicgstabStep aa r0hat (BICGSTAB x r p) = BICGSTAB xj1 rj1 pj1 where
  aap = aa #> p
  alphaj = (r `dot` r0hat) / (aap `dot` r0hat)
  sj = r ^-^ (alphaj .* aap)
  aasj = aa #> sj
  omegaj = (aasj `dot` sj) / (aasj `dot` aasj)
  xj1 = x ^+^ (alphaj .* p) ^+^ (omegaj .* sj)    -- updated solution
  rj1 = sj ^-^ (omegaj .* aasj)
  betaj = (rj1 `dot` r0hat)/(r `dot` r0hat) * alphaj / omegaj
  pj1 = rj1 ^+^ (betaj .* (p ^-^ (omegaj .* aap)))

data BICGSTAB =
  BICGSTAB { _xBicgstab, _rBicgstab, _pBicgstab :: SpVector Double} deriving Eq

-- | iterate solver until convergence or until max # of iterations is reached
bicgstab
  :: SpMatrix Double
     -> SpVector Double
     -> SpVector Double
     -> SpVector Double
     -> BICGSTAB
bicgstab aa b x0 r0hat =
  execState (untilConverged _xBicgstab (bicgstabStep aa r0hat)) bicgsInit where
   r0 = b ^-^ (aa #> x0)    -- residual of initial guess solution
   p0 = r0
   bicgsInit = BICGSTAB x0 r0 p0

instance Show BICGSTAB where
  show (BICGSTAB x r p) = "x = " ++ show x ++ "\n" ++
                                "r = " ++ show r ++ "\n" ++
                                "p = " ++ show p ++ "\n"





-- * Moore-Penrose pseudoinverse
-- | Least-squares approximation of a rectangular system of equaitons. Uses <\\> for the linear solve
pinv :: SpMatrix Double -> SpVector Double -> SpVector Double
pinv aa b = aa #^# aa <\> atb where
  atb = transposeSM aa #> b







-- * Linear solver interface

data LinSolveMethod = CGNE_ | TFQMR_ | BCG_ | CGS_ | BICGSTAB_ deriving (Eq, Show) 

-- -- | Linear solve with _random_ starting vector
-- linSolveM ::
--   PrimMonad m =>
--     LinSolveMethod -> SpMatrix Double -> SpVector Double -> m (SpVector Double)
-- linSolveM method aa b = do
--   let (m, n) = dim aa
--       nb     = dim b
--   if n /= nb then error "linSolve : operand dimensions mismatch" else do
--     x0 <- randVec nb
--     case method of CGS_ -> return $ _xBicgstab (bicgstab aa b x0 x0)
--                    BICGSTAB_ -> return $ _x (cgs aa b x0 x0)

-- | Linear solve with _deterministic_ starting vector (every component at 0.1) 
linSolve ::
  LinSolveMethod -> SpMatrix Double -> SpVector Double -> SpVector Double
linSolve method aa b
  | n /= nb = error "linSolve : operand dimensions mismatch"
  | otherwise = solve aa b where
      solve aa' b' | isDiagonalSM aa' = reciprocal aa' #> b' -- diagonal solve is easy
                   | otherwise = solveWith aa' b' 
      solveWith aa' b' = case method of
        CGNE_ -> _xCgne (cgne aa' b' x0)
        TFQMR_ -> _xTfq (tfqmr aa' b' x0)
        BCG_ -> _xBcg (bcg aa' b' x0)
        BICGSTAB_ ->  _xBicgstab (bicgstab aa' b' x0 x0)
        CGS_ -> _x (cgs aa' b' x0 x0)
      x0 = mkSpVectorD n $ replicate n 0.1 
      (m, n) = dim aa
      nb     = dim b

-- | <\\> : linSolve using the BiCGSTAB method as default
(<\>) :: SpMatrix Double -> SpVector Double -> SpVector Double      
(<\>) = linSolve BICGSTAB_ 
  









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

-- | Keep a moving window buffer (length 2) of state `x` to assess convergence, stop when either a condition on that list is satisfied or when max # of iterations is reached  
loopUntilAcc :: Int -> ([t] -> Bool) -> (t -> t)  -> t -> t
loopUntilAcc nitermax q f x = go 0 [] x where
  go i ll xx | length ll < 2 = go (i + 1) (y : ll) y 
             | otherwise = if q ll || i == nitermax
                           then xx
                           else go (i + 1) (take 2 $ y:ll) y
                where y = f xx

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


-- helper functions for estimating convergence
meanl :: (Foldable t, Fractional a) => t a -> a
meanl xx = 1/fromIntegral (length xx) * sum xx

norm2l :: (Foldable t, Functor t, Floating a) => t a -> a
norm2l xx = sqrt $ sum (fmap (**2) xx)

diffSqL :: Floating a => [a] -> a
diffSqL xx = (x1 - x2)**2 where [x1, x2] = [head xx, xx!!1]







-- | iterate until convergence is verified or we run out of a fixed iteration budget
untilConverged :: MonadState a m => (a -> SpVector Double) -> (a -> a) -> m a
untilConverged fproj = modifyInspectN 200 (normDiffConverged fproj)

-- | convergence check (FIXME)
normDiffConverged :: (Foldable t, Functor t) =>
     (a -> SpVector Double) -> t a -> Bool
normDiffConverged fp xx = nearZero $ normSq (foldrMap fp (^-^) (zeroSV 0) xx)


  



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

  












-- * Random arrays

randArray :: PrimMonad m => Int -> Double -> Double -> m [Double]
randArray n mu sig = do
  g <- MWC.create
  replicateM n (MWC.normal mu sig g)
  



-- * Random matrices and vectors

-- |Dense SpMatrix
randMat :: PrimMonad m => Int -> m (SpMatrix Double)
randMat n = do
  g <- MWC.create
  aav <- replicateM (n^2) (MWC.normal 0 1 g)
  let ii_ = [0 .. n-1]
      (ix_,iy_) = unzip $ concatMap (zip ii_ . replicate n) ii_
  return $ fromListSM (n,n) $ zip3 ix_ iy_ aav

-- | Dense SpVector  
randVec :: PrimMonad m => Int -> m (SpVector Double)
randVec n = do
  g <- MWC.create
  bv <- replicateM n (MWC.normal 0 1 g)
  let ii_ = [0..n-1]
  return $ fromListSV n $ zip ii_ bv



-- | Sparse SpMatrix
randSpMat :: Int -> Int -> IO (SpMatrix Double)
randSpMat n nsp | nsp > n^2 = error "randSpMat : nsp must be < n^2 "
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  jj <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSM (n,n) $ zip3 ii jj aav

-- | Sparse SpVector
randSpVec :: Int -> Int -> IO (SpVector Double)
randSpVec n nsp | nsp > n = error "randSpVec : nsp must be < n"
                | otherwise = do
  g <- MWC.create
  aav <- replicateM nsp (MWC.normal 0 1 g)
  ii <- replicateM nsp (MWC.uniformR (0, n-1) g :: IO Int)
  return $ fromListSV n $ zip ii aav



