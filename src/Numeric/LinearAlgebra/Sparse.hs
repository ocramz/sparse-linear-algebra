{-# LANGUAGE FlexibleContexts, TypeFamilies, MultiParamTypeClasses, FlexibleInstances  #-}
-- {-# OPTIONS_GHC -O2 -rtsopts -with-rtsopts=-K32m -prof#-}
module Numeric.LinearAlgebra.Sparse
       (
         -- * Matrix factorizations
         qr, lu,
         -- * Incomplete LU
         ilu0,
         -- * Condition number
         conditionNumberSM,
         -- * Householder reflection
         hhMat, hhRefl,
         -- * Givens' rotation
         givens,
         -- * Eigensolvers
         eigsQR, eigRayleigh,
         -- * Linear solvers
         linSolve, LinSolveMethod, (<\>),
         -- ** Methods
         cgne, tfqmr, bicgstab, cgs, bcg,
         _xCgne, _xTfq, _xBicgstab, _x, _xBcg,
         cgsStep, bicgstabStep,
         CGNE, TFQMR, BICGSTAB, CGS, BCG,
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
         modifyInspectN, runAppendN',
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





  
-- * Sparsify : remove almost-0 elements (|x| < eps)
-- | Sparsify an SpVector
sparsifySV :: SpVector Double -> SpVector Double
sparsifySV (SV d im) = SV d $ IM.filter (\x -> abs x >= eps) im





-- * Matrix condition number

-- |uses the R matrix from the QR factorization
conditionNumberSM :: SpMatrix Double -> Double
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
hhRefl :: SpVector Double -> SpMatrix Double
hhRefl = hhMat 2.0











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
-}

givens :: SpMatrix Double -> IxRow -> IxCol -> SpMatrix Double
givens mm i j 
  | isValidIxSM mm (i,j) && isSquareSM mm =
       sparsifySM $ fromListSM' [(i,i,c),(j,j,c),(j,i,-s),(i,j,s)] (eye (nrows mm))
  | otherwise = error "givens : indices out of bounds"      
  where
    (c, s, _) = givensCoef a b
    i' = head $ fromMaybe (error $ "givens: no compatible rows for entry " ++ show (i,j)) (candidateRows (immSM mm) i j)
    a = mm @@ (i', j)
    b = mm @@ (i, j)   -- element to zero out

-- |Is the `k`th the first nonzero column in the row?
firstNonZeroColumn :: IM.IntMap a -> IxRow -> Bool
firstNonZeroColumn mm k = isJust (IM.lookup k mm) &&
                          isNothing (IM.lookupLT k mm)

-- |Returns a set of rows {k} that satisfy QR.C1
candidateRows :: IM.IntMap (IM.IntMap a) -> IxRow -> IxCol -> Maybe [IM.Key]
candidateRows mm i j | IM.null u = Nothing
                     | otherwise = Just (IM.keys u) where
  u = IM.filterWithKey (\irow row -> irow /= i &&
                                     firstNonZeroColumn row j) mm






-- * QR decomposition


-- | Applies Givens rotation iteratively to zero out sub-diagonal elements
qr :: SpMatrix Double -> (SpMatrix Double, SpMatrix Double)
qr mm = (transposeSM qmatt, rmat)  where
  qmatt = F.foldl' (#~#) ee $ gmats mm -- Q^T = (G_n * G_n-1 ... * G_1)
  rmat = qmatt #~# mm                  -- R = Q^T A
  ee = eye (nrows mm)
      
-- | Givens matrices in order [G1, G2, .. , G_N ]
gmats :: SpMatrix Double -> [SpMatrix Double]
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
eigsQR :: Int -> SpMatrix Double -> SpVector Double
eigsQR nitermax m = extractDiagDense $ execState (convergtest eigsStep) m where
  eigsStep m = r #~# q where (q, r) = qr m
  convergtest g = modifyInspectN nitermax f g where
    f [m1, m2] = let dm1 = extractDiagDense m1
                     dm2 = extractDiagDense m2
                 in norm2 (dm1 ^-^ dm2) <= eps






-- ** Rayleigh iteration

-- | `eigsRayleigh n mm` performs `n` iterations of the Rayleigh algorithm on matrix `mm` and returns the eigenpair closest to the initialization. It displays cubic-order convergence, but it also requires an educated guess on the initial eigenpair
eigRayleigh :: Int                -- max # iterations
     -> SpMatrix Double           -- matrix
     -> (SpVector Double, Double) -- initial guess of (eigenvector, eigenvalue)
     -> (SpVector Double, Double) -- final estimate of (eigenvector, eigenvalue)
eigRayleigh nitermax m = execState (convergtest (rayleighStep m)) where
  convergtest g = modifyInspectN nitermax f g where
    f [(b1, _), (b2, _)] = norm2 (b2 ^-^ b1) <= eps 
  rayleighStep aa (b, mu) = (b', mu') where
      ii = eye (nrows aa)
      nom = (aa ^-^ (mu `matScale` ii)) <\> b
      b' = normalize 2 nom
      mu' = b' `dot` (aa #> b') / (b' `dot` b')




-- * Householder vector 

-- (Golub & Van Loan, Alg. 5.1.1, function `house`)
hhV :: SpVector Double -> (SpVector Double, Double)
hhV x = (v, beta) where
  n = dim x
  tx = tailSV x
  sigma = tx `dot` tx
  vtemp = singletonSV 1 `concatSV` tx
  (v, beta) | sigma <= eps = (vtemp, 0)
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















-- * LU factorization
-- ** Doolittle algorithm
{- Doolittle algorithm for factoring A' = P A, where P is a permutation matrix such that A' has a nonzero as its (0, 0) entry -}

-- | LU factors
lu :: SpMatrix Double -> (SpMatrix Double, SpMatrix Double)
lu aa = (lfin, ufin) where
  (ixf,lf,uf) = execState (modifyUntil q (luUpd aa)) (luInit aa)
  lfin = lf
  ufin = uUpd aa (ixf, lf, uf)
  q (i, _, _) = i == (nrows aa - 1)

-- | First iteration of LU
luInit ::
  (Num t, Fractional a) => SpMatrix a -> (t, SpMatrix a, SpMatrix a)
luInit aa = (1, l0, u0) where
  n = nrows aa
  l0 = insertCol (eye n) ((1/u00) .* extractSubCol aa 0 (1,n - 1)) 0  -- initial L
  u0 = insertRow (zeroSM n n) (extractRow aa 0) 0               -- initial U
  u00 = u0 @@ (0,0)  -- make sure this is non-zero by applying permutation

-- | LU update step
luUpd :: SpMatrix Double
     -> (Int, SpMatrix Double, SpMatrix Double)
     -> (Int, SpMatrix Double, SpMatrix Double)
luUpd aa (i, l, u) = (i', l', u') where
  n = nrows aa  
  u' = uUpdSparse aa (i, l, u)  -- update U
  l' = lUpdSparse aa (i, l, u') -- update L
  i' = i + 1     -- increment i


uUpd' ::
  Num a =>
  ([(Int, a)] -> [(Int, a)]) ->
  SpMatrix a ->
  (Rows, SpMatrix a, SpMatrix a) ->
  SpMatrix a
uUpd' ff amat (ix, lmat, umat) = insertRow umat uv ix where
  n = nrows amat
  colsix = [ix .. n - 1]
  us = ff $ zip colsix $ map (solveForUij amat lmat umat ix) colsix
  uv = fromListSV n us

uUpd :: Num a => SpMatrix a -> (Rows, SpMatrix a, SpMatrix a) -> SpMatrix a
uUpd = uUpd' id

-- update U while sparsifying
uUpdSparse ::
  SpMatrix Double -> (Rows, SpMatrix Double, SpMatrix Double) -> SpMatrix Double
uUpdSparse = uUpd' (filter (isNz . snd))




-- solve for element Uij
solveForUij ::
  Num a => SpMatrix a -> SpMatrix a -> SpMatrix a -> IxRow -> IxCol -> a
solveForUij amat lmat umat i j = a - p where
  a = amat @@! (i, j)
  p = contractSub lmat umat i j (i - 1)


-- solve for element Lij
solveForLij ::
  SpMatrix Double -> SpMatrix Double -> SpMatrix Double -> IxRow -> IxCol -> Double
solveForLij amat lmat umat i j
  | isNz ujj = (a - p)/ujj
  | otherwise =
     error $ unwords ["solveForLij : U",
                      show (j ,j ),
                      "is close to 0. Permute rows in order to have a nonzero diagonal of U"]
  where
   a = amat @@! (i, j)
   ujj = umat @@! (j , j)   -- NB this must be /= 0
   p = contractSub lmat umat i j (i - 1)




lUpd' :: ([(Rows, Double)] -> [(Int, Double)])
     -> SpMatrix Double
     -> (Rows, SpMatrix Double, SpMatrix Double)
     -> SpMatrix Double
lUpd' ff amat (ix, lmat, umat) = insertCol lmat lv ix where
  n = nrows amat
  rowsix = [ix + 1 .. n - 1]
  ls = ff $ zip rowsix $ map (\i -> solveForLij amat lmat umat i ix) rowsix
  lv = fromListSV n ls

lUpd :: SpMatrix Double -> (Rows, SpMatrix Double, SpMatrix Double) -> SpMatrix Double
lUpd = lUpd' id

lUpdSparse ::
  SpMatrix Double -> (Rows, SpMatrix Double, SpMatrix Double) -> SpMatrix Double
lUpdSparse = lUpd' (filter (isNz . snd))

















-- Produces the permutation matrix necessary to have a nonzero in position (iref, jref). This is used in the LU factorization
permutAA :: Num b => IxRow -> IxCol -> SpMatrix a -> Maybe (SpMatrix b)
permutAA iref jref (SM (nro,_) mm) 
  | isJust (lookupIM2 iref jref mm) = Nothing -- eye nro
  | otherwise = Just $ permutationSM nro [head u] where
      u = IM.keys (ifilterIM2 ff mm)
      ff i j _ = i /= iref &&
                 j == jref











-- * Incomplete LU

-- | used for Incomplete LU : remove entries in `m` corresponding to zero entries in `m2`



ilu0 aa = (lh, uh) where
  (l, u) = lu aa
  lh = sparsifyLU l aa
  uh = sparsifyLU u aa
  sparsifyLU m m2 = SM (dim m) $ ifilterIM2 f (dat m) where
    f i j _ = isJust (lookupSM m2 i j)










-- * Preconditioning

-- | Partition a matrix into strictly subdiagonal, diagonal and strictly superdiagonal parts
diagPartitions :: SpMatrix a -> (SpMatrix a, SpMatrix a, SpMatrix a)
diagPartitions aa = (e,d,f) where
  e = extractSubDiag aa
  d = extractDiag aa
  f = extractSuperDiag aa


-- ** SSOR

-- | `mSsor aa omega` : if `omega = 1` it returns the diagonal of `aa`, 
mSsor :: Fractional a => SpMatrix a -> a -> SpMatrix a
mSsor aa omega = l ## r where
  (e, d, f) = diagPartitions aa
  n = nrows e
  l = d ^-^ scale omega e
  r = eye n ^-^ scale omega (reciprocal d ## f)








-- * Iterative linear solvers




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
        CGS_ ->  _xBicgstab (bicgstab aa' b' x0 x0)
        BICGSTAB_ -> _x (cgs aa' b' x0 x0)
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
normDiffConverged fp xx = normSq (foldrMap fp (^-^) (zeroSV 0) xx) <= eps


  



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



