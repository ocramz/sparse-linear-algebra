{-# language TypeFamilies, FlexibleContexts #-}
{-# OPTIONS_GHC -Wno-missing-signatures -Wno-unused-local-binds -Wno-name-shadowing #-}
module Numeric.LinearAlgebra.EigenSolvers.Experimental where

import Control.Exception.Common
import Control.Iterative
import Data.Sparse.Common

import Control.Monad.Catch
import Control.Monad.State.Strict



-- | `eigRayleigh n mm` performs `n` iterations of the Rayleigh algorithm on matrix `mm` and returns the eigenpair closest to the initialization. It displays cubic-order convergence, but it also requires an educated guess on the initial eigenpair.

-- eigRayleigh nitermax debq prntf m = untilConvergedGM "eigRayleigh" config (const True) (rayStep m)
--   where
--     ii = eye (nrows m)
--     config = IterConf nitermax debq fst prntf
--     rayStep aa (b, mu) = do
--       nom <- (m ^-^ (mu `matScale` ii)) <\> b
--       let b' = normalize2' nom
--           mu' = (b' <.> (aa #> b')) / (b' <.> b')
--       return (b', mu')

      

-- | Golub-Kahan-Lanczos bidiagonalization (see "Restarted Lanczos Bidiagonalization for the SVD", SLEPc STR-8, http://slepc.upv.es/documentation/manual.htm )
gklBidiag aa q1nn | dim q1nn == n = return (pp, bb, qq)
                  | otherwise = throwM (MatVecSizeMismatchException "hhBidiag" (dim aa) (dim q1nn))
  where
  (m,n) = (nrows aa, ncols aa)
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
