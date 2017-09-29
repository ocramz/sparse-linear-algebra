{-# language TypeFamilies, FlexibleContexts #-}
module Numeric.LinearAlgebra.LinearSolvers.Experimental where

import Control.Exception.Common
import Control.Iterative
import Data.Sparse.Common

import Control.Monad.Catch




-- ** TFQMR

tfqmrInit aa b x0 = TFQMR x0 w0 u0 v0 d0 m tau0 theta0 eta0 rho0 alpha0 where
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
